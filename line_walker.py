#!/usr/bin/env python3
"""
line_walker.py — Iterative LINE reconstruction from a SINE tail seed.

Given a SINE tail sequence and a genome, walks along genomic LINE copies
to reconstruct the full LINE element step by step.

Pipeline per step:
    (optional) extended-flank fast path:
             prev-step hits → bedtools getfasta (--extended-flank bp)
             → MAFFT align → vsearch cluster → consensus  [early exit if OK]
    fallback: seed → sear -s 15 -k  (stop at 15 hits, reuse genome splits)
       → sort by bitscore, take top 10
       → extract 150 bp directional flank (strand-aware)
       → MAFFT align flanks
       → vsearch cluster at 80 %
       → ≥3 members per cluster → majority-rule consensus → extend
       → 1 cluster: continue;  N clusters: fork branches

Stop conditions:
  - fewer than --branch-min hits / flanks / cluster members
  - consensus shorter than 30 bp
  - --steps limit reached

Dependencies: sear, samtools, bedtools, mafft, vsearch
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

# ═══════════════════════════════════════════════════════════════════════
# FASTA I/O
# ═══════════════════════════════════════════════════════════════════════

def read_fasta(path):
    """Read single-sequence FASTA → (name, sequence)."""
    name, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    break
                name = line[1:].split()[0]
            elif name:
                parts.append(line.upper())
    return name, ''.join(parts)


def read_multi_fasta(path):
    """Read multi-sequence FASTA → [(name, seq), ...]."""
    entries, name, parts = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    entries.append((name, ''.join(parts).upper()))
                name = line[1:].split()[0]
                parts = []
            elif name is not None:
                parts.append(line)
    if name is not None:
        entries.append((name, ''.join(parts).upper()))
    return entries


def write_fasta(path, name, seq, wrap=80):
    with open(path, 'w') as fh:
        fh.write(f'>{name}\n')
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + '\n')


def write_multi_fasta(path, entries, wrap=80):
    with open(path, 'w') as fh:
        for name, seq in entries:
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + '\n')


# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════

def run(cmd, cwd=None):
    """Run shell command; raise on failure."""
    r = subprocess.run(cmd, shell=True, cwd=cwd,
                       executable='/bin/bash',
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       text=True)
    if r.returncode != 0:
        sys.stderr.write(f"FAIL: {cmd}\n{r.stderr}\n")
        raise RuntimeError(cmd)
    return r.stdout


def load_chrom_sizes(fai_path):
    """Parse .fai → {chrom: length}."""
    sizes = {}
    with open(fai_path) as fh:
        for line in fh:
            p = line.split('\t')
            sizes[p[0]] = int(p[1])
    return sizes


def majority_consensus(aligned_seqs, min_depth=2):
    """Majority-rule consensus from aligned sequence strings.

    Positions where gaps outnumber bases are skipped (natural trimming).
    """
    if not aligned_seqs:
        return ''
    maxlen = max(len(s) for s in aligned_seqs)
    cons = []
    for i in range(maxlen):
        counts = {}
        gaps = 0
        for s in aligned_seqs:
            b = s[i].upper() if i < len(s) else '-'
            if b in 'ACGT':
                counts[b] = counts.get(b, 0) + 1
            else:
                gaps += 1
        total = sum(counts.values())
        if total >= min_depth and total > gaps:
            cons.append(max(counts, key=counts.get))
    return ''.join(cons)


# ═══════════════════════════════════════════════════════════════════════
# sear integration
# ═══════════════════════════════════════════════════════════════════════

def setup_sear_workdir(outdir, genome_abs):
    """Create shared sear working directory with genome symlink.

    Genome splits (*.2k.part_*.bnk) are created once here and reused
    across all walking steps thanks to sear -k.
    """
    wd = os.path.join(outdir, '_sear_work')
    os.makedirs(wd, exist_ok=True)
    link     = os.path.join(wd, 'genome.fa')
    fai_link = os.path.join(wd, 'genome.fa.fai')
    if not os.path.exists(link):
        os.symlink(genome_abs, link)
    if not os.path.exists(fai_link) and os.path.exists(genome_abs + '.fai'):
        os.symlink(genome_abs + '.fai', fai_link)
    return wd


def run_sear(sear_wd, query_fa, search_n, threads):
    """Run sear in the shared workdir.  Returns path to output BED."""
    qname = os.path.basename(query_fa)
    q_dest = os.path.join(sear_wd, qname)
    shutil.copy(query_fa, q_dest)

    cmd = f"THREADS={threads} sear -s {search_n} -k {qname} genome.fa"
    run(cmd, cwd=sear_wd)

    stem = os.path.splitext(qname)[0]
    bed  = os.path.join(sear_wd, f'gen-{stem}.bed')
    bnk  = os.path.join(sear_wd, f'gen-{stem}.bnk')

    # .bnk not needed (we extract our own flanks)
    if os.path.exists(bnk):
        os.remove(bnk)
    os.remove(q_dest)

    return bed if os.path.exists(bed) else None


def parse_bed7(path):
    """Parse 7-column sear BED → list of dicts, sorted by bitscore desc."""
    hits = []
    with open(path) as fh:
        for line in fh:
            p = line.strip().split('\t')
            if len(p) < 7:
                continue
            strand = p[5].strip().replace(',', '')
            if strand not in ('+', '-'):
                strand = '+'
            hits.append(dict(
                chrom=p[0], start=int(p[1]), end=int(p[2]),
                homology=float(p[3]), length=int(p[4]),
                strand=strand, bitscore=float(p[6])
            ))
    hits.sort(key=lambda h: h['bitscore'], reverse=True)
    return hits


def write_bed7(hits, out_path):
    """Write list of hit dicts as 7-column BED."""
    with open(out_path, 'w') as fh:
        for h in hits:
            fh.write(
                f"{h['chrom']}\t{h['start']}\t{h['end']}\t"
                f"{h['homology']}\t{h['length']}\t{h['strand']}\t"
                f"{h['bitscore']}\n"
            )


def _interval_distance(a_start, a_end, b_start, b_end):
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0


def filter_hits_by_previous_loci(all_hits, prev_hits, max_jump):
    """Keep hits close to previous-step loci on same chrom/strand.

    This preserves branch continuity and prevents branch jumping to unrelated
    genomic copies when a fresh sear search is required.
    """
    if not prev_hits:
        return all_hits

    kept = []
    for h in all_hits:
        for p in prev_hits:
            if h['chrom'] != p['chrom'] or h['strand'] != p['strand']:
                continue
            if _interval_distance(h['start'], h['end'],
                                  p['start'], p['end']) <= max_jump:
                kept.append(h)
                break
    return kept


# ═══════════════════════════════════════════════════════════════════════
# Flank coordinate logic
# ═══════════════════════════════════════════════════════════════════════

def write_flank_bed(hits, direction, flank_size, csizes, out_path):
    """Write flank-only BED.

    direction='downstream'  →  3′ end of the hit (on its strand)
    direction='upstream'    →  5′ end of the hit (on its strand)
    direction='right'       →  higher genomic coordinates (absolute)
    direction='left'        →  lower genomic coordinates (absolute)

    Returns number of valid flanks written.
    """
    n = 0
    with open(out_path, 'w') as fh:
        for i, h in enumerate(hits):
            clen = csizes.get(h['chrom'], 10**9)
            if direction == 'right':
                fs, fe = h['end'], min(h['end'] + flank_size, clen)
            elif direction == 'left':
                fs, fe = max(0, h['start'] - flank_size), h['start']
            elif direction == 'downstream':
                if h['strand'] == '+':
                    fs, fe = h['end'], min(h['end'] + flank_size, clen)
                else:
                    fs, fe = max(0, h['start'] - flank_size), h['start']
            else:  # upstream
                if h['strand'] == '+':
                    fs, fe = max(0, h['start'] - flank_size), h['start']
                else:
                    fs, fe = h['end'], min(h['end'] + flank_size, clen)
            if fe - fs < 20:
                continue
            fh.write(f"{h['chrom']}\t{fs}\t{fe}\thit{i}\t0\t{h['strand']}\n")
            n += 1
    return n


# ═══════════════════════════════════════════════════════════════════════
# One walking step
# ═══════════════════════════════════════════════════════════════════════

def _dump_stats(stats, step_dir):
    with open(os.path.join(step_dir, 'stats.json'), 'w') as fh:
        json.dump(stats, fh, indent=2)


def _try_extended_flanks(step_num, tag, prev_hits_bed, genome, sd, args,
                         csizes):
    """Attempt extended flank extraction from previous step's sear hits.

    Uses the same genomic hit locations discovered in the prior step but
    extracts a larger window (args.extended_flank bp) in the walk direction,
    then runs the full align → cluster → consensus pipeline.

    Returns (extensions, hits_bed7) or ([], None) if failed.
    """
    try:
        all_hits = parse_bed7(prev_hits_bed)
    except Exception:
        return [], None

    top = all_hits[:args.top_hits]
    if len(top) < args.branch_min:
        return [], None

    ext_sd = os.path.join(sd, 'ext')
    os.makedirs(ext_sd, exist_ok=True)

    ext_flank_bed = os.path.join(ext_sd, 'flanks.bed')
    nf = write_flank_bed(top, args.direction, args.extended_flank, csizes,
                         ext_flank_bed)
    if nf < args.branch_min:
        print(f"  [{tag}] extended: only {nf} usable flanks — skip")
        return [], None

    ext_flanks_fa = os.path.join(ext_sd, 'flanks.fa')
    try:
        run(f"bedtools getfasta -s -nameOnly -fi {genome}"
            f" -bed {ext_flank_bed} > {ext_flanks_fa}", cwd=ext_sd)
    except RuntimeError:
        return [], None

    seqs = read_multi_fasta(ext_flanks_fa)
    if len(seqs) < args.branch_min:
        return [], None

    ext_aligned_fa = os.path.join(ext_sd, 'aligned.fa')
    print(f"  [{tag}] extended: MAFFT aligning {len(seqs)} sequences "
          f"({args.extended_flank} bp flanks) ...")
    try:
        run(f"mafft --thread {args.threads} --threadtb {args.threads} "
            f"--localpair --maxiterate 1000 --ep 0.123 --nuc --reorder "
            f"--quiet {ext_flanks_fa} > {ext_aligned_fa}", cwd=ext_sd)
    except RuntimeError:
        return [], None

    clu_dir = os.path.join(ext_sd, 'clusters')
    os.makedirs(clu_dir, exist_ok=True)
    ctr = os.path.join(clu_dir, 'centroids.fa')
    uc = os.path.join(clu_dir, 'clusters.uc')
    pfx = os.path.join(clu_dir, 'cluster_')
    try:
        run(f"vsearch --cluster_fast {ext_flanks_fa} --id {args.cluster_id} "
            f"--centroids {ctr} --uc {uc} --clusters {pfx} --quiet",
            cwd=ext_sd)
    except RuntimeError:
        return [], None

    cluster_files = sorted(
        f for f in os.listdir(clu_dir)
        if f.startswith('cluster_')
        and f not in ('centroids.fa', 'clusters.uc')
        and not f.endswith('.uc') and not f.endswith('.fa')
    )
    qualifying = []
    for cf in cluster_files:
        members = read_multi_fasta(os.path.join(clu_dir, cf))
        if len(members) >= args.branch_min:
            qualifying.append(members)

    if not qualifying:
        print(f"  [{tag}] extended: no clusters ≥{args.branch_min} — stop")
        return [], None

    ext_hits_bed = os.path.join(ext_sd, 'extended_hits.bed')
    write_bed7(top, ext_hits_bed)

    extensions = []
    for ci, members in enumerate(qualifying):
        label = chr(65 + ci)
        cfa = os.path.join(ext_sd, f'cluster_{label}.fa')
        calign = os.path.join(ext_sd, f'cluster_{label}_aligned.fa')
        ccons = os.path.join(ext_sd, f'consensus_{label}.fa')

        write_multi_fasta(cfa, members)

        if len(qualifying) == 1 and os.path.exists(ext_aligned_fa):
            aligned = read_multi_fasta(ext_aligned_fa)
        elif len(members) > 1:
            try:
                run(f"mafft --thread {args.threads} --localpair "
                    f"--maxiterate 1000 --ep 0.123 --nuc --reorder --quiet "
                    f"{cfa} > {calign}", cwd=ext_sd)
            except RuntimeError:
                continue
            aligned = read_multi_fasta(calign)
        else:
            aligned = members

        cons = majority_consensus([s for _, s in aligned])
        if len(cons) < 30:
            print(f"  [{tag}] extended cluster {label}: consensus only "
                  f"{len(cons)} bp — skip")
            continue

        write_fasta(ccons, f'step{step_num}_{label}_ext', cons)
        extensions.append((cons, label))
        print(f"  [{tag}] extended cluster {label}: {len(members)} seqs → "
              f"{len(cons)} bp consensus")

    return extensions, ext_hits_bed


def walk_step(step_num, query_fa, genome, outdir, sear_wd, args, csizes,
              branch='', prev_hits_bed=None):
    """Execute one walking step.

    If args.try_extended_first is True and prev_hits_bed is available, first
    attempts to extract extended flanks (args.extended_flank bp) from the
    previous step's sear hits before falling back to a full genome search.

    Returns (extensions, next_hits_bed) where extensions is a list of
    (consensus_seq, cluster_label) for next steps and next_hits_bed is the
    path to the sear BED file to pass to the following step (None when the
    extended-extraction fast path was taken). Empty extensions list means
    this branch stops.
    """
    tag = f"step_{step_num:03d}" + (f"_branch_{branch}" if branch else "")
    sd  = os.path.join(outdir, tag)
    os.makedirs(sd, exist_ok=True)
    shutil.copy(query_fa, os.path.join(sd, 'query.fa'))

    stats = dict(step=step_num, branch=branch, tag=tag)

    if (args.try_extended_first
            and prev_hits_bed is not None
            and os.path.exists(prev_hits_bed)):
        print(f"  [{tag}] trying extended flanks "
              f"({args.extended_flank} bp, direction={args.direction}) ...")
        ext_exts, ext_hits_bed = _try_extended_flanks(
            step_num, tag, prev_hits_bed, genome, sd, args, csizes
        )
        if ext_exts:
            stats['status'] = 'extended_from_prev_hits'
            stats['extended_flank_size'] = args.extended_flank
            stats['extensions'] = [
                dict(branch=l, length=len(s)) for s, l in ext_exts
            ]
            _dump_stats(stats, sd)
            return ext_exts, ext_hits_bed
        print(f"  [{tag}] extended extraction failed — falling back to sear")

    # 1 ── sear search ─────────────────────────────────────────────────
    print(f"  [{tag}] sear -s {args.search_hits} ...")
    bed_src = run_sear(sear_wd, query_fa, args.search_hits, args.threads)
    if bed_src is None or os.path.getsize(bed_src) == 0:
        print(f"  [{tag}] no hits")
        stats['status'] = 'no_hits'
        _dump_stats(stats, sd)
        return [], None

    bed_dst = os.path.join(sd, 'sear_hits.bed')
    shutil.move(bed_src, bed_dst)

    # 2 ── select top hits by bitscore ─────────────────────────────────
    all_hits = parse_bed7(bed_dst)
    if prev_hits_bed is not None and os.path.exists(prev_hits_bed):
        prev_hits = parse_bed7(prev_hits_bed)
        filtered = filter_hits_by_previous_loci(all_hits, prev_hits,
                                                args.max_jump)
        if filtered:
            stats['hits_before_continuity_filter'] = len(all_hits)
            stats['hits_after_continuity_filter'] = len(filtered)
            all_hits = filtered
        else:
            stats['status'] = 'no_continuous_hits'
            stats['hits_before_continuity_filter'] = len(all_hits)
            stats['hits_after_continuity_filter'] = 0
            print(f"  [{tag}] continuity filter kept 0 hits — stop")
            _dump_stats(stats, sd)
            return [], bed_dst

    top = all_hits[:args.top_hits]
    stats['hits_total'] = len(all_hits)
    stats['hits_used']  = len(top)
    print(f"  [{tag}] {len(all_hits)} hits → top {len(top)} by bitscore")

    if len(top) < args.branch_min:
        stats['status'] = 'too_few_hits'
        _dump_stats(stats, sd)
        return [], bed_dst

    # 3 ── extract directional flanks ──────────────────────────────────
    flank_bed = os.path.join(sd, 'flanks.bed')
    nf = write_flank_bed(top, args.direction, args.flank, csizes, flank_bed)
    if nf < args.branch_min:
        print(f"  [{tag}] only {nf} usable flanks (need {args.branch_min})")
        stats['status'] = 'too_few_flanks'
        _dump_stats(stats, sd)
        return [], bed_dst

    flanks_fa = os.path.join(sd, 'flanks.fa')
    run(f"bedtools getfasta -s -nameOnly -fi {genome} -bed {flank_bed}"
        f" > {flanks_fa}", cwd=sd)

    seqs = read_multi_fasta(flanks_fa)
    if len(seqs) < args.branch_min:
        stats['status'] = 'too_few_flanks'
        _dump_stats(stats, sd)
        return [], bed_dst
    stats['flanks'] = len(seqs)

    # 4 ── MAFFT align flanks ──────────────────────────────────────────
    aligned_fa = os.path.join(sd, 'aligned.fa')
    print(f"  [{tag}] MAFFT aligning {len(seqs)} flanks ...")
    run(f"mafft --thread {args.threads} --threadtb {args.threads} "
        f"--localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet "
        f"{flanks_fa} > {aligned_fa}", cwd=sd)

    # 5 ── vsearch cluster ─────────────────────────────────────────────
    clu_dir = os.path.join(sd, 'clusters')
    os.makedirs(clu_dir, exist_ok=True)
    ctr = os.path.join(clu_dir, 'centroids.fa')
    uc  = os.path.join(clu_dir, 'clusters.uc')
    pfx = os.path.join(clu_dir, 'cluster_')

    print(f"  [{tag}] vsearch cluster (id={args.cluster_id}) ...")
    run(f"vsearch --cluster_fast {flanks_fa} --id {args.cluster_id} "
        f"--centroids {ctr} --uc {uc} --clusters {pfx} --quiet", cwd=sd)

    # vsearch --clusters creates files: cluster_0, cluster_1, …
    cluster_files = sorted(
        f for f in os.listdir(clu_dir)
        if f.startswith('cluster_')
        and f not in ('centroids.fa', 'clusters.uc')
        and not f.endswith('.uc') and not f.endswith('.fa')
    )
    qualifying = []
    for cf in cluster_files:
        members = read_multi_fasta(os.path.join(clu_dir, cf))
        if len(members) >= args.branch_min:
            qualifying.append(members)

    stats['clusters_total']      = len(cluster_files)
    stats['clusters_qualifying'] = len(qualifying)
    stats['cluster_sizes'] = [
        len(read_multi_fasta(os.path.join(clu_dir, cf)))
        for cf in cluster_files
    ]
    print(f"  [{tag}] {len(cluster_files)} clusters, "
          f"{len(qualifying)} with ≥{args.branch_min} members")

    # 6 ── build consensus ─────────────────────────────────────────────
    if not qualifying:
        print(f"  [{tag}] no clusters ≥{args.branch_min} — stop")
        stats['status'] = 'no_qualifying_cluster'
        _dump_stats(stats, sd)
        return [], bed_dst

    extensions = []
    for ci, members in enumerate(qualifying):
        label  = chr(65 + ci)                    # A, B, C, …
        cfa    = os.path.join(sd, f'cluster_{label}.fa')
        calign = os.path.join(sd, f'cluster_{label}_aligned.fa')
        ccons  = os.path.join(sd, f'consensus_{label}.fa')

        write_multi_fasta(cfa, members)

        if len(qualifying) == 1 and os.path.exists(aligned_fa):
            # Reuse the full MAFFT alignment already computed in step 4
            aligned = read_multi_fasta(aligned_fa)
        elif len(members) > 1:
            run(f"mafft --thread {args.threads} --localpair "
                f"--maxiterate 1000 --ep 0.123 --nuc --reorder --quiet "
                f"{cfa} > {calign}", cwd=sd)
            aligned = read_multi_fasta(calign)
        else:
            aligned = members

        cons = majority_consensus([s for _, s in aligned])

        if len(cons) < 30:
            print(f"  [{tag}] cluster {label}: consensus only "
                  f"{len(cons)} bp — skip")
            continue

        write_fasta(ccons, f'step{step_num}_{label}', cons)
        extensions.append((cons, label))
        print(f"  [{tag}] cluster {label}: {len(members)} seqs → "
              f"{len(cons)} bp consensus")

    if extensions:
        stats['status'] = 'extended'
        stats['extensions'] = [
            dict(branch=l, length=len(s)) for s, l in extensions
        ]
    else:
        stats['status'] = 'no_usable_consensus'
    _dump_stats(stats, sd)
    return extensions, bed_dst


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description='LINE_walker — reconstruct full LINE from SINE tail seed')

    ap.add_argument('-s', '--seed', required=True,
                    help='Seed FASTA (SINE tail or LINE fragment)')
    ap.add_argument('-g', '--genome', required=True,
                    help='Genome FASTA (will be indexed if needed)')
    ap.add_argument('-o', '--outdir', required=True,
                    help='Output directory')
    ap.add_argument('-d', '--direction', default='upstream',
                    choices=('upstream', 'downstream', 'left', 'right'),
                    help='Walk direction: upstream/downstream (strand-aware) '
                        'or left/right (absolute genomic coordinates). '
                        '(default: upstream)')

    g = ap.add_argument_group('search tuning')
    g.add_argument('--steps',       type=int,   default=30,
                   help='Max walking steps (default: 30)')
    g.add_argument('--search-hits', type=int,   default=50,
                   help='sear: stop after N hits (default: 50)')
    g.add_argument('--top-hits',    type=int,   default=20,
                   help='Keep N best hits by bitscore (default: 20)')
    g.add_argument('--flank',       type=int,   default=150,
                   help='Flank extraction size in bp (default: 150)')
    g.add_argument('--seed-window', type=int,   default=200,
                   help='Use last N bp of accumulated seq as query '
                        '(default: 200)')
    g.add_argument('--max-jump', type=int, default=5000,
                   help='Max genomic distance (bp) from previous-step hits '
                        'allowed when continuity filtering fresh sear hits '
                        '(default: 5000)')

    g = ap.add_argument_group('extended-flank optimization')
    g.add_argument('--extended-flank', type=int, default=500,
                   help='Flank size in bp for extended extraction from '
                        'previous step\'s hits before falling back to '
                        'genome search (default: 500)')
    g.add_argument('--try-extended-first', dest='try_extended_first',
                   action='store_true', default=True,
                   help='Attempt extended flank extraction from previous '
                        'step\'s hits before running genome search '
                        '(default: enabled)')
    g.add_argument('--no-try-extended-first', dest='try_extended_first',
                   action='store_false',
                   help='Disable extended flank optimization; always run '
                        'genome search (reverts to original behaviour)')

    g = ap.add_argument_group('clustering')
    g.add_argument('--cluster-id',  type=float, default=0.80,
                   help='vsearch cluster identity (default: 0.80)')
    g.add_argument('--branch-min',  type=int,   default=3,
                   help='Min seqs per cluster to keep (default: 3)')
    g.add_argument('--max-variants', type=int,  default=3,
                   help='Max LINE variants (branches) to pursue '
                        '(default: 3). Extra clusters are discarded '
                        'by ascending member count.')

    g = ap.add_argument_group('system')
    g.add_argument('--threads',     type=int,   default=4,
                   help='Threads for mafft / sear (default: 4)')

    args = ap.parse_args()

    # ── setup ─────────────────────────────────────────────────────────
    genome = os.path.abspath(args.genome)
    args.outdir = os.path.abspath(args.outdir)
    if not os.path.exists(genome):
        sys.exit(f"Genome not found: {genome}")
    if not os.path.exists(genome + '.fai'):
        print("Indexing genome …")
        run(f"samtools faidx {genome}")

    csizes = load_chrom_sizes(genome + '.fai')
    seed_name, seed_seq = read_fasta(args.seed)
    if len(seed_seq) < 30:
        sys.exit("Seed sequence too short (< 30 bp)")

    os.makedirs(args.outdir, exist_ok=True)
    sear_wd = setup_sear_workdir(args.outdir, genome)

    print(f"LINE_walker")
    print(f"  Seed      : {seed_name} ({len(seed_seq)} bp)")
    print(f"  Genome    : {genome}")
    print(f"  Direction : {args.direction}")
    print(f"  Pipeline  : sear -s {args.search_hits} → top {args.top_hits} "
          f"→ {args.flank} bp flank → cluster @{args.cluster_id} "
          f"→ branch ≥{args.branch_min}  max-variants: {args.max_variants}")
    if args.try_extended_first:
        print(f"  Extended  : try {args.extended_flank} bp from prev hits "
              f"before genome search (--try-extended-first)")
    print()

    log = dict(
        seed=seed_name, seed_len=len(seed_seq), genome=genome,
        direction=args.direction, config=vars(args), branches={}
    )

    # ── walking BFS ───────────────────────────────────────────────────
    # queue items: (accumulated_seq, extension_only, branch_path, next_step,
    #               prev_hits_bed)
    queue    = [(seed_seq, '', '', 1, None)]
    finished = []                       # (acc, ext, branch, reason)

    while queue:
        acc, ext, bp, step, prev_hits_bed = queue.pop(0)

        if step > args.steps:
            print(f"  [branch {bp or 'main'}] max steps reached")
            finished.append((acc, ext, bp, 'max_steps'))
            continue

        # query = tail of accumulated sequence
        qseq = acc[-args.seed_window:] if len(acc) > args.seed_window else acc
        qfa  = os.path.join(args.outdir, f'_q_{bp or "main"}_s{step}.fa')
        write_fasta(qfa, f'q_s{step}', qseq)

        exts, new_hits_bed = walk_step(step, qfa, genome, args.outdir,
                           sear_wd, args, csizes, bp,
                           prev_hits_bed)
        os.remove(qfa)

        if not exts:
            finished.append((acc, ext, bp, 'stopped'))
            continue

        # ── cap total branches at --max-variants ──────────────────
        active_branches = len(queue) + len(finished)
        if len(exts) > 1:
            slots = max(1, args.max_variants - active_branches)
            if slots < len(exts):
                # keep the clusters with the most members (most robust)
                exts = exts[:slots]
                print(f"  [branch cap] keeping {slots} of original "
                      f"clusters (--max-variants {args.max_variants})")

        if len(exts) == 1:
            new_acc = acc + exts[0][0]
            new_ext = ext + exts[0][0]
            queue.append((new_acc, new_ext, bp, step + 1, new_hits_bed))
            print(f"  → {len(new_acc)} bp  (branch {bp or 'main'})")
        else:
            for seq, lbl in exts:
                new_bp  = (bp + lbl) if bp else lbl
                new_acc = acc + seq
                new_ext = ext + seq
                queue.append((new_acc, new_ext, new_bp, step + 1,
                              new_hits_bed))
                print(f"  → fork {new_bp}: {len(new_acc)} bp")
        print()

    # ── output ────────────────────────────────────────────────────────
    cand_path = os.path.join(args.outdir, 'LINE_candidates.fa')
    ext_path  = os.path.join(args.outdir, 'LINE_extensions.fa')
    with open(cand_path, 'w') as fc, open(ext_path, 'w') as fe:
        for acc, ext, bp, reason in finished:
            lbl = bp or 'main'
            fc.write(f'>LINE_{lbl} len={len(acc)} stop={reason}\n')
            for i in range(0, len(acc), 80):
                fc.write(acc[i:i + 80] + '\n')
            if ext:
                fe.write(f'>LINE_{lbl}_ext len={len(ext)} stop={reason}\n')
                for i in range(0, len(ext), 80):
                    fe.write(ext[i:i + 80] + '\n')

    log['branches'] = {
        (b or 'main'): dict(length=len(a), ext_length=len(e), reason=r)
        for a, e, b, r in finished
    }
    log_path = os.path.join(args.outdir, 'walk_log.json')
    with open(log_path, 'w') as fh:
        json.dump(log, fh, indent=2)

    print('=' * 60)
    print('LINE_walker done')
    print(f'  Full sequences : {cand_path}')
    print(f'  Extensions only: {ext_path}')
    for acc, ext, bp, reason in finished:
        print(f'    {bp or "main"}: {len(acc)} bp total, '
              f'{len(ext)} bp novel  ({reason})')
    print(f'  Log: {log_path}')
    print('=' * 60)


if __name__ == '__main__':
    main()
