#!/usr/bin/env python3
"""
line_walker.py — Iterative LINE reconstruction from a SINE tail seed.

Given a SINE tail sequence and a genome, walks along genomic LINE copies
to reconstruct the full LINE element step by step.

Pipeline per step:
  seed → sear -s 15 -k  (stop at 15 hits, reuse genome splits)
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


# ═══════════════════════════════════════════════════════════════════════
# Flank coordinate logic
# ═══════════════════════════════════════════════════════════════════════

def write_flank_bed(hits, direction, flank_size, csizes, out_path):
    """Write strand-aware flank-only BED.

    direction='downstream'  →  3′ end of the hit (on its strand)
    direction='upstream'    →  5′ end of the hit (on its strand)

    Returns number of valid flanks written.
    """
    n = 0
    with open(out_path, 'w') as fh:
        for i, h in enumerate(hits):
            clen = csizes.get(h['chrom'], 10**9)
            if direction == 'downstream':
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


def walk_step(step_num, query_fa, genome, outdir, sear_wd, args, csizes,
              branch=''):
    """Execute one walking step.

    Returns list of (consensus_seq, cluster_label) for next steps.
    Empty list → this branch stops.
    """
    tag = f"step_{step_num:03d}" + (f"_branch_{branch}" if branch else "")
    sd  = os.path.join(outdir, tag)
    os.makedirs(sd, exist_ok=True)
    shutil.copy(query_fa, os.path.join(sd, 'query.fa'))

    stats = dict(step=step_num, branch=branch, tag=tag)

    # 1 ── sear search ─────────────────────────────────────────────────
    print(f"  [{tag}] sear -s {args.search_hits} ...")
    bed_src = run_sear(sear_wd, query_fa, args.search_hits, args.threads)
    if bed_src is None or os.path.getsize(bed_src) == 0:
        print(f"  [{tag}] no hits")
        stats['status'] = 'no_hits'
        _dump_stats(stats, sd)
        return []

    bed_dst = os.path.join(sd, 'sear_hits.bed')
    shutil.move(bed_src, bed_dst)

    # 2 ── select top hits by bitscore ─────────────────────────────────
    all_hits = parse_bed7(bed_dst)
    top = all_hits[:args.top_hits]
    stats['hits_total'] = len(all_hits)
    stats['hits_used']  = len(top)
    print(f"  [{tag}] {len(all_hits)} hits → top {len(top)} by bitscore")

    if len(top) < args.branch_min:
        stats['status'] = 'too_few_hits'
        _dump_stats(stats, sd)
        return []

    # 3 ── extract directional flanks ──────────────────────────────────
    flank_bed = os.path.join(sd, 'flanks.bed')
    nf = write_flank_bed(top, args.direction, args.flank, csizes, flank_bed)
    if nf < args.branch_min:
        print(f"  [{tag}] only {nf} usable flanks (need {args.branch_min})")
        stats['status'] = 'too_few_flanks'
        _dump_stats(stats, sd)
        return []

    flanks_fa = os.path.join(sd, 'flanks.fa')
    run(f"bedtools getfasta -s -nameOnly -fi {genome} -bed {flank_bed}"
        f" > {flanks_fa}", cwd=sd)

    seqs = read_multi_fasta(flanks_fa)
    if len(seqs) < args.branch_min:
        stats['status'] = 'too_few_flanks'
        _dump_stats(stats, sd)
        return []
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
    # If no cluster has ≥ branch_min members, fall back to treating
    # ALL flanks as a single group (use the MAFFT alignment already done).
    if not qualifying:
        print(f"  [{tag}] no clusters ≥{args.branch_min} — "
              f"using all {len(seqs)} flanks as single group")
        qualifying = [seqs]
        stats['fallback_single_group'] = True

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
    return extensions


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
                    choices=('upstream', 'downstream'),
                    help='Walk direction relative to seed on its strand '
                         '(default: upstream — toward LINE 5′ end)')

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

    g = ap.add_argument_group('clustering')
    g.add_argument('--cluster-id',  type=float, default=0.80,
                   help='vsearch cluster identity (default: 0.80)')
    g.add_argument('--branch-min',  type=int,   default=3,
                   help='Min seqs per cluster to keep (default: 3)')

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
          f"→ branch ≥{args.branch_min}")
    print()

    log = dict(
        seed=seed_name, seed_len=len(seed_seq), genome=genome,
        direction=args.direction, config=vars(args), branches={}
    )

    # ── walking BFS ───────────────────────────────────────────────────
    # queue items: (accumulated_seq, extension_only, branch_path, next_step)
    queue    = [(seed_seq, '', '', 1)]
    finished = []                       # (acc, ext, branch, reason)

    while queue:
        acc, ext, bp, step = queue.pop(0)

        if step > args.steps:
            print(f"  [branch {bp or 'main'}] max steps reached")
            finished.append((acc, ext, bp, 'max_steps'))
            continue

        # query = tail of accumulated sequence
        qseq = acc[-args.seed_window:] if len(acc) > args.seed_window else acc
        qfa  = os.path.join(args.outdir, f'_q_{bp or "main"}_s{step}.fa')
        write_fasta(qfa, f'q_s{step}', qseq)

        exts = walk_step(step, qfa, genome, args.outdir, sear_wd,
                         args, csizes, bp)
        os.remove(qfa)

        if not exts:
            finished.append((acc, ext, bp, 'stopped'))
            continue

        if len(exts) == 1:
            new_acc = acc + exts[0][0]
            new_ext = ext + exts[0][0]
            queue.append((new_acc, new_ext, bp, step + 1))
            print(f"  → {len(new_acc)} bp  (branch {bp or 'main'})")
        else:
            for seq, lbl in exts:
                new_bp  = (bp + lbl) if bp else lbl
                new_acc = acc + seq
                new_ext = ext + seq
                queue.append((new_acc, new_ext, new_bp, step + 1))
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
