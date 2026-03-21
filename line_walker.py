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
    → ≥5 members per cluster → majority-rule consensus → extend
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


def trim_anchor_overlap(consensus, direction, anchor_overlap, min_extension=30):
    """Trim the reused anchor portion from a consensus flank sequence.

    Flanks can include a short already-known anchor segment to stabilize
    clustering when sear hit endpoints jitter.  That anchor is trimmed away
    before the extension is appended to the growing LINE consensus.
    """
    if anchor_overlap <= 0 or not consensus:
        return consensus

    trim_bp = min(anchor_overlap, max(0, len(consensus) - min_extension))
    if trim_bp == 0:
        return consensus

    if direction in ('5prime', 'upstream', 'left'):
        return consensus[:-trim_bp]
    return consensus[trim_bp:]


# ═══════════════════════════════════════════════════════════════════════
# sear integration
# ═══════════════════════════════════════════════════════════════════════

def setup_sear_workdir(outdir, genome_abs, work_subdir='_sear_work'):
    """Create shared sear working directory with genome symlink.

    Genome splits (*.2k.part_*.bnk) are created once here and reused
    across all walking steps thanks to sear -k.
    """
    wd = os.path.join(outdir, work_subdir)
    os.makedirs(wd, exist_ok=True)
    link     = os.path.join(wd, 'genome.fa')
    fai_link = os.path.join(wd, 'genome.fa.fai')
    if not os.path.exists(link):
        os.symlink(genome_abs, link)
    if not os.path.exists(fai_link) and os.path.exists(genome_abs + '.fai'):
        os.symlink(genome_abs + '.fai', fai_link)
    return wd


def run_sear(sear_wd, query_fa, search_n, threads, extra_flags=''):
    """Run sear in the shared workdir.  Returns path to output BED."""
    qname = os.path.basename(query_fa)
    q_dest = os.path.join(sear_wd, qname)
    shutil.copy(query_fa, q_dest)

    flags = f' {extra_flags}' if extra_flags else ''
    cmd = f"THREADS={threads} sear{flags} -s {search_n} -k {qname} genome.fa"
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


def deduplicate_hits_by_locus(hits, window_bp):
    """Collapse artificial duplicate hits from same locus.

    Keeps the highest-bitscore hit first and drops later hits that are on the
    same chrom+strand and overlap (or lie within window_bp) of an already
    retained hit.
    """
    if window_bp < 0:
        window_bp = 0
    if not hits:
        return []

    kept = []
    for h in hits:
        is_dup = False
        for k in kept:
            if h['chrom'] != k['chrom'] or h['strand'] != k['strand']:
                continue
            if _interval_distance(h['start'], h['end'],
                                  k['start'], k['end']) <= window_bp:
                is_dup = True
                break
        if not is_dup:
            kept.append(h)
    return kept


def select_hits_for_clustering(all_hits, cluster_hits, min_bitscore_frac):
    """Select a broader, score-filtered pool of hits for flank clustering.

    The old logic used only the first --top-hits entries, which can miss a
    coherent family if the high-scoring tail is spread across many near-ties.
    """
    if not all_hits:
        return []

    top_bitscore = all_hits[0]['bitscore']
    min_bitscore = top_bitscore * min_bitscore_frac
    selected = [h for h in all_hits if h['bitscore'] >= min_bitscore]
    return selected[:cluster_hits]


def build_seed_hit_bank(seed_fa, source_genome, source_sear_wd, outdir, csizes,
                        args):
    """Build pseudo-genome bank from initial seed hits ±bank_flank bp.

    Returns (bank_genome_fa, bank_stats) or (None, bank_stats) on failure.
    """
    bank_dir = os.path.join(outdir, '_seed_hit_bank')
    os.makedirs(bank_dir, exist_ok=True)

    seed_bed_src = run_sear(source_sear_wd, seed_fa, args.search_hits,
                            args.threads)
    if seed_bed_src is None or os.path.getsize(seed_bed_src) == 0:
        return None, dict(status='no_seed_hits')

    seed_bed = os.path.join(bank_dir, 'seed_hits.bed')
    shutil.move(seed_bed_src, seed_bed)
    all_hits = parse_bed7(seed_bed)
    all_hits = deduplicate_hits_by_locus(all_hits, args.dedup_locus_window)
    if not all_hits:
        return None, dict(status='no_seed_hits')

    selected_hits = select_hits_for_clustering(
        all_hits, args.bank_max_seed_hits, args.min_bitscore_frac
    )
    if len(selected_hits) < args.branch_min:
        return None, dict(
            status='too_few_seed_hits',
            seed_hits_total=len(all_hits),
            seed_hits_selected=len(selected_hits)
        )

    regions_bed = os.path.join(bank_dir, 'bank_regions.bed')
    unique_regions = {}
    with open(regions_bed, 'w') as fh:
        for i, h in enumerate(selected_hits):
            clen = csizes.get(h['chrom'])
            if clen is None:
                continue
            rs = max(0, h['start'] - args.bank_flank)
            re = min(clen, h['end'] + args.bank_flank)
            if re - rs < args.seed_window:
                continue
            key = (h['chrom'], rs, re)
            if key in unique_regions:
                continue
            unique_regions[key] = True
            name = f"{h['chrom']}:{rs}-{re}:h{i}"
            fh.write(f"{h['chrom']}\t{rs}\t{re}\t{name}\n")

    if not unique_regions:
        return None, dict(
            status='no_usable_regions',
            seed_hits_total=len(all_hits),
            seed_hits_selected=len(selected_hits)
        )

    bank_genome = os.path.join(bank_dir, 'bank.fa')
    run(f"bedtools getfasta -nameOnly -fi {source_genome} -bed {regions_bed}"
        f" > {bank_genome}")
    run(f"samtools faidx {bank_genome}")

    return bank_genome, dict(
        status='ok',
        seed_hits_total=len(all_hits),
        seed_hits_selected=len(selected_hits),
        bank_regions=len(unique_regions),
        bank_flank=args.bank_flank,
        bank_max_seed_hits=args.bank_max_seed_hits
    )


# ═══════════════════════════════════════════════════════════════════════
# Flank coordinate logic
# ═══════════════════════════════════════════════════════════════════════

def write_flank_bed(hits, direction, flank_size, csizes, out_path,
                    offset=0, anchor_overlap=0):
    """Write flank-only BED.

    direction='3prime'      →  3′ end of the seed-aligned hit (on strand)
    direction='5prime'      →  5′ end of the seed-aligned hit (on strand)
    direction='right'       →  higher genomic coordinates (absolute)
    direction='left'        →  lower genomic coordinates (absolute)

    offset shifts the extraction anchor in the walk direction so that
    successive extended-flank steps extract *new* sequence instead of
    re-extracting the same window.  For the normal sear path the offset
    is always 0 (the hits are already at the current position).

    anchor_overlap includes a short already-known segment across the anchor
    into the extracted sequence.  This stabilizes clustering when sear places
    the same true locus with slightly different local hit endpoints.

    Returns number of valid flanks written.
    """
    n = 0
    with open(out_path, 'w') as fh:
        for i, h in enumerate(hits):
            clen = csizes.get(h['chrom'], 10**9)
            hit_len = max(0, h['end'] - h['start'])
            eff_overlap = anchor_overlap if offset > 0 else min(anchor_overlap,
                                                                 hit_len)
            if direction == 'right':
                anchor = h['end'] + offset
                fs = max(0, anchor - eff_overlap)
                fe = min(anchor + flank_size, clen)
            elif direction == 'left':
                anchor = h['start'] - offset
                fs = max(0, anchor - flank_size)
                fe = min(anchor + eff_overlap, clen)
            elif direction in ('downstream', '3prime'):
                if h['strand'] == '+':
                    anchor = h['end'] + offset
                    fs = max(0, anchor - eff_overlap)
                    fe = min(anchor + flank_size, clen)
                else:
                    anchor = h['start'] - offset
                    fs = max(0, anchor - flank_size)
                    fe = min(anchor + eff_overlap, clen)
            else:  # 5prime / upstream
                if h['strand'] == '+':
                    anchor = h['start'] - offset
                    fs = max(0, anchor - flank_size)
                    fe = min(anchor + eff_overlap, clen)
                else:
                    anchor = h['end'] + offset
                    fs = max(0, anchor - eff_overlap)
                    fe = min(anchor + flank_size, clen)
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
                         csizes, ext_offset=0):
    """Attempt extended flank extraction from previous step's sear hits.

    Uses the same genomic hit locations discovered in the prior step but
    extracts a larger window (args.extended_flank bp) in the walk direction,
    offset by ext_offset bp so that each step advances to new sequence.

    Returns (extensions, hits_bed7) or ([], None) if failed.
    """
    try:
        all_hits = parse_bed7(prev_hits_bed)
    except Exception:
        return [], None

    selected_hits = select_hits_for_clustering(
        all_hits, args.cluster_hits, args.min_bitscore_frac
    )
    if len(selected_hits) < args.branch_min:
        return [], None

    ext_sd = os.path.join(sd, 'ext')
    os.makedirs(ext_sd, exist_ok=True)

    ext_flank_bed = os.path.join(ext_sd, 'flanks.bed')
    nf = write_flank_bed(selected_hits, args.direction, args.extended_flank,
                         csizes, ext_flank_bed, offset=ext_offset,
                         anchor_overlap=args.anchor_overlap)
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
    cluster_members_all = []
    for cf in cluster_files:
        members = read_multi_fasta(os.path.join(clu_dir, cf))
        cluster_members_all.append(members)
        if len(members) >= args.branch_min:
            qualifying.append(members)

    if not qualifying:
        largest = max(cluster_members_all, key=len) if cluster_members_all else []
        if len(largest) >= args.rescue_min_members:
            qualifying = [largest]
            print(f"  [{tag}] extended: no clusters ≥{args.branch_min}; "
                  f"rescue largest cluster ({len(largest)} seqs)")
        elif os.path.exists(ext_aligned_fa):
            aligned = read_multi_fasta(ext_aligned_fa)
            cons = majority_consensus([s for _, s in aligned])
            cons = trim_anchor_overlap(cons, args.direction,
                                       args.anchor_overlap)
            if len(cons) >= 30:
                ext_hits_bed = os.path.join(ext_sd, 'extended_hits.bed')
                write_bed7(selected_hits, ext_hits_bed)
                ccons = os.path.join(ext_sd, 'consensus_aln.fa')
                write_fasta(ccons, f'step{step_num}_aln_ext', cons)
                print(f"  [{tag}] extended alignment fallback: "
                      f"{len(aligned)} seqs → {len(cons)} bp consensus")
                return [(cons, 'A')], ext_hits_bed
            print(f"  [{tag}] extended: alignment fallback too short "
                  f"({len(cons)} bp) — skip")
            return [], None
        else:
            print(f"  [{tag}] extended: no clusters ≥{args.branch_min} — stop")
            return [], None

    ext_hits_bed = os.path.join(ext_sd, 'extended_hits.bed')
    write_bed7(selected_hits, ext_hits_bed)

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
        cons = trim_anchor_overlap(cons, args.direction, args.anchor_overlap)
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
              branch='', prev_hits_bed=None, ext_offset=0):
    """Execute one walking step.

    If args.try_extended_first is True and prev_hits_bed is available, first
    attempts to extract extended flanks (args.extended_flank bp) from the
    previous step's sear hits before falling back to a full genome search.

    ext_offset is the number of bp already walked past the reference hits
    in prev_hits_bed; used to shift the extended-flank window so each step
    extracts new sequence.

    Returns (extensions, next_hits_bed, used_extended) where extensions is
    a list of (consensus_seq, cluster_label), next_hits_bed is the BED to
    pass to the following step, and used_extended indicates whether the
    extended-flank fast path was taken (True) or sear ran (False).
    """
    tag = f"step_{step_num:03d}" + (f"_branch_{branch}" if branch else "")
    sd  = os.path.join(outdir, tag)
    os.makedirs(sd, exist_ok=True)
    shutil.copy(query_fa, os.path.join(sd, 'query.fa'))

    stats = dict(step=step_num, branch=branch, tag=tag)
    fell_back_from_extended = False

    if (args.try_extended_first
            and prev_hits_bed is not None
            and os.path.exists(prev_hits_bed)):
        print(f"  [{tag}] trying extended flanks "
              f"({args.extended_flank} bp, direction={args.direction}, "
              f"offset={ext_offset}) ...")
        ext_exts, ext_hits_bed = _try_extended_flanks(
            step_num, tag, prev_hits_bed, genome, sd, args, csizes,
            ext_offset=ext_offset
        )
        if ext_exts:
            stats['status'] = 'extended_from_prev_hits'
            stats['extended_flank_size'] = args.extended_flank
            stats['ext_offset'] = ext_offset
            stats['extensions'] = [
                dict(branch=l, length=len(s)) for s, l in ext_exts
            ]
            _dump_stats(stats, sd)
            return ext_exts, ext_hits_bed, True
        fell_back_from_extended = True
        print(f"  [{tag}] extended extraction failed — falling back to sear")

    # 1 ── sear search ─────────────────────────────────────────────────
    print(f"  [{tag}] sear -s {args.search_hits} ...")
    bed_src = run_sear(sear_wd, query_fa, args.search_hits, args.threads)
    if bed_src is None or os.path.getsize(bed_src) == 0:
        print(f"  [{tag}] no hits")
        stats['status'] = 'no_hits'
        _dump_stats(stats, sd)
        return [], None, False

    bed_dst = os.path.join(sd, 'sear_hits.bed')
    shutil.move(bed_src, bed_dst)

    # 2 ── select top hits by bitscore ─────────────────────────────────
    all_hits = parse_bed7(bed_dst)
    stats['hits_before_dedup'] = len(all_hits)
    all_hits = deduplicate_hits_by_locus(all_hits, args.dedup_locus_window)
    stats['hits_after_dedup'] = len(all_hits)
    if not fell_back_from_extended and prev_hits_bed is not None and os.path.exists(prev_hits_bed):
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
            return [], bed_dst, False
    elif fell_back_from_extended:
        print(f"  [{tag}] skipping continuity filter (sear fallback after extended mode)")

    top = all_hits[:args.top_hits]
    cluster_hits = select_hits_for_clustering(
        all_hits, args.cluster_hits, args.min_bitscore_frac
    )
    stats['hits_total'] = len(all_hits)
    stats['hits_used']  = len(top)
    stats['cluster_hits_used'] = len(cluster_hits)
    stats['cluster_min_bitscore_frac'] = args.min_bitscore_frac
    print(f"  [{tag}] {len(all_hits)} hits → top {len(top)} by bitscore")
    print(f"  [{tag}] clustering pool: {len(cluster_hits)} hits "
          f"(bitscore ≥ top×{args.min_bitscore_frac})")

    if len(cluster_hits) < args.branch_min:
        stats['status'] = 'too_few_hits'
        _dump_stats(stats, sd)
        return [], bed_dst, False

    # 3 ── extract directional flanks ──────────────────────────────────
    flank_bed = os.path.join(sd, 'flanks.bed')
    nf = write_flank_bed(cluster_hits, args.direction, args.flank, csizes,
                         flank_bed, anchor_overlap=args.anchor_overlap)
    if nf < args.branch_min:
        print(f"  [{tag}] only {nf} usable flanks (need {args.branch_min})")
        stats['status'] = 'too_few_flanks'
        _dump_stats(stats, sd)
        return [], bed_dst, False

    flanks_fa = os.path.join(sd, 'flanks.fa')
    run(f"bedtools getfasta -s -nameOnly -fi {genome} -bed {flank_bed}"
        f" > {flanks_fa}", cwd=sd)

    seqs = read_multi_fasta(flanks_fa)
    if len(seqs) < args.branch_min:
        stats['status'] = 'too_few_flanks'
        _dump_stats(stats, sd)
        return [], bed_dst, False
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
    cluster_members_all = []
    for cf in cluster_files:
        members = read_multi_fasta(os.path.join(clu_dir, cf))
        cluster_members_all.append(members)
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
        largest = max(cluster_members_all, key=len) if cluster_members_all else []
        if len(largest) >= args.rescue_min_members:
            qualifying = [largest]
            stats['rescue_single_cluster'] = True
            stats['rescue_cluster_size'] = len(largest)
            print(f"  [{tag}] no clusters ≥{args.branch_min}; "
                  f"rescue largest cluster ({len(largest)} seqs)")
        elif os.path.exists(aligned_fa):
            aln_entries = read_multi_fasta(aligned_fa)
            aln_cons = majority_consensus([s for _, s in aln_entries])
            aln_cons = trim_anchor_overlap(aln_cons, args.direction,
                                           args.anchor_overlap)
            if len(aln_cons) >= 30:
                ccons = os.path.join(sd, 'consensus_aln.fa')
                write_fasta(ccons, f'step{step_num}_aln', aln_cons)
                extensions = [(aln_cons, 'A')]
                stats['status'] = 'alignment_fallback'
                stats['alignment_consensus_len'] = len(aln_cons)
                stats['extensions'] = [
                    dict(branch='A', length=len(aln_cons))
                ]
                print(f"  [{tag}] alignment fallback: "
                      f"{len(aln_entries)} seqs → {len(aln_cons)} bp consensus")
                _dump_stats(stats, sd)
                return extensions, bed_dst, False
            print(f"  [{tag}] alignment fallback too short "
                  f"({len(aln_cons)} bp)")
            stats['status'] = 'no_qualifying_cluster'
            _dump_stats(stats, sd)
            return [], bed_dst, False
        else:
            print(f"  [{tag}] no clusters ≥{args.branch_min} — stop")
            stats['status'] = 'no_qualifying_cluster'
            _dump_stats(stats, sd)
            return [], bed_dst, False

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
        cons = trim_anchor_overlap(cons, args.direction, args.anchor_overlap)

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
    return extensions, bed_dst, False


def collect_spanning_loci(run_outdir, branch, merge_window):
    """Collect step-1 sear hit loci and count how many later steps also
    placed a hit within merge_window of each anchor.

    Only step-1 sear hits are used as anchors (they define the genuine
    seed-matching loci).  Every subsequent step's hits are checked against
    those anchors; any hit within merge_window is credited to the anchor.

    Returns groups sorted by step_count desc then bitscore desc.
    """
    bp_suffix = f'_branch_{branch}' if branch else ''

    # Step 1 anchors
    step1_bed = os.path.join(run_outdir, f'step_001{bp_suffix}', 'sear_hits.bed')
    if not os.path.exists(step1_bed) or os.path.getsize(step1_bed) == 0:
        return []
    anchors = parse_bed7(step1_bed)   # already sorted by bitscore desc
    if not anchors:
        return []

    groups = [
        {'chrom': h['chrom'], 'strand': h['strand'],
         'min_start': h['start'], 'max_end': h['end'],
         'bitscore': h['bitscore'], 'steps': {1}, 'n_hits': 1}
        for h in anchors
    ]

    step = 2
    while True:
        tag = f'step_{step:03d}{bp_suffix}'
        sd = os.path.join(run_outdir, tag)
        if not os.path.exists(sd):
            break
        for bed_cand in [
            os.path.join(sd, 'ext', 'extended_hits.bed'),
            os.path.join(sd, 'sear_hits.bed'),
        ]:
            if os.path.exists(bed_cand) and os.path.getsize(bed_cand) > 0:
                for h in parse_bed7(bed_cand):
                    for g in groups:
                        if (g['chrom'] != h['chrom']
                                or g['strand'] != h['strand']):
                            continue
                        if _interval_distance(h['start'], h['end'],
                                              g['min_start'],
                                              g['max_end']) <= merge_window:
                            g['max_end'] = max(g['max_end'], h['end'])
                            g['min_start'] = min(g['min_start'], h['start'])
                            g['steps'].add(step)
                            g['n_hits'] += 1
                            break
                break
        step += 1

    for g in groups:
        g['step_count'] = len(g['steps'])

    groups.sort(key=lambda g: (-g['step_count'], -g['bitscore']))
    return groups


def export_spanning_evidence(finished, run_outdir, walk_genome, walk_csizes,
                             args, walking_5prime=True):
    """For each branch, find loci with hits spanning the most walk steps,
    extract the LINE body region from the genome, and align with consensus.

    Extraction is direction-aware: cons_len is added only into the walk
    direction from the seed anchor; a small buffer is added on the other side.

    Writes per branch:
      spanning_evidence/<branch>.bed6
      spanning_evidence/<branch>.fa
      spanning_evidence/<branch>_with_consensus_aligned.fa

    Returns dict: label → (n_loci, aligned_fa_path) or None.
    """
    evid_dir = os.path.join(run_outdir, 'spanning_evidence')
    os.makedirs(evid_dir, exist_ok=True)

    # merge_window = max_jump: same guarantee used by continuity filter
    merge_window = args.max_jump

    results = {}
    for acc, _ext, bp, _reason in finished:
        lbl = bp or 'main'
        cons_len = len(acc)
        groups = collect_spanning_loci(run_outdir, bp, merge_window)
        if not groups:
            print(f"  Spanning [{lbl}]: no step-1 loci found")
            results[lbl] = None
            continue

        top = groups[:args.best_loci]
        print(f"  Spanning [{lbl}]: top {len(top)} loci "
              f"(best: {top[0]['step_count']} steps covered, "
              f"anchor bitscore={top[0]['bitscore']:.1f})")

        bed6 = os.path.join(evid_dir, f'{lbl}.bed6')
        FLANK = 300  # buffer on the seed (non-walk) side
        with open(bed6, 'w') as fh:
            for i, g in enumerate(top, start=1):
                clen = walk_csizes.get(g['chrom'], 10**9)
                # Extend cons_len into the walk direction from the seed
                # anchor, FLANK bp buffer on the other side.
                # 5prime walk on + strand: LINE is to the LEFT of seed
                # 5prime walk on - strand: LINE is to the RIGHT of seed
                if walking_5prime:
                    if g['strand'] == '+':
                        s = max(0, g['min_start'] - cons_len - FLANK)
                        e = min(clen, g['max_end'] + FLANK)
                    else:
                        s = max(0, g['min_start'] - FLANK)
                        e = min(clen, g['max_end'] + cons_len + FLANK)
                else:  # 3prime
                    if g['strand'] == '+':
                        s = max(0, g['min_start'] - FLANK)
                        e = min(clen, g['max_end'] + cons_len + FLANK)
                    else:
                        s = max(0, g['min_start'] - cons_len - FLANK)
                        e = min(clen, g['max_end'] + FLANK)
                name = f"locus{i}_{g['chrom']}({s}-{e})"
                fh.write(f"{g['chrom']}\t{s}\t{e}\t{name}\t0\t{g['strand']}\n")

        loci_fa = os.path.join(evid_dir, f'{lbl}.fa')
        run(f"bedtools getfasta -s -nameOnly -fi {walk_genome}"
            f" -bed {bed6} > {loci_fa}")

        loci_entries = read_multi_fasta(loci_fa)
        if not loci_entries:
            print(f"  Spanning [{lbl}]: bedtools returned no sequences")
            results[lbl] = None
            continue

        to_align = os.path.join(evid_dir, f'{lbl}_with_consensus.fa')
        aligned  = os.path.join(evid_dir, f'{lbl}_with_consensus_aligned.fa')
        # Seed region uppercase, extensions lowercase — aids visual inspection.
        if walking_5prime:
            mixed_acc = acc[:len(_ext)].lower() + acc[len(_ext):]
        else:
            mixed_acc = acc[:len(acc) - len(_ext)] + acc[len(acc) - len(_ext):].lower()
        write_multi_fasta(to_align, [(f'LINE_{lbl}_consensus', mixed_acc)]
                          + loci_entries)

        print(f"  Spanning [{lbl}]: MAFFT aligning "
              f"{len(loci_entries)} loci + consensus ...")
        run(f"mafft --thread {args.threads} --localpair --maxiterate 1000 "
            f"--ep 0.123 --nuc --reorder --quiet {to_align} > {aligned}")

        print(f"  Spanning [{lbl}]: alignment → {aligned}")
        results[lbl] = (len(loci_entries), aligned)

    return results


def export_best_loci_for_candidates(finished, spanning_results, source_genome,
                                    run_outdir, args):
    """For each final consensus, export the top spanning loci (already found
    during the walk) and align them with the consensus.

    Re-uses loci from *spanning_results* so no expensive genome re-search
    is needed.

    Writes, per branch:
      - best_loci_<branch>.tsv   (from spanning evidence)
      - best_loci_<branch>.fa
      - best_loci_<branch>_with_consensus_aligned.fa
    """
    loci_dir = os.path.join(run_outdir, 'best_loci')
    os.makedirs(loci_dir, exist_ok=True)

    evid_dir = os.path.join(run_outdir, 'spanning_evidence')
    summary = {}
    for acc, _ext, bp, _reason in finished:
        lbl = bp or 'main'
        tag = f'LINE_{lbl}'

        span = spanning_results.get(lbl)
        if not span:
            print(f"  Best loci [{lbl}]: no spanning loci available")
            summary[lbl] = dict(status='no_spanning_loci')
            continue

        n_loci, span_aligned = span

        # Copy spanning evidence files into best_loci dir
        span_fa = os.path.join(evid_dir, f'{lbl}.fa')
        span_bed6 = os.path.join(evid_dir, f'{lbl}.bed6')
        if not os.path.exists(span_fa):
            summary[lbl] = dict(status='no_spanning_fa')
            continue

        loci_fa = os.path.join(loci_dir, f'best_loci_{lbl}.fa')
        shutil.copy(span_fa, loci_fa)
        if os.path.exists(span_bed6):
            shutil.copy(span_bed6,
                        os.path.join(loci_dir, f'best_loci_{lbl}.bed6'))

        # Align consensus with the top loci
        aligned = os.path.join(
            loci_dir, f'best_loci_{lbl}_with_consensus_aligned.fa'
        )
        to_align = os.path.join(
            loci_dir, f'best_loci_{lbl}_with_consensus.fa'
        )
        loci_entries = read_multi_fasta(loci_fa)
        write_multi_fasta(
            to_align,
            [(f'{tag}_consensus', acc)] + loci_entries
        )
        run(f"mafft --thread {args.threads} --localpair --maxiterate 1000 "
            f"--ep 0.123 --nuc --reorder --quiet {to_align} > {aligned}")

        print(f"  Best loci [{lbl}]: {n_loci} loci from spanning evidence")
        print(f"  Best loci [{lbl}]: alignment → {aligned}")

        summary[lbl] = dict(
            status='ok',
            top_hits=n_loci,
        )

    return summary


def run_walk_direction(direction, seed_name, seed_seq, walk_genome, walk_csizes,
                       sear_wd, source_genome, source_sear_wd, args, run_outdir,
                       bank_stats=None):
    """Run one complete LINE walking pass for one direction."""
    os.makedirs(run_outdir, exist_ok=True)

    old_direction = args.direction
    if direction == 'upstream':
        effective_direction = '5prime'
    elif direction == 'downstream':
        effective_direction = '3prime'
    else:
        effective_direction = direction
    args.direction = effective_direction
    walking_5prime = effective_direction == '5prime'
    try:
        print(f"LINE_walker")
        print(f"  Seed      : {seed_name} ({len(seed_seq)} bp)")
        print(f"  Genome    : {walk_genome}")
        print(f"  Direction : {effective_direction} (seed-relative)")
        print(f"  Pipeline  : sear -s {args.search_hits} → top {args.top_hits} "
              f"→ {args.flank} bp flank (+{args.anchor_overlap} bp anchor) "
              f"→ cluster @{args.cluster_id} "
              f"→ branch ≥{args.branch_min}  max-variants: {args.max_variants}")
        if args.try_extended_first:
            print(f"  Extended  : try {args.extended_flank} bp from prev hits "
                  f"before genome search (--try-extended-first)")
        if args.use_seed_hit_bank:
            print(f"  Seed bank : enabled (±{args.bank_flank} bp from up to "
                  f"{args.bank_max_seed_hits} near-top seed hits)")
            if bank_stats is not None:
                print(f"              selected {bank_stats.get('seed_hits_selected', 0)} "
                      f"hits → {bank_stats.get('bank_regions', 0)} regions")
        print()

        log = dict(
            seed=seed_name, seed_len=len(seed_seq), genome=walk_genome,
            direction=effective_direction, config=vars(args), branches={}
        )
        if bank_stats is not None:
            log['seed_hit_bank'] = bank_stats

        queue = [(seed_seq, '', '', 1, None, 0)]
        finished = []

        while queue:
            acc, ext, bp, step, prev_hits_bed, ext_offset = queue.pop(0)

            if step > args.steps:
                print(f"  [branch {bp or 'main'}] max steps reached")
                finished.append((acc, ext, bp, 'max_steps'))
                continue

            if walking_5prime:
                qseq = acc[:args.seed_window]
            else:
                qseq = acc[-args.seed_window:]
            qfa = os.path.join(run_outdir, f'_q_{bp or "main"}_s{step}.fa')
            write_fasta(qfa, f'q_s{step}', qseq)

            exts, new_hits_bed, used_extended = walk_step(
                step, qfa, walk_genome, run_outdir, sear_wd, args,
                walk_csizes, bp, prev_hits_bed, ext_offset=ext_offset)
            os.remove(qfa)

            if not exts:
                finished.append((acc, ext, bp, 'stopped'))
                continue

            active_branches = len(queue) + len(finished)
            if len(exts) > 1:
                slots = max(1, args.max_variants - active_branches)
                if slots < len(exts):
                    exts = exts[:slots]
                    print(f"  [branch cap] keeping {slots} of original "
                          f"clusters (--max-variants {args.max_variants})")

            if len(exts) == 1:
                if walking_5prime:
                    new_acc = exts[0][0] + acc
                    new_ext = exts[0][0] + ext
                else:
                    new_acc = acc + exts[0][0]
                    new_ext = ext + exts[0][0]
                if used_extended:
                    new_offset = ext_offset + len(exts[0][0])
                else:
                    new_offset = len(exts[0][0])
                queue.append((new_acc, new_ext, bp, step + 1, new_hits_bed,
                              new_offset))
                print(f"  → {len(new_acc)} bp  (branch {bp or 'main'})")
            else:
                for seq, lbl in exts:
                    new_bp = (bp + lbl) if bp else lbl
                    if walking_5prime:
                        new_acc = seq + acc
                        new_ext = seq + ext
                    else:
                        new_acc = acc + seq
                        new_ext = ext + seq
                    if used_extended:
                        new_offset = ext_offset + len(seq)
                    else:
                        new_offset = len(seq)
                    queue.append((new_acc, new_ext, new_bp, step + 1,
                                  new_hits_bed, new_offset))
                    print(f"  → fork {new_bp}: {len(new_acc)} bp")
            print()

        # ── export spanning genomic evidence ────────────────────────────
        spanning_results = export_spanning_evidence(
            finished, run_outdir, walk_genome, walk_csizes, args, walking_5prime
        )

        cand_path = os.path.join(run_outdir, 'LINE_candidates.fa')
        ext_path = os.path.join(run_outdir, 'LINE_extensions.fa')
        with open(cand_path, 'w') as fc, open(ext_path, 'w') as fe:
            for acc, ext, bp, reason in finished:
                lbl = bp or 'main'
                # Seed region uppercase, novel extensions lowercase.
                if walking_5prime:
                    mixed = acc[:len(ext)].lower() + acc[len(ext):]
                else:
                    mixed = acc[:len(acc) - len(ext)] + acc[len(acc) - len(ext):].lower()
                fc.write(f'>LINE_{lbl} len={len(mixed)} stop={reason}\n')
                for i in range(0, len(mixed), 80):
                    fc.write(mixed[i:i + 80] + '\n')
                if ext:
                    fe.write(f'>LINE_{lbl}_ext len={len(ext)} stop={reason}\n')
                    for i in range(0, len(ext), 80):
                        fe.write(ext[i:i + 80] + '\n')

        log['branches'] = {
            (b or 'main'): dict(length=len(a), ext_length=len(e), reason=r)
            for a, e, b, r in finished
        }
        log['best_loci'] = export_best_loci_for_candidates(
            finished, spanning_results, source_genome, run_outdir, args
        )

        log_path = os.path.join(run_outdir, 'walk_log.json')
        with open(log_path, 'w') as fh:
            json.dump(log, fh, indent=2)

        print('=' * 60)
        print('LINE_walker done')
        print(f'  Full sequences : {cand_path}')
        print(f'  Extensions only: {ext_path}')
        print(f'  Best loci      : {os.path.join(run_outdir, "best_loci")}')
        for acc, ext, bp, reason in finished:
            print(f'    {bp or "main"}: {len(acc)} bp total, '
                  f'{len(ext)} bp novel  ({reason})')
        for lbl, res in spanning_results.items():
            if res:
                n, aln = res
                print(f'  Spanning [{lbl}] : {n} loci aligned → {aln}')
        print(f'  Log: {log_path}')
        print('=' * 60)

        return finished, log_path
    finally:
        args.direction = old_direction


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
                    choices=('5prime', '3prime', 'upstream', 'downstream',
                             'left', 'right', 'both'),
                    help='Walk direction. Preferred seed-relative modes: '
                        '5prime/3prime (upstream/downstream are aliases). '
                        'left/right are absolute genomic-coordinate modes. '
                        'Use both to run 5prime and 3prime in separate '
                        'subdirectories. (default: upstream)')

    g = ap.add_argument_group('search tuning')
    g.add_argument('--steps',       type=int,   default=30,
                   help='Max walking steps (default: 30)')
    g.add_argument('--search-hits', type=int,   default=50,
                   help='sear: stop after N hits (default: 50)')
    g.add_argument('--top-hits',    type=int,   default=20,
                   help='Keep N best hits by bitscore for reporting and '
                        'inspection (default: 20)')
    g.add_argument('--best-loci', type=int, default=5,
                   help='For each final consensus, report top N best genomic '
                        'loci by bitscore and create consensus+loci '
                        'alignment (default: 5)')
    g.add_argument('--best-loci-search-hits', type=int, default=200,
                   help='sear hit cap when mapping final consensus back to '
                        'genome for best-loci reporting (default: 200)')
    g.add_argument('--cluster-hits', type=int, default=50,
                   help='Max number of near-top hits used for flank '
                        'clustering/consensus (default: 50)')
    g.add_argument('--min-bitscore-frac', type=float, default=0.90,
                   help='Keep hits with bitscore at least top_bitscore × F '
                        'for clustering pool selection (default: 0.90)')
    g.add_argument('--flank',       type=int,   default=150,
                   help='Flank extraction size in bp (default: 150)')
    g.add_argument('--anchor-overlap', type=int, default=60,
                   help='Include up to N bp of already-matched sequence '
                        'across the flank anchor to stabilize clustering; '
                        'trimmed off before extension (default: 60)')
    g.add_argument('--seed-window', type=int,   default=200,
                   help='Use last N bp of accumulated seq as query '
                        '(default: 200)')
    g.add_argument('--max-jump', type=int, default=5000,
                   help='Max genomic distance (bp) from previous-step hits '
                        'allowed when continuity filtering fresh sear hits '
                        '(default: 5000)')
    g.add_argument('--dedup-locus-window', type=int, default=30,
                   help='Collapse same-locus near-duplicate sear hits when '
                        'distance/overlap is within N bp (default: 30)')
    g.add_argument('--use-seed-hit-bank', action='store_true',
                   help='Build a pseudo-genome from initial seed-hit loci '
                        '(± --bank-flank) and perform all walking only in '
                        'that extracted bank')
    g.add_argument('--bank-flank', type=int, default=20000,
                   help='Half-window size in bp around each selected seed hit '
                        'when building pseudo-genome bank (default: 20000)')
    g.add_argument('--bank-max-seed-hits', type=int, default=100,
                   help='Max near-top initial seed hits used to build '
                        'pseudo-genome bank (default: 100)')

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
    g.add_argument('--branch-min',  type=int,   default=5,
                   help='Min seqs per cluster to keep (default: 5)')
    g.add_argument('--rescue-min-members', type=int, default=5,
                   help='If no cluster reaches --branch-min, rescue only the '
                        'single largest cluster when it has at least this '
                        'many members (default: 5)')
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
    source_sear_wd = setup_sear_workdir(args.outdir, genome,
                                        '_sear_work_genome')
    sear_wd = source_sear_wd

    walk_genome = genome
    walk_csizes = csizes
    bank_stats = None
    if args.use_seed_hit_bank:
        print("Preparing seed-hit pseudo-genome bank ...")
        bank_genome, bank_stats = build_seed_hit_bank(
            args.seed, genome, sear_wd, args.outdir, csizes, args
        )
        if bank_genome is None:
            reason = bank_stats.get('status', 'unknown') if bank_stats else 'unknown'
            sys.exit(f"Failed to build seed-hit bank: {reason}")
        walk_genome = os.path.abspath(bank_genome)
        walk_csizes = load_chrom_sizes(walk_genome + '.fai')
        sear_wd = setup_sear_workdir(args.outdir, walk_genome,
                                     '_sear_work_bank')

    if args.direction == 'both':
        directions = ['5prime', '3prime']
    else:
        directions = [args.direction]

    all_results = {}
    for direction in directions:
        if args.direction == 'both':
            run_outdir = os.path.join(args.outdir, direction)
            print(f"\n=== Running direction: {direction} ===\n")
        else:
            run_outdir = args.outdir

        finished, log_path = run_walk_direction(
            direction=direction,
            seed_name=seed_name,
            seed_seq=seed_seq,
            walk_genome=walk_genome,
            walk_csizes=walk_csizes,
            sear_wd=sear_wd,
            source_genome=genome,
            source_sear_wd=source_sear_wd,
            args=args,
            run_outdir=run_outdir,
            bank_stats=bank_stats
        )
        all_results[direction] = finished

    # ── merge 5' and 3' walks for -d both ─────────────────────────────
    if args.direction == 'both':
        fivep = all_results.get('5prime', [])
        threep = all_results.get('3prime', [])
        if fivep and threep:
            merged_dir = os.path.join(args.outdir, 'merged')
            os.makedirs(merged_dir, exist_ok=True)
            merged_entries = []
            for acc5, ext5, bp5, reason5 in fivep:
                for acc3, ext3, bp3, reason3 in threep:
                    lbl5 = bp5 or 'main'
                    lbl3 = bp3 or 'main'
                    merged_name = f'LINE_5p{lbl5}_3p{lbl3}'
                    # acc5 = [5'_extensions][seed]
                    # acc3 = [seed][3'_extensions]
                    # ext3 = 3' extensions only
                    merged_seq = acc5 + ext3
                    merged_entries.append((merged_name, merged_seq))
                    print(f"  merged {merged_name}: {len(ext5)} bp 5' ext + "
                          f"{len(seed_seq)} bp seed + {len(ext3)} bp 3' ext "
                          f"= {len(merged_seq)} bp")

            merged_fa = os.path.join(merged_dir, 'LINE_merged.fa')
            write_multi_fasta(merged_fa, [
                (f'{name} len={len(seq)}', seq)
                for name, seq in merged_entries
            ])
            print(f"\n  Merged sequences : {merged_fa}")
            print(f"  Combinations     : {len(merged_entries)}")


if __name__ == '__main__':
    main()
