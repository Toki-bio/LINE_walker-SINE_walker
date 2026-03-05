#!/usr/bin/env python3
"""
sine_walker.py — Iterative SINE consensus reconstruction from a seed.

Given a SINE fragment (seed) and a genome, walks bidirectionally to
reconstruct the full SINE element step by step.

Pipeline per arm per step:
  current fragment → sear (reuse genome splits)
  → sort by bitscore, take top --top-hits
  → extract --fragment bp directional flank (strand-aware)
  → MAFFT align flanks
  → majority-rule consensus → next fragment → grow arm

Runs 5′ and 3′ arms independently.

Stop condition per arm:
  - fewer than --branch-min usable flanks / consensus too short
  → fallback: align accumulated top-10 whole hits against current consensus

Final output:
  - SINE_reconstruction.fa  — [5′ extension] + [seed] + [3′ extension]
  - SINE_top_hits.fa        — top genomic hits vs reconstruction
  - SINE_top_hits_aligned.fa — mafft alignment of top hits + reconstruction
  - walk_log.json

Dependencies: sear, samtools, bedtools, mafft
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


def majority_consensus(aligned_seqs, min_depth=2, min_freq=0.0):
    """Majority-rule consensus from aligned sequence strings.

    Positions where gaps outnumber bases are skipped (natural trimming).
    If min_freq > 0, positions where the top base frequency is below
    min_freq are also skipped (trims noisy flanking from padded extractions).
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
            best = max(counts, key=counts.get)
            if min_freq > 0 and counts[best] / total < min_freq:
                continue
            cons.append(best)
    return ''.join(cons)


# ═══════════════════════════════════════════════════════════════════════
# sear integration
# ═══════════════════════════════════════════════════════════════════════

def setup_sear_workdir(outdir, genome_abs):
    """Shared sear working directory with genome symlink (splits reused)."""
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
    """Run sear in shared workdir. Returns path to output BED, or None."""
    qname  = os.path.basename(query_fa)
    q_dest = os.path.join(sear_wd, qname)
    shutil.copy(query_fa, q_dest)

    cmd = f"THREADS={threads} sear -s {search_n} -k {qname} genome.fa"
    run(cmd, cwd=sear_wd)

    stem = os.path.splitext(qname)[0]
    bed  = os.path.join(sear_wd, f'gen-{stem}.bed')
    bnk  = os.path.join(sear_wd, f'gen-{stem}.bnk')
    if os.path.exists(bnk):
        os.remove(bnk)
    os.remove(q_dest)

    return bed if os.path.exists(bed) and os.path.getsize(bed) > 0 else None


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


def filter_hits(hits, min_hit_len, min_bitscore_frac):
    """Apply noise guards to a sorted hit list.

    Guards:
      1. min_hit_len        — drop hits shorter than this (absolute bp).
      2. min_bitscore_frac  — drop hits below top_bitscore * frac.
         Relative filter adapts to each search; removes the long weak tail.

    Returns filtered list (order preserved).
    """
    if not hits:
        return hits
    bs_floor = hits[0]['bitscore'] * min_bitscore_frac
    return [
        h for h in hits
        if h['length'] >= min_hit_len and h['bitscore'] >= bs_floor
    ]


def consensus_quality(aligned_seqs):
    """Fraction of alignment columns where any base has strict majority (>50%).

    A low value means the flanks are too noisy to build a reliable consensus.
    """
    if not aligned_seqs:
        return 0.0
    maxlen = max(len(s) for s in aligned_seqs)
    n_total = n_confident = 0
    for i in range(maxlen):
        counts = {}
        for s in aligned_seqs:
            b = s[i].upper() if i < len(s) else '-'
            if b in 'ACGT':
                counts[b] = counts.get(b, 0) + 1
        total = sum(counts.values())
        if total == 0:
            continue
        n_total += 1
        if max(counts.values()) / total > 0.5:
            n_confident += 1
    return n_confident / n_total if n_total else 0.0


# ═══════════════════════════════════════════════════════════════════════
# Auto-trim — trim reconstruction to the region supported by genomic hits
# ═══════════════════════════════════════════════════════════════════════

def compute_coverage_profile(aligned_fa):
    """From a MAFFT alignment file, compute per-reconstruction-position coverage.

    Coverage = number of other aligned sequences that have a base (not gap)
    at the same column as each non-gap position of the reconstruction.

    Returns list of int, length == number of bases in the reconstruction.
    """
    entries = read_multi_fasta(aligned_fa)
    recon_aln = None
    hit_alns = []
    for name, seq in entries:
        if name == 'SINE_reconstruction':
            recon_aln = seq
        else:
            hit_alns.append(seq)

    if recon_aln is None or not hit_alns:
        return []

    coverage = []
    for i in range(len(recon_aln)):
        if recon_aln[i].upper() not in 'ACGT':
            continue  # skip gap columns in reconstruction row
        n = sum(1 for s in hit_alns
                if i < len(s) and s[i].upper() in 'ACGT')
        coverage.append(n)
    return coverage


def auto_trim(reconstruction_seq, coverage, min_cov=3):
    """Trim reconstruction to the region supported by >= min_cov hits.

    Scans inward from both ends until coverage meets the threshold.
    Returns (trimmed_seq, trim_5p_bp, trim_3p_bp).
    """
    n = len(coverage)
    if n == 0 or n != len(reconstruction_seq):
        return reconstruction_seq, 0, 0

    left = 0
    while left < n and coverage[left] < min_cov:
        left += 1

    right = n - 1
    while right >= 0 and coverage[right] < min_cov:
        right -= 1

    if left > right:
        return reconstruction_seq, 0, 0  # can't trim

    return reconstruction_seq[left:right + 1], left, n - 1 - right


# ═══════════════════════════════════════════════════════════════════════
# Subfamily detection — iterative bank-depletion (sear-scored)
# ═══════════════════════════════════════════════════════════════════════

def detect_subfamilies(seed_fa_path, genome, outdir, sear_wd, csizes, args):
    """Iterative bank-depletion subfamily detection — lightweight variant.

    All scoring is done by sear (fast); only small MAFFT alignments (~30
    sequences) are run per round.

    Algorithm:
    1. Build bank: sear seed against real genome → pad ±subfam_pad →
       extract → bank.fa, indexed as a mini-genome.
    2. Round 1: sear seed against bank → top ~30 hits → MAFFT align →
       strict majority consensus.  Sear consensus against bank → assign
       members (bitscore ≥ threshold).  Remove from bank.
    3. Round 2+: sear each found consensus against remaining bank; pick
       sequence with lowest max-bitscore (most distant from all found).
       Sear that sequence against bank → top ~30 → MAFFT → consensus.
       Assign members → remove.  Repeat.

    Returns list of dicts [{rank, n_members, consensus_len}, ...] or [].
    """
    sfdir = os.path.join(outdir, 'subfamilies')
    os.makedirs(sfdir, exist_ok=True)

    # ── 1. Build bank: search seed → pad → extract ────────────────
    print("  [subfamilies] building sequence bank ...")
    bed_path = run_sear(sear_wd, seed_fa_path, args.search_hits, args.threads)
    if bed_path is None:
        print("  [subfamilies] no seed hits — skipping")
        return []

    bed_dst = os.path.join(sfdir, 'seed_hits.bed')
    shutil.move(bed_path, bed_dst)

    all_hits = parse_bed7(bed_dst)
    usable = [h for h in all_hits if h['length'] >= args.min_hit_len]
    top = usable[:args.subfam_top_hits]

    if len(top) < args.branch_min * 2:
        print(f"  [subfamilies] only {len(top)} usable hits — need ≥ "
              f"{args.branch_min * 2}")
        return []

    print(f"  [subfamilies] {len(all_hits)} raw → {len(usable)} usable → "
          f"taking top {len(top)}")

    pad = args.subfam_pad
    padded_bed = os.path.join(sfdir, 'padded_hits.bed')
    n_written = 0
    with open(padded_bed, 'w') as fh:
        for i, h in enumerate(top):
            clen = csizes.get(h['chrom'], 10**9)
            s = max(0, h['start'] - pad)
            e = min(h['end'] + pad, clen)
            if e - s < 50:
                continue
            fh.write(f"{h['chrom']}\t{s}\t{e}\thit{i}\t0\t{h['strand']}\n")
            n_written += 1

    bank_fa = os.path.join(sfdir, 'bank.fa')
    run(f"bedtools getfasta -s -nameOnly -fi {genome} -bed {padded_bed}"
        f" > {bank_fa}", cwd=sfdir)

    bank = {name: seq for name, seq in read_multi_fasta(bank_fa)}
    remaining = set(bank.keys())

    print(f"  [subfamilies] bank: {len(bank)} sequences "
          f"(±{pad} bp around each seed hit)")

    # ── 2. Index bank as mini-genome for sear ─────────────────────
    bank_abs = os.path.abspath(bank_fa)
    run(f"samtools faidx {bank_abs}")
    bank_sear_wd = setup_sear_workdir(sfdir, bank_abs)

    def score_bank(query_fa_path):
        """sear query against bank → {seq_name: max_bitscore} for remaining."""
        bed = run_sear(bank_sear_wd, query_fa_path,
                       args.search_hits, args.threads)
        if bed is None:
            return {}
        hits = parse_bed7(bed)
        scores = {}
        for h in hits:
            if h['chrom'] in remaining:
                scores[h['chrom']] = max(
                    scores.get(h['chrom'], 0), h['bitscore'])
        return scores

    # ── 3. Iterative depletion ────────────────────────────────────
    found_cons = []   # [(name, seq)]
    results    = []
    current_query_fa = seed_fa_path   # round 1 uses original seed
    round_num  = 0

    while len(remaining) >= args.branch_min:
        round_num += 1
        rdir = os.path.join(sfdir, f'round_{round_num}')
        os.makedirs(rdir, exist_ok=True)

        print(f"\n  [subfamilies] ── round {round_num} "
              f"({len(remaining)} remaining) ──")

        # ── sear query against bank → select top group ────────────
        scores = score_bank(current_query_fa)
        scored = [(n, scores[n]) for n in remaining
                  if n in scores and scores[n] > 0]
        scored.sort(key=lambda x: x[1], reverse=True)

        top_n = min(args.final_top_hits, len(scored))
        if top_n < args.branch_min:
            print(f"  [subfamilies] round {round_num}: only "
                  f"{top_n} bank hits — done")
            break

        top_names = [n for n, _ in scored[:top_n]]
        print(f"  [subfamilies] {len(scored)} bank hits → "
              f"aligning top {len(top_names)}")

        # ── MAFFT align small group → consensus ───────────────────
        group_fa = os.path.join(rdir, 'group.fa')
        write_multi_fasta(group_fa, [(n, bank[n]) for n in top_names])

        aligned_fa = os.path.join(rdir, 'aligned.fa')
        run(f"mafft --thread {args.threads} --localpair "
            f"--maxiterate 1000 --ep 0.123 --nuc --reorder --quiet "
            f"{group_fa} > {aligned_fa}", cwd=rdir)

        aligned = read_multi_fasta(aligned_fa)
        cons = majority_consensus([s for _, s in aligned],
                                  min_depth=2, min_freq=0.6)

        if len(cons) < 30:
            print(f"  [subfamilies] round {round_num}: consensus only "
                  f"{len(cons)} bp — done")
            break

        cons_name = f'SINE_subfamily_{round_num}'
        cons_fa = os.path.join(rdir, 'consensus.fa')
        write_fasta(cons_fa, cons_name, cons)

        # ── sear consensus against bank → assign members ──────────
        member_scores = score_bank(cons_fa)
        if member_scores:
            top_bs = max(member_scores.values())
            threshold = top_bs * args.subfam_id
            members = {n for n in remaining
                       if member_scores.get(n, 0) >= threshold}
        else:
            members = set(top_names) & remaining

        if len(members) < args.branch_min:
            members = set(top_names) & remaining

        if len(members) < args.branch_min:
            print(f"  [subfamilies] round {round_num}: only "
                  f"{len(members)} members — done")
            break

        # ── Save outputs ──────────────────────────────────────────
        write_fasta(os.path.join(outdir, f'SINE_subfamily_{round_num}.fa'),
                    cons_name, cons)

        remaining -= members
        found_cons.append((cons_name, cons))
        results.append(dict(rank=round_num, n_members=len(members),
                            consensus_len=len(cons)))

        print(f"  [subfamilies] subfamily {round_num}: {len(members)} members "
              f"→ {len(cons)} bp consensus  "
              f"({len(remaining)} remaining)")

        # ── Find most distant remaining for next round ────────────
        if len(remaining) < args.branch_min:
            break

        # sear each found consensus against remaining → per-seq max bitscore
        max_bs = {}   # {name: max bitscore to any found consensus}
        for fc_name, fc_seq in found_cons:
            fc_fa = os.path.join(rdir, f'_score_{fc_name}.fa')
            write_fasta(fc_fa, fc_name, fc_seq)
            fs = score_bank(fc_fa)
            for n in remaining:
                max_bs[n] = max(max_bs.get(n, 0), fs.get(n, 0))

        scorable = {n: s for n, s in max_bs.items() if n in remaining}
        if not scorable:
            break

        centroid = min(scorable, key=scorable.get)
        print(f"  [subfamilies] most distant remaining: "
              f"bitscore {scorable[centroid]:.0f} to nearest found")

        next_fa = os.path.join(sfdir, f'_next_query_{round_num}.fa')
        write_fasta(next_fa, f'centroid_r{round_num}', bank[centroid])
        current_query_fa = next_fa

    # ── 4. Align all found subfamily consensuses ──────────────────
    if len(found_cons) >= 2:
        combined_fa = os.path.join(sfdir, 'all_consensuses.fa')
        write_multi_fasta(combined_fa, found_cons)

        subfam_aln = os.path.join(outdir, 'SINE_subfamilies_aligned.fa')
        run(f"mafft --thread {args.threads} --localpair "
            f"--maxiterate 1000 --ep 0.123 --nuc --reorder --quiet "
            f"{combined_fa} > {subfam_aln}", cwd=sfdir)
        print(f"\n  [subfamilies] subfamily alignment: "
              f"SINE_subfamilies_aligned.fa")

    if not results:
        print("  [subfamilies] no subfamilies detected")
    else:
        print(f"\n  [subfamilies] {len(results)} subfamilies detected")

    return results


# ═══════════════════════════════════════════════════════════════════════
# Hit pool — deduplicated accumulator across steps of one arm
# ═══════════════════════════════════════════════════════════════════════

class HitPool:
    """Accumulate sear hits across steps; keep unique loci sorted by bitscore."""

    def __init__(self):
        self._seen  = set()   # (chrom, start, end, strand)
        self._hits  = []

    def add(self, hits):
        for h in hits:
            key = (h['chrom'], h['start'], h['end'], h['strand'])
            if key not in self._seen:
                self._seen.add(key)
                self._hits.append(h)
        self._hits.sort(key=lambda h: h['bitscore'], reverse=True)

    def top(self, n):
        return self._hits[:n]

    def __len__(self):
        return len(self._hits)


# ═══════════════════════════════════════════════════════════════════════
# Flank extraction
# ═══════════════════════════════════════════════════════════════════════

def write_flank_bed(hits, direction, flank_size, csizes, out_path):
    """Write strand-aware flank BED (extension zone beyond each hit).

    direction='upstream'    → 5′ direction relative to the hit strand
    direction='downstream'  → 3′ direction relative to the hit strand
    Returns number of valid entries written.
    """
    n = 0
    with open(out_path, 'w') as fh:
        for i, h in enumerate(hits):
            clen = csizes.get(h['chrom'], 10**9)
            if direction == 'upstream':
                if h['strand'] == '+':
                    fs, fe = max(0, h['start'] - flank_size), h['start']
                else:
                    fs, fe = h['end'], min(h['end'] + flank_size, clen)
            else:  # downstream
                if h['strand'] == '+':
                    fs, fe = h['end'], min(h['end'] + flank_size, clen)
                else:
                    fs, fe = max(0, h['start'] - flank_size), h['start']
            if fe - fs < 10:
                continue
            fh.write(f"{h['chrom']}\t{fs}\t{fe}\thit{i}\t0\t{h['strand']}\n")
            n += 1
    return n


def write_hit_bed(hits, out_path):
    """Write BED for the full hit intervals (used for fallback extraction)."""
    with open(out_path, 'w') as fh:
        for i, h in enumerate(hits):
            fh.write(f"{h['chrom']}\t{h['start']}\t{h['end']}\t"
                     f"hit{i}\t0\t{h['strand']}\n")


# ═══════════════════════════════════════════════════════════════════════
# Single arm walking step
# ═══════════════════════════════════════════════════════════════════════

def _dump_stats(stats, step_dir):
    with open(os.path.join(step_dir, 'stats.json'), 'w') as fh:
        json.dump(stats, fh, indent=2)


def walk_step(step_num, arm_name, fragment_fa, genome, outdir,
              sear_wd, args, csizes, hit_pool):
    """Execute one step of one arm.

    arm_name: '5prime' or '3prime'
    direction for flank: 'upstream' (5′ arm) or 'downstream' (3′ arm)

    Returns:
      (new_fragment_seq, status)
      new_fragment_seq = None if step fails
      status: 'extended' | 'no_hits' | 'too_few_flanks' | 'short_consensus'
    """
    direction = 'upstream' if arm_name == '5prime' else 'downstream'
    tag  = f"step_{step_num:03d}_{arm_name}"
    sd   = os.path.join(outdir, tag)
    os.makedirs(sd, exist_ok=True)
    shutil.copy(fragment_fa, os.path.join(sd, 'query.fa'))

    stats = dict(step=step_num, arm=arm_name, tag=tag, direction=direction)

    # 1 ── sear search ─────────────────────────────────────────────────
    print(f"  [{tag}] sear -s {args.search_hits} ...")
    bed_src = run_sear(sear_wd, fragment_fa, args.search_hits, args.threads)

    if bed_src is None:
        print(f"  [{tag}] no hits")
        stats['status'] = 'no_hits'
        _dump_stats(stats, sd)
        return None, 'no_hits'

    bed_dst = os.path.join(sd, 'sear_hits.bed')
    shutil.move(bed_src, bed_dst)

    # 2 ── select top hits; apply noise guards; add to pool ──────────
    all_hits      = parse_bed7(bed_dst)
    filtered_hits = filter_hits(all_hits, args.min_hit_len,
                                args.min_bitscore_frac)
    top           = filtered_hits[:args.top_hits]
    hit_pool.add(filtered_hits)

    stats['hits_raw']      = len(all_hits)
    stats['hits_filtered'] = len(filtered_hits)
    stats['hits_used']     = len(top)
    print(f"  [{tag}] {len(all_hits)} raw → {len(filtered_hits)} after "
          f"noise filter → top {len(top)} | pool: {len(hit_pool)}")

    if len(top) < args.branch_min:
        stats['status'] = 'too_few_hits'
        _dump_stats(stats, sd)
        return None, 'too_few_hits'

    # 3 ── extract directional flanks ──────────────────────────────────
    flank_bed = os.path.join(sd, 'flanks.bed')
    nf = write_flank_bed(top, direction, args.fragment, csizes, flank_bed)

    if nf < args.branch_min:
        print(f"  [{tag}] only {nf} usable flanks (need {args.branch_min})")
        stats['status'] = 'too_few_flanks'
        _dump_stats(stats, sd)
        return None, 'too_few_flanks'

    flanks_fa = os.path.join(sd, 'flanks.fa')
    run(f"bedtools getfasta -s -nameOnly -fi {genome} -bed {flank_bed}"
        f" > {flanks_fa}", cwd=sd)

    seqs = read_multi_fasta(flanks_fa)
    if len(seqs) < args.branch_min:
        stats['status'] = 'too_few_flanks'
        _dump_stats(stats, sd)
        return None, 'too_few_flanks'
    stats['flanks'] = len(seqs)

    # 4 ── MAFFT align ─────────────────────────────────────────────────
    aligned_fa = os.path.join(sd, 'aligned.fa')
    print(f"  [{tag}] MAFFT aligning {len(seqs)} flanks ({args.fragment} bp)...")
    run(f"mafft --thread {args.threads} --threadtb {args.threads} "
        f"--localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet "
        f"{flanks_fa} > {aligned_fa}", cwd=sd)

    aligned   = read_multi_fasta(aligned_fa)
    aln_seqs  = [s for _, s in aligned]
    qual      = consensus_quality(aln_seqs)
    cons      = majority_consensus(aln_seqs)

    stats['cons_quality'] = round(qual, 3)

    if qual < args.min_cons_quality:
        print(f"  [{tag}] consensus quality {qual:.2f} < "
              f"{args.min_cons_quality} — flanks too noisy, stop")
        stats['status'] = 'low_quality'
        _dump_stats(stats, sd)
        return None, 'low_quality'

    if len(cons) < 10:
        print(f"  [{tag}] consensus only {len(cons)} bp — stop")
        stats['status'] = 'short_consensus'
        _dump_stats(stats, sd)
        return None, 'short_consensus'

    cons_fa = os.path.join(sd, 'fragment_consensus.fa')
    write_fasta(cons_fa, f'{tag}_cons', cons)

    stats['status']          = 'extended'
    stats['consensus_len']   = len(cons)
    _dump_stats(stats, sd)

    print(f"  [{tag}] → {len(cons)} bp consensus fragment")
    return cons, 'extended'


# ═══════════════════════════════════════════════════════════════════════
# Fallback: align accumulated top hits against current reconstruction
# ═══════════════════════════════════════════════════════════════════════

def do_fallback_alignment(arm_name, hit_pool, reconstruction_seq,
                          genome, outdir, args):
    """Extract top-N accumulated hits, align with current consensus, save."""
    fbdir = os.path.join(outdir, f'fallback_{arm_name}')
    os.makedirs(fbdir, exist_ok=True)

    top = hit_pool.top(args.fallback_hits)
    if not top:
        print(f"  [fallback_{arm_name}] no hits in pool — skipping alignment")
        return None

    print(f"  [fallback_{arm_name}] aligning top {len(top)} accumulated hits "
          f"against current reconstruction ...")

    # Extract hit sequences
    hit_bed = os.path.join(fbdir, 'top_hits.bed')
    write_hit_bed(top, hit_bed)
    hits_fa = os.path.join(fbdir, 'top_hits.fa')
    run(f"bedtools getfasta -s -nameOnly -fi {genome} -bed {hit_bed}"
        f" > {hits_fa}", cwd=fbdir)

    # Write current reconstruction as reference
    recon_fa = os.path.join(fbdir, 'reconstruction_ref.fa')
    write_fasta(recon_fa, f'reconstruction_{arm_name}', reconstruction_seq)

    # Combine ref + hits for alignment
    combined_fa = os.path.join(fbdir, 'combined.fa')
    run(f"cat {recon_fa} {hits_fa} > {combined_fa}", cwd=fbdir)

    # MAFFT align
    aligned_fa = os.path.join(fbdir, 'fallback_aligned.fa')
    run(f"mafft --thread {args.threads} --localpair --maxiterate 1000 "
        f"--ep 0.123 --nuc --reorder --quiet "
        f"{combined_fa} > {aligned_fa}", cwd=fbdir)

    print(f"  [fallback_{arm_name}] alignment: {aligned_fa}")
    return aligned_fa


# ═══════════════════════════════════════════════════════════════════════
# Final alignment of top hits vs completed reconstruction
# ═══════════════════════════════════════════════════════════════════════

def do_final_alignment(reconstruction_seq, genome, outdir,
                       sear_wd, args):
    """Search genome with full reconstruction; extract & align top hits.

    Instead of using the fragment-sized pool hits (which are only ~40 bp each),
    we re-search the genome with the complete reconstruction so sear finds
    full-element-length copies.  Then auto-trim the reconstruction to the
    region actually supported by genomic hits.

    Returns (aligned_fa_path, trimmed_seq_or_None).
    """
    fdir = os.path.join(outdir, 'final_alignment')
    os.makedirs(fdir, exist_ok=True)

    # Write reconstruction as sear query
    recon_fa = os.path.join(fdir, 'reconstruction.fa')
    write_fasta(recon_fa, 'SINE_reconstruction', reconstruction_seq)

    print(f"  [final_alignment] sear search with {len(reconstruction_seq)} bp "
          f"reconstruction ...")
    bed_path = run_sear(sear_wd, recon_fa, args.search_hits, args.threads)

    if bed_path is None:
        print("  [final_alignment] no hits — skipping")
        return None, None

    bed_dst = os.path.join(fdir, 'sear_hits.bed')
    shutil.move(bed_path, bed_dst)

    all_hits  = parse_bed7(bed_dst)
    filtered  = filter_hits(all_hits, args.min_hit_len,
                            args.min_bitscore_frac)
    top_hits  = filtered[:args.final_top_hits]

    if not top_hits:
        print("  [final_alignment] no hits after filtering — skipping")
        return None, None

    print(f"  [final_alignment] {len(all_hits)} raw → {len(filtered)} "
          f"filtered → top {len(top_hits)} for alignment ...")

    hit_bed = os.path.join(fdir, 'top_hits.bed')
    write_hit_bed(top_hits, hit_bed)
    hits_fa = os.path.join(fdir, 'top_hits.fa')
    run(f"bedtools getfasta -s -nameOnly -fi {genome} -bed {hit_bed}"
        f" > {hits_fa}", cwd=fdir)

    combined_fa = os.path.join(fdir, 'combined.fa')
    run(f"cat {recon_fa} {hits_fa} > {combined_fa}", cwd=fdir)

    aligned_fa = os.path.join(fdir, 'final_aligned.fa')
    run(f"mafft --thread {args.threads} --localpair --maxiterate 1000 "
        f"--ep 0.123 --nuc --reorder --quiet "
        f"{combined_fa} > {aligned_fa}", cwd=fdir)

    # Copy top-level outputs
    shutil.copy(hits_fa,    os.path.join(outdir, 'SINE_top_hits.fa'))
    shutil.copy(aligned_fa, os.path.join(outdir, 'SINE_top_hits_aligned.fa'))

    # ── auto-trim ─────────────────────────────────────────────────────
    trimmed_seq = None
    coverage = compute_coverage_profile(aligned_fa)
    if coverage:
        trimmed, t5, t3 = auto_trim(
            reconstruction_seq, coverage, min_cov=args.trim_min_cov)
        if t5 > 0 or t3 > 0:
            trimmed_path = os.path.join(outdir,
                                        'SINE_reconstruction_trimmed.fa')
            write_fasta(trimmed_path, 'SINE_trimmed', trimmed)
            trimmed_seq = trimmed
            print(f"  [auto-trim] {len(reconstruction_seq)} → "
                  f"{len(trimmed)} bp  "
                  f"(cut {t5} bp 5′, {t3} bp 3′)")
        else:
            print("  [auto-trim] full reconstruction supported — "
                  "no trimming needed")

    print(f"  [final_alignment] done: {aligned_fa}")
    return aligned_fa, trimmed_seq


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description='SINE_walker — reconstruct full SINE from a seed fragment')

    ap.add_argument('-s', '--seed', required=True,
                    help='Seed FASTA (SINE fragment; any region of the element)')
    ap.add_argument('-g', '--genome', required=True,
                    help='Genome FASTA (will be indexed if needed)')
    ap.add_argument('-o', '--outdir', required=True,
                    help='Output directory')
    ap.add_argument('-d', '--direction', default='both',
                    choices=('both', '5prime', '3prime'),
                    help='Walk direction(s) (default: both)')

    g = ap.add_argument_group('search tuning')
    g.add_argument('--steps',        type=int,   default=10,
                   help='Max walking steps per arm (default: 10)')
    g.add_argument('--search-hits',  type=int,   default=50,
                   help='sear: stop after N hits per step (default: 50)')
    g.add_argument('--top-hits',     type=int,   default=20,
                   help='Top hits used for flank extraction per step (default: 20)')
    g.add_argument('--fragment',     type=int,   default=40,
                   help='Fragment size in bp for each extension step (default: 40)')
    g.add_argument('--fallback-hits', type=int,  default=10,
                   help='Top accumulated hits for fallback alignment (default: 10)')
    g.add_argument('--final-top-hits', type=int, default=30,
                   help='Top hits used in full-reconstruction final alignment '
                        '(default: 30). More hits → better auto-trim coverage.')

    g = ap.add_argument_group('quality / noise guards')
    g.add_argument('--branch-min',       type=int,   default=3,
                   help='Min sequences per alignment to build consensus (default: 3)')
    g.add_argument('--min-hit-len',      type=int,   default=20,
                   help='Drop sear hits shorter than this bp (default: 20). '
                        'Absolute guard against sub-fragment noise.')
    g.add_argument('--min-bitscore-frac', type=float, default=0.5,
                   help='Drop hits below top_bitscore × frac (default: 0.5). '
                        'Relative filter removes the weak-hit tail adaptively.')
    g.add_argument('--min-cons-quality', type=float, default=0.5,
                   help='Min fraction of alignment columns with strict majority '
                        'vote >50%% (default: 0.5). Low value = noisy flanks → '
                        'step fails rather than extending on garbage.')
    g.add_argument('--trim-min-cov',    type=int,   default=3,
                   help='Auto-trim: min number of hits covering each position '
                        '(default: 3). Positions at the edges with fewer hits '
                        'are trimmed from the reconstruction.')

    g = ap.add_argument_group('subfamily detection')
    g.add_argument('--no-subfam',       action='store_true',
                   help='Skip subfamily detection entirely.')
    g.add_argument('--subfam-id',       type=float, default=0.60,
                   help='Bitscore fraction for member assignment: bank '
                        'sequences scoring ≥ best_match × this value are '
                        'assigned to the subfamily (default: 0.60).')
    g.add_argument('--subfam-pad',      type=int,   default=500,
                   help='Pad each seed hit by this many bp on each side '
                        'to build the sequence bank (default: 500).')
    g.add_argument('--subfam-top-hits', type=int,   default=100,
                   help='Max seed hits used for the sequence bank '
                        '(default: 100).')

    g = ap.add_argument_group('system')
    g.add_argument('--threads',      type=int,   default=4,
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

    csizes    = load_chrom_sizes(genome + '.fai')
    seed_name, seed_seq = read_fasta(args.seed)
    if len(seed_seq) < 10:
        sys.exit("Seed sequence too short (< 10 bp)")

    os.makedirs(args.outdir, exist_ok=True)
    sear_wd = setup_sear_workdir(args.outdir, genome)

    arms = []
    if args.direction in ('both', '5prime'):
        arms.append('5prime')
    if args.direction in ('both', '3prime'):
        arms.append('3prime')

    print(f"SINE_walker")
    print(f"  Seed      : {seed_name} ({len(seed_seq)} bp)")
    print(f"  Genome    : {genome}")
    print(f"  Direction : {args.direction}  arms: {arms}")
    print(f"  Fragment  : {args.fragment} bp  |  steps: {args.steps}  |  "
          f"top-hits: {args.top_hits}  |  branch-min: {args.branch_min}")
    print()

    log = dict(
        seed=seed_name, seed_len=len(seed_seq), genome=genome,
        direction=args.direction, config=vars(args),
        arms={}
    )

    # ── arm state ─────────────────────────────────────────────────────
    # 5′ arm: fragments stored innermost-first (index 0 = closest to seed)
    #   reconstruction prefix = reversed → fragments[-1] + ... + fragments[0]
    # 3′ arm: fragments stored innermost-first
    #   reconstruction suffix = fragments[0] + ... + fragments[-1]
    arm_frags  = {arm: [] for arm in arms}   # list of fragment strings
    arm_query  = {arm: seed_seq for arm in arms}  # current query seq
    arm_pool   = {arm: HitPool() for arm in arms}
    arm_active = {arm: True for arm in arms}
    arm_stop   = {arm: 'running' for arm in arms}

    # ── walking loop (alternates arms each step for fair progression) ──
    for step in range(1, args.steps + 1):
        any_active = False
        for arm in arms:
            if not arm_active[arm]:
                continue
            any_active = True

            # Write current query fragment
            qfa = os.path.join(args.outdir, f'_q_{arm}_s{step}.fa')
            write_fasta(qfa, f'q_{arm}_s{step}', arm_query[arm])

            new_frag, status = walk_step(
                step, arm, qfa, genome, args.outdir,
                sear_wd, args, csizes, arm_pool[arm]
            )
            os.remove(qfa)

            if new_frag is None:
                arm_active[arm] = False
                arm_stop[arm]   = status
                print(f"  [{arm}] stopped: {status}")
            else:
                arm_frags[arm].append(new_frag)
                arm_query[arm] = new_frag    # next step queries this fragment
                total = (sum(len(f) for f in arm_frags[arm]) + len(seed_seq))
                print(f"  [{arm}] +{len(new_frag)} bp  "
                      f"(arm {sum(len(f) for f in arm_frags[arm])} bp, "
                      f"reconstruction ~{total} bp)")
            print()

        if not any_active:
            print("  All arms stopped.")
            break

    for arm in arms:
        if arm_active[arm]:
            arm_stop[arm] = 'max_steps'

    # ── build reconstruction ───────────────────────────────────────────
    # 5′ arm: fragments[0] is closest to seed → print reversed for 5→3 order
    prefix = ''.join(reversed(arm_frags.get('5prime', [])))
    suffix = ''.join(arm_frags.get('3prime', []))
    reconstruction = prefix + seed_seq + suffix

    print()
    print(f"Reconstruction: {len(prefix)} bp 5′ ext  +  "
          f"{len(seed_seq)} bp seed  +  {len(suffix)} bp 3′ ext  "
          f"= {len(reconstruction)} bp total")

    recon_path = os.path.join(args.outdir, 'SINE_reconstruction.fa')
    write_fasta(recon_path, f'SINE_{seed_name}', reconstruction)
    print(f"  Wrote: {recon_path}")

    # ── fallback alignments for stopped arms ──────────────────────────
    for arm in arms:
        # Build partial reconstruction for this arm's perspective
        current_recon = prefix + seed_seq + suffix
        do_fallback_alignment(
            arm, arm_pool[arm], current_recon, genome, args.outdir, args
        )

    # ── final alignment: top hits vs full reconstruction ───────────────
    print()
    _, trimmed_seq = do_final_alignment(
        reconstruction, genome, args.outdir, sear_wd, args
    )

    # ── subfamily detection (re-search with original seed) ────────────
    subfam_results = []
    if not args.no_subfam:
        print()
        subfam_seed_fa = os.path.join(args.outdir, '_subfam_seed.fa')
        write_fasta(subfam_seed_fa, seed_name, seed_seq)
        subfam_results = detect_subfamilies(
            subfam_seed_fa, genome, args.outdir, sear_wd, csizes, args)
        if os.path.exists(subfam_seed_fa):
            os.remove(subfam_seed_fa)

    # ── log ───────────────────────────────────────────────────────────
    for arm in arms:
        log['arms'][arm] = dict(
            steps=len(arm_frags[arm]),
            extension_bp=sum(len(f) for f in arm_frags[arm]),
            stop_reason=arm_stop[arm],
            pool_size=len(arm_pool[arm])
        )
    log['reconstruction_len'] = len(reconstruction)
    log['prefix_len']  = len(prefix)
    log['suffix_len']  = len(suffix)
    if trimmed_seq:
        log['trimmed_len'] = len(trimmed_seq)
    if subfam_results:
        log['subfamilies'] = subfam_results

    log_path = os.path.join(args.outdir, 'walk_log.json')
    with open(log_path, 'w') as fh:
        json.dump(log, fh, indent=2)

    print()
    print('=' * 60)
    print('SINE_walker done')
    print(f'  Reconstruction : {recon_path}  ({len(reconstruction)} bp)')
    if trimmed_seq:
        print(f'  Trimmed        : '
              f'{os.path.join(args.outdir, "SINE_reconstruction_trimmed.fa")}  '
              f'({len(trimmed_seq)} bp)')
    print(f'  Top hits       : {os.path.join(args.outdir, "SINE_top_hits.fa")}')
    print(f'  Aligned hits   : {os.path.join(args.outdir, "SINE_top_hits_aligned.fa")}')
    if subfam_results:
        print(f'  Subfamilies    : {len(subfam_results)} detected')
        for sf in subfam_results:
            print(f'    subfamily {sf["rank"]}: {sf["n_members"]} members, '
                  f'{sf["consensus_len"]} bp')
        if len(subfam_results) >= 2:
            print(f'  Subfam align   : '
                  f'{os.path.join(args.outdir, "SINE_subfamilies_aligned.fa")}')
    print(f'  Log            : {log_path}')
    for arm in arms:
        info = log['arms'][arm]
        print(f'    {arm}: {info["steps"]} steps, '
              f'{info["extension_bp"]} bp extended  ({info["stop_reason"]})')
    print('=' * 60)


if __name__ == '__main__':
    main()
