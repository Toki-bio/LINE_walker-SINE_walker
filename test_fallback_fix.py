#!/usr/bin/env python3
"""
Test: walk_step continues past stale continuity filter after extended-flank fallback.

Scenario reproduced from the Nephrurus_levis run:
  - Steps 2-7 use extended flanks (no sear, prev_hits_bed = step-1 sear BED on bank_seq0)
  - Step 8: extended flanks fail → sear fallback
  - Sear returns hits on bank_seq1 (different bank fragment)
  - OLD CODE: continuity filter compares bank_seq1 hits against bank_seq0 prev_hits → 0 kept → STOP
  - NEW CODE: fell_back_from_extended=True → filter skipped → walk continues
"""
import os
import sys
import argparse
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import line_walker as lw


def make_args(**kwargs):
    defaults = dict(
        branch_min=5, rescue_min_members=5, cluster_hits=50,
        min_bitscore_frac=0.90, top_hits=20, flank=150,
        anchor_overlap=60,
        extended_flank=500, try_extended_first=True, direction='5prime',
        search_hits=50, threads=4, cluster_id=0.80, max_jump=5000,
        dedup_locus_window=30, steps=30, max_variants=3,
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def make_run_mock(seq_count=50, seq_len=150, cluster_size=20):
    """Return a mock for lw.run that fakes mafft, vsearch, and bedtools."""
    def mock_run(cmd, cwd=None):
        if 'mafft' in cmd and '>' in cmd:
            out = cmd.split('>')[-1].strip()
            with open(out, 'w') as f:
                for i in range(seq_count):
                    f.write(f'>hit{i}\n' + 'A' * seq_len + '\n')
        elif 'vsearch' in cmd and '--clusters' in cmd:
            pfx = cmd.split('--clusters')[1].strip().split()[0]
            with open(pfx + '0', 'w') as f:
                for i in range(cluster_size):
                    f.write(f'>hit{i}\n' + 'A' * seq_len + '\n')
        elif 'bedtools' in cmd and '>' in cmd:
            out = cmd.split('>')[-1].strip()
            with open(out, 'w') as f:
                for i in range(seq_count):
                    f.write(f'>hit{i}\n' + 'A' * seq_len + '\n')
        return ''
    return mock_run


class TestExtendedFallbackContinuityFix(unittest.TestCase):

    def test_trim_anchor_overlap_5prime(self):
        trimmed = lw.trim_anchor_overlap('A' * 150 + 'C' * 60, '5prime', 60)
        self.assertEqual(trimmed, 'A' * 150)

    def test_trim_anchor_overlap_3prime(self):
        trimmed = lw.trim_anchor_overlap('C' * 60 + 'G' * 150, '3prime', 60)
        self.assertEqual(trimmed, 'G' * 150)

    def test_write_flank_bed_includes_anchor_overlap(self):
        hits = [dict(chrom='chr1', start=1000, end=1140, strand='+',
                     homology=0.9, length=140, bitscore=100.0)]
        csizes = {'chr1': 100000}

        with tempfile.TemporaryDirectory() as tmpdir:
            bed = os.path.join(tmpdir, 'flanks.bed')
            count = lw.write_flank_bed(
                hits, '5prime', 150, csizes, bed, anchor_overlap=60
            )
            self.assertEqual(count, 1)
            with open(bed) as fh:
                line = fh.readline().rstrip('\n')

        self.assertEqual(line, 'chr1\t850\t1060\thit0\t0\t+')

    # ── helper: write two sets of hits on different bank chromosomes ──────────

    def _write_prev_bed(self, path, chrom='bank_seq0'):
        """Simulate step-1 sear hits (the stale hits in prev_hits_bed)."""
        lw.write_bed7([
            dict(chrom=chrom, start=1000 + i * 100, end=1200 + i * 100,
                 strand='+', homology=0.9, length=200, bitscore=float(100 - i))
            for i in range(10)
        ], path)

    def _write_sear_bed(self, path, chrom='bank_seq1'):
        """Simulate fresh sear hits on a DIFFERENT bank chromosome."""
        lw.write_bed7([
            dict(chrom=chrom, start=5000 + i * 200, end=5150 + i * 200,
                 strand='+', homology=0.9, length=150,
                 bitscore=float(100 - i * 0.1))
            for i in range(50)
        ], path)

    # ── Test 1: demonstrate the root cause ───────────────────────────────────

    def test_cross_chrom_filter_returns_empty(self):
        """
        filter_hits_by_previous_loci kills all hits when chromosomes differ.
        This is what caused the stall: step-1 hits on bank_seq0 vs fresh
        sear hits on bank_seq1 → empty → walk stops.
        """
        prev = [dict(chrom='bank_seq0', start=1000, end=1200, strand='+',
                     homology=0.9, length=200, bitscore=100.0)] * 5
        new  = [dict(chrom='bank_seq1', start=5000, end=5200, strand='+',
                     homology=0.9, length=150, bitscore=100.0)] * 10
        result = lw.filter_hits_by_previous_loci(new, prev, max_jump=5000)
        self.assertEqual(result, [],
            "Cross-chrom hits are entirely filtered — this is the stall root cause")

    # ── Test 2: the fix — walk continues after extended fallback ─────────────

    def test_walk_step_continues_after_extended_fallback(self):
        """
        After fix: fell_back_from_extended=True bypasses the stale continuity
        filter and walk_step produces extensions from fresh sear hits.
        """
        args = make_args()
        csizes = {f'bank_seq{i}': 100000 for i in range(10)}

        with tempfile.TemporaryDirectory() as tmpdir:
            qfa = os.path.join(tmpdir, 'query.fa')
            lw.write_fasta(qfa, 'query', 'A' * 200)

            prev_bed = os.path.join(tmpdir, 'prev_sear_hits.bed')
            self._write_prev_bed(prev_bed, chrom='bank_seq0')

            sear_bed = os.path.join(tmpdir, 'sear_out.bed')
            self._write_sear_bed(sear_bed, chrom='bank_seq1')  # different chrom!

            with patch('line_walker._try_extended_flanks', return_value=([], None)), \
                 patch('line_walker.run_sear', return_value=sear_bed), \
                 patch('line_walker.run', side_effect=make_run_mock()):

                exts, bed_out, used_ext = lw.walk_step(
                    step_num=8,
                    query_fa=qfa,
                    genome='/fake/genome.fa',
                    outdir=tmpdir,
                    sear_wd=tmpdir,
                    args=args,
                    csizes=csizes,
                    branch='',
                    prev_hits_bed=prev_bed,
                    ext_offset=3000,
                )

        self.assertFalse(used_ext,
            "sear fallback path must be taken (not extended)")
        self.assertIsNotNone(bed_out,
            "bed_out must not be None")
        self.assertGreater(len(exts), 0,
            "walk_step must produce extensions — continuity filter must have been skipped")

    # ── Test 3: normal sear path still applies the continuity filter ─────────

    def test_normal_sear_path_still_filters(self):
        """
        Without extended-flank mode, a normal sear step on a different bank
        chromosome should still be caught by the continuity filter and stop.
        This verifies the fix doesn't disable the filter in the normal case.
        """
        args = make_args(try_extended_first=False)   # no extended path
        csizes = {f'bank_seq{i}': 100000 for i in range(10)}

        with tempfile.TemporaryDirectory() as tmpdir:
            qfa = os.path.join(tmpdir, 'query.fa')
            lw.write_fasta(qfa, 'query', 'A' * 200)

            prev_bed = os.path.join(tmpdir, 'prev_sear_hits.bed')
            self._write_prev_bed(prev_bed, chrom='bank_seq0')

            sear_bed = os.path.join(tmpdir, 'sear_out.bed')
            self._write_sear_bed(sear_bed, chrom='bank_seq1')

            with patch('line_walker.run_sear', return_value=sear_bed):
                exts, bed_out, used_ext = lw.walk_step(
                    step_num=8,
                    query_fa=qfa,
                    genome='/fake/genome.fa',
                    outdir=tmpdir,
                    sear_wd=tmpdir,
                    args=args,
                    csizes=csizes,
                    branch='',
                    prev_hits_bed=prev_bed,
                    ext_offset=0,
                )

        self.assertEqual(exts, [],
            "Normal sear path must still be blocked by continuity filter for cross-chrom hits")
        self.assertFalse(used_ext)

    # ── Test 4: same-chrom hits still pass continuity filter normally ─────────

    def test_same_chrom_hits_pass_continuity_filter(self):
        """
        Sear hits on the same chrom as prev_hits, within max_jump, must still
        pass the continuity filter and produce extensions.
        """
        args = make_args(try_extended_first=False)
        csizes = {'bank_seq0': 100000}

        with tempfile.TemporaryDirectory() as tmpdir:
            qfa = os.path.join(tmpdir, 'query.fa')
            lw.write_fasta(qfa, 'query', 'A' * 200)

            prev_bed = os.path.join(tmpdir, 'prev_sear_hits.bed')
            self._write_prev_bed(prev_bed, chrom='bank_seq0')

            # Sear hits on SAME chrom, close to prev hits (within max_jump=5000)
            sear_bed = os.path.join(tmpdir, 'sear_out.bed')
            lw.write_bed7([
                dict(chrom='bank_seq0', start=800 + i * 100, end=1000 + i * 100,
                     strand='+', homology=0.9, length=200,
                     bitscore=float(100 - i * 0.1))
                for i in range(50)
            ], sear_bed)

            with patch('line_walker.run_sear', return_value=sear_bed), \
                 patch('line_walker.run', side_effect=make_run_mock()):

                exts, bed_out, used_ext = lw.walk_step(
                    step_num=5,
                    query_fa=qfa,
                    genome='/fake/genome.fa',
                    outdir=tmpdir,
                    sear_wd=tmpdir,
                    args=args,
                    csizes=csizes,
                    branch='',
                    prev_hits_bed=prev_bed,
                    ext_offset=0,
                )

        self.assertGreater(len(exts), 0,
            "Same-chrom hits within max_jump must pass the continuity filter")


if __name__ == '__main__':
    unittest.main(verbosity=2)
