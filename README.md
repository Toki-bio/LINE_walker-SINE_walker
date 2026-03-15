# LINE_walker & SINE_walker

Iterative seed-extension tools for reconstructing full transposable element sequences from a single genomic fragment.

Both tools use the same core approach:
**sear** (Smith-Waterman / minimap2 search) → **bedtools** (flanking sequence extraction) → **mafft** (multiple alignment) → **majority-vote consensus** → repeat.

---

## Tools

### LINE_walker

Reconstructs a full LINE element starting from a known SINE tail (or any LINE fragment). Walks in one direction (default: upstream = toward LINE 5′ end) following the LINE body step by step. Clusters flanks with vsearch at each step to detect variant branches; by default pursues at most 3 LINE variants (`--max-variants 3`).

```
LINE_walker -s sine_tail.fa -g genome.fa -o outdir/ [options]
```

### SINE_walker

Reconstructs a full SINE element from any internal seed. Walks **bidirectionally** (5′ and 3′ arms simultaneously, independent searches). Each arm extends by short fragments (default 40 bp). When no significant hits are found, fallback alignment of accumulated hits vs the current reconstruction is produced for manual inspection.

```
SINE_walker -s seed.fa -g genome.fa -o outdir/ [options]
```

---

## Dependencies

| Tool | Purpose |
|---|---|
| `sear` | Fast Smith-Waterman genome search (part of SINEderella) |
| `minimap2` | Long-query genome search (auto-used by sear for queries ≥1000 bp) |
| `samtools` | Genome indexing (`faidx`) |
| `bedtools` | Flanking sequence extraction (`getfasta`) |
| `mafft` | Multiple sequence alignment |
| `vsearch` | Clustering (LINE_walker only) |
| `python3` | ≥ 3.8 |

Activate the appropriate conda environment before running.

---

## sear — minimap2 mode

When the query sequence is **≥ 1000 bp**, sear automatically switches from ssearch36 (fragment-scanning) to **minimap2** for faster, whole-genome alignment. This is particularly important for LINE_walker, where accumulated queries grow beyond what ssearch36 handles efficiently.

**Sensitivity cascade:** sear tries minimap2 presets in order — `asm5` → `asm10` → `asm20` — stopping at the first preset that finds hits. This ensures the most exact matches are preferred while still finding divergent copies.

**`-b` auto-sizing:** When `-b` is used without an explicit N:
- Query < 5000 bp → `-b 10` (10 best hits aligned)
- Query ≥ 5000 bp → `-b 5` (5 best hits aligned)

An explicit value (e.g. `-b 20`) always overrides the auto-sizing.

**Override:** Use `--force-ssearch` to disable minimap2 and force ssearch36 for all query lengths.

---

## LINE_walker options

```
required:
  -s / --seed          Seed FASTA (SINE tail or LINE fragment)
  -g / --genome        Genome FASTA (auto-indexed if needed)
  -o / --outdir        Output directory

search tuning:
  --steps        N     Max walking steps (default: 30)
  --search-hits  N     sear: stop after N hits per step (default: 50)
  --top-hits     N     Keep N best hits by bitscore (default: 20)
  --flank        N     Flank extraction size in bp (default: 150)
  --seed-window  N     Use last N bp of accumulated seq as query (default: 200)
  --direction          upstream | downstream | left | right (default: upstream)
  --max-jump     N     Continuity filter: max distance from previous-step
                       loci for fresh sear hits (default: 5000)

extended-flank optimization:
  --extended-flank N   Flank size in bp when reusing previous-step hit
                       coordinates before fallback genome search (default: 500)
  --try-extended-first Attempt extended-flank extraction before genome search
                       (default: enabled)
  --no-try-extended-first
                       Disable extended-flank fast path; always run genome
                       search each step (original behavior)

clustering:
  --cluster-id   F     vsearch cluster identity (default: 0.80)
  --branch-min   N     Min seqs per cluster to continue (default: 3)
  --rescue-min-members N
                       If no cluster reaches --branch-min, rescue only the
                       single largest cluster when it has at least N members
                       (default: 2)
  --max-variants N     Max LINE variants (branches) to pursue (default: 3).
                       Extra clusters are discarded by ascending member count.

system:
  --threads      N     Threads for mafft / sear (default: 4)
```

**Output:**
- `LINE_candidates.fa` — full accumulated sequences (seed + all extensions)
- `LINE_extensions.fa` — novel extension only
- `walk_log.json` — per-step stats
- `step_NNN/` — per-step intermediate files

**Extended-flank fast path (LINE_walker):**
- By default, each step first tries reusing the **previous step's hits** and
  extracting larger directional flanks (`--extended-flank`, default 500 bp).
- If that path fails (too few usable flanks/consensus), LINE_walker
  automatically falls back to a fresh genome search using the normal
  `--flank` size (default 150 bp).

**Integrity behavior (LINE_walker):**
- No "all-flanks as single group" fallback is used anymore.
- If no cluster reaches `--branch-min`, only the **single largest cluster** can
  be rescued (controlled by `--rescue-min-members`), which preserves integrity
  better than mixed-flank consensus.
- After fallback fresh `sear`, hits are continuity-filtered against the
  previous-step loci (`--max-jump`) to reduce branch jumping across unrelated
  genomic copies.
- `upstream/downstream` are strand-aware (TE-oriented), while `left/right`
  are absolute genomic-coordinate directions.

---

## SINE_walker options

```
required:
  -s / --seed          Seed FASTA (any SINE fragment)
  -g / --genome        Genome FASTA (auto-indexed if needed)
  -o / --outdir        Output directory

direction:
  -d / --direction     both | 5prime | 3prime (default: both)

search tuning:
  --steps        N     Max steps per arm (default: 10)
  --search-hits  N     sear: stop after N hits per step (default: 50)
  --top-hits     N     Top hits used for flank extraction (default: 20)
  --fragment     N     Extension step size in bp (default: 40)
  --fallback-hits N    Accumulated hits used in fallback / final alignment (default: 10)

noise guards:
  --branch-min       N    Min sequences to build consensus; fewer → stop (default: 3)
  --min-hit-len      N    Drop sear hits shorter than N bp (default: 20)
  --min-bitscore-frac F   Drop hits below top_bitscore × F (default: 0.5)
  --min-cons-quality  F   Stop if < F fraction of alignment columns have majority vote (default: 0.5)

system:
  --threads      N     Threads for mafft / sear (default: 4)
```

**Output:**
- `SINE_reconstruction.fa` — `[5′ extension] + [seed] + [3′ extension]`
- `SINE_top_hits.fa` — top accumulated genomic hits (both arms)
- `SINE_top_hits_aligned.fa` — mafft alignment of top hits + reconstruction
- `fallback_5prime/` / `fallback_3prime/` — alignment at stop point, per arm
- `walk_log.json` — per-step stats including consensus quality scores
- `step_NNN_5prime/` / `step_NNN_3prime/` — per-step intermediate files

---

## Noise guards (SINE_walker)

Short seeds (30–50 bp) can produce many noisy sear hits. Three guards are applied each step before building a consensus:

1. **`--min-hit-len`** — absolute length cutoff; drops sub-fragment matches
2. **`--min-bitscore-frac`** — relative to the best hit in that step; removes the weak-score tail adaptively
3. **`--min-cons-quality`** — after MAFFT alignment, if fewer than this fraction of columns have a clear majority base, the flanks are too inconsistent and the step stops cleanly

To recover from over-aggressive filtering, lower `--min-bitscore-frac` (e.g. 0.3) or `--min-cons-quality` (e.g. 0.4). To reduce extension on noise, raise them or increase `--branch-min`.

---

## Example

```bash
# Reconstruct a SINE from a known internal fragment
SINE_walker \
  -s my_sine_seed.fa \
  -g species_genome.fa \
  -o sine_walk_out/ \
  --fragment 40 \
  --steps 15 \
  --threads 8

# Reconstruct a LINE from a SINE tail (walk toward the LINE 5′ end)
LINE_walker \
  -s sine_tail.fa \
  -g species_genome.fa \
  -o line_walk_out/ \
  --direction upstream \
  --flank 150 \
  --steps 30 \
  --max-variants 3 \
  --threads 8
```

---

## How it works

```
SINE_walker (bidirectional):

         ←  5′ arm extends ←          → 3′ arm extends →
[step N] [step 2] [step 1] [ SEED ] [step 1] [step 2] [step N]
                             sear
                              ↓
                        top hits (filtered)
                              ↓
                   extract N-bp flanks (bedtools)
                              ↓
                        mafft alignment
                              ↓
                   quality check → majority consensus
                              ↓
                        new fragment = query for next step
```

Each arm maintains its own accumulated **hit pool** across all steps. When an arm stops, the pool's top-N hits are aligned against the current reconstruction as a fallback inspection output.
