# Arm 1: design of the controlled perturbations (a simulation)

Written and committed before any assembly was fragmented or edited. Arm 1
alters real assemblies at chosen positions. Its truth is PCAn's own call on
the unaltered assembly, so every measure here is self-consistency, written as
"recall relative to the unfragmented call" or "relative to the unedited
call", and never as accuracy.

## 1a. Which assemblies are fragmented

An arm 0 assembly is eligible when all of these hold:

1. its species has the role `primary` in `config/species_arm0.tsv`;
2. NCBI gives its assembly level as Chromosome or Complete Genome;
3. its arm 0 PCAn run finished;
4. PCAn's call count equals the number of nuclear chromosome sequences, counted
   from the NCBI sequence report as assembled molecules of type Chromosome;
5. every call lies on a different nuclear chromosome sequence.

From the eligible assemblies, at most one per genus is taken: the one with the
highest contig N50 as computed in `results/arm0/assembly_stats.tsv` (sequences
split at every run of N), ties broken by the lower accession in ASCII order.
The realised set and its size are written to `results/arm1/genomes.tsv`.

## 1a. How they are fragmented

Two breakpoint models, both clean cuts between adjacent bases with no
sequence removed (`scripts/lib/fragment.py`):

- `uniform`: a homogeneous Poisson process along each sequence.
- `at_weighted`: an inhomogeneous Poisson process whose intensity in each 500 bp
  window is proportional to exp(gamma × (a − a_genome)), with a the window's
  AT fraction, a_genome the assembly's, and gamma = ln(10) / 0.30, so that a
  window 30 percentage points more AT-rich than the genome is ten times as
  likely to break. The weights are normalised so that both models place the
  same expected number of cuts at the same rate. The functional form and the
  constant are assumptions; arms 2 and 3 measure where real assemblies break.

Target contig N50: 1,000, 500, 200, 100, 50, 20, 10 and 5 kb. For each genome,
model and level, the cut rate is calibrated by bisection so that the mean
realised contig N50 of 12 calibration draws is within 2 per cent of the
target. Calibration uses its own seed. Ten replicates are then drawn per
genome, model and level, each with a seed derived from the accession, model,
level and replicate number (the first 8 hexadecimal digits of the SHA-256 of
`accession|model|level|replicate`, read as an integer). Realised contig N50,
L50 and contig count are recorded for every replicate.

A level at or above the assembly's own contig N50 needs no cut. Its replicates
are recorded as unperturbed, are not run through PCAn (arm 0 shows PCAn is
deterministic on an unchanged input), and count as full recall.

Each replicate is generated, run through `scripts/04_run_pcan.sh`, scored and
deleted. The call table, the candidate table from the filter trace, the cut
positions and the fragment statistics are kept.

## 1a. Outcomes

For every call PCAn makes on the unfragmented assembly (the truth set), in
every replicate:

- `retained_exact`: a call at the same source coordinates with the same CDEII
  length;
- `retained_altered`: an overlapping call that is not exact;
- `lost_geometric`: no overlapping call, and a cut predicted the loss;
- `lost_pipeline`: no overlapping call, and no cut predicted it; the filter
  trace names the step that removed the candidate, or records that no
  candidate was formed.

The geometric prediction follows PCAn's code, not the paper's description. A
forward-strand call is predicted lost when a cut falls inside its extraction
window, the CDEIII motif plus 249 bp upstream (275 bp), because a window that
starts before the start of a sequence comes back empty. A reverse-strand call
is predicted lost when a cut falls between its CDEI and the end of its CDEIII
motif; a cut further upstream only shortens its window. Whether any cut falls
in the window, for either strand, is recorded as well.

Calls in a replicate that overlap no truth call are `gained`.

Recall relative to the unfragmented call is the fraction of truth calls that
are retained, exact or altered. Exact recall counts only `retained_exact`. The
synteny-checkable fraction is the fraction of truth calls whose fragment
extends 10 kb beyond both ends of the call, the span PCAn's synteny step reads
(`Range = 10000`), computed whether or not the call was retained.

Summaries are the median and interquartile range across replicates, per
genome and pooled over genomes, per model and level. Intervals on pooled
quantities resample genomes (10,000 resamples, seed 20261009). Nothing is
bootstrapped over centromeres.

## 1a. Sensitivity of the AT-weighted model

gamma is re-run at ln(5) / 0.30 and ln(20) / 0.30, at target N50 of 50, 20 and
10 kb, with five replicates per genome and level.

## 1b. Which genomes are edited

S288C (GCA_000146045.2), plus, from the arm 1a eligible assemblies, the one
whose PCAn calls have the shortest median CDEII length and the one whose calls
have the longest, ties broken by the lower accession, skipping S288C if it is
itself the shortest or the longest. The realised set is written to
`results/arm1/variant_genomes.tsv`.

## 1b. Edits

`scripts/lib/variants.py`. An insertion duplicates the s bases of CDEII that
end at the insertion point and inserts the copy immediately after them, a
tandem duplication with both copies inside CDEII. A deletion removes s
contiguous bases from inside CDEII. CDEI and the CDEIII motif, as PCAn
reported them on the unedited genome, are never touched. Positions are drawn
uniformly from the valid range with a seed derived as in 1a from
`accession|contig|size`.

Single variants: every centromere PCAn calls on each genome, one at a time,
at -10, -5, +5, +10, +15, +20, +30 and +40 bp.

Progressive transitions, S288C only: 1, 2, 4, 8, 12 and all 16 centromeres
carry a +10 bp insertion, the centromeres chosen at random without
replacement, five replicates per level with seeds derived from `k|replicate`.

## 1b. Outcomes

For each edited centromere: whether PCAn still calls it (a call overlapping
its expected span), and whether the reported CDEII length equals the planted
length (the unedited length plus s). For the progressive design: the total
call count, the number of edited centromeres called, and the number called at
the planted length. Per genome, fractions carry Clopper-Pearson intervals over
centromeres, labelled as within-genome; pooled fractions resample the three
genomes, which gives only a rough interval and is described as such.

## Budget

If the measured time per PCAn run projects arm 1 beyond the time budget in
`project.conf`, replicates are cut before genomes, never below five per level.
What was cut is written to `results/arm1/cuts.tsv` and stated in the README.
