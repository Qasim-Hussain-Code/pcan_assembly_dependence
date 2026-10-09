# Analysis plan for arms 2 and 3

Registered by commit before any arm 2 data were examined: no arm 2 assembly
had been downloaded, run through PCAn or aligned when this file was
committed, and the list of strain pairs had not been built. The commit
timestamp is the registration. Any later change is added below in a dated
section; the text above it stays as it is.

## Arm 2 question

For *S. cerevisiae* strains with both a published short-read and a published
long-read assembly: do PCAn's calls on the two agree, on centromere count and
on each CDEII length? Where they disagree, is the region missing, broken, or
present but uncalled? Do short-read assemblies break near centromeres more
often than in random sequence of the same AT content?

Only my own PCAn runs are compared with each other. The paper's published
calls play no part in arm 2.

## Pairs

The independent unit is the strain, with one pair per strain.

A strain qualifies when one long-read assembly and one short-read assembly of
it exist and both are public. Long-read means the assembly was built from
PacBio or Oxford Nanopore reads, with or without short-read polishing.
Short-read means it was built from Illumina reads alone. Assemblies built from
Sanger or 454 reads, or whose technology the source does not state, are not
used.

Strains are matched by identifier only: the same strain name or code in the
metadata of both sources, or a documented synonym published by a source (for
example the standardised names of Peter et al. 2018). The evidence for every
match is written beside the pair in `config/arm2_pairs.tsv`. A name that is
merely similar is not a match.

If a strain has several candidate assemblies of one kind, the long-read
assembly with the highest contig N50 is used, and for the short-read side the
Peter et al. 2018 assembly if the strain is in that study (it is the
short-read source the PCAn paper itself used), otherwise the short-read
assembly with the highest contig N50. Ties go to the lower accession.

Every candidate strain that is dropped is listed with its reason.

## PCAn

PCAn v1.0 as frozen in arm 0: genus *Saccharomyces*, CDEIII motif
`CDEIII_Saccharomyces_MEME.txt` at 1e-6, CDEI motif `CDEI_ZT_MEME.txt` at 1e-2,
through `scripts/04_run_pcan.sh`. Nothing is changed for any strain or
technology.

## Locating each long-read centromere in the short-read assembly

For every PCAn call on the long-read assembly (the element runs from the
start of CDEI to the end of the CDEIII motif), the left flank is the 2,000 bp
immediately before the element and the right flank the 2,000 bp immediately
after it. A call with less than 2,000 bp of sequence on either side on its
long-read contig cannot be lifted over and is reported as such, not scored.

Both flanks are aligned to the short-read assembly with minimap2 2.31
(`-x asm10 -c`). A flank is placed only if exactly one primary alignment
reaches identity at least 0.95 (matching bases over alignment block length),
query coverage at least 0.90 and mapping quality at least 20. Otherwise it is
unplaced. Nothing is forced.

Each long-read centromere then receives one status (`scripts/lib/status.py`):

- `intact_called`: both flanks on one short-read sequence, in order and
  orientation, at most the element length plus 1,000 bp apart and overlapping
  by at most 50 bp, no N between them, and a short-read PCAn call within 10 bp
  of the region between them;
- `intact_uncalled`: the same, with no short-read call there;
- `n_run`: the same, with an N between the flanks;
- `split`: one flank placed, or both placed but on different sequences, in
  opposite orientations, or further apart than allowed;
- `absent`: neither flank placed.

The `n_run` status cannot occur in a contig assembly, which carries no N.

The CDEII length difference is the short-read call's CDEII length minus the
long-read call's, for `intact_called` centromeres. Short-read calls that
correspond to no long-read centromere are counted per strain.

A break near a centromere is any of `split`, `n_run` or `absent`: the short-read
assembly does not hold the element and its 2 kb flanks as one gap-free piece.

## Ploidy and heterozygosity

The ploidy and zygosity that each source reports are recorded per strain. In a
heterozygous strain a collapsed short-read assembly and a phased or
haploid-resolved long-read assembly can disagree at a centromere without
either being wrong. Every disagreement in a heterozygous strain is reported in
its own column and is not called an assembly error. Where Supplementary Data
4 of the PCAn paper lists two different centromere sequences for that strain
and chromosome, the disagreement is labelled as possibly allelic. Excluding
heterozygous strains is one of the sensitivity analyses.

## Null model

To ask whether short-read assemblies break near centromeres more often than
elsewhere, null windows are drawn from each long-read assembly. Each
centromere gets 20 null windows of the same length, with AT content within
plus or minus 2.5 percentage points of the centromere element's (a bin 5
points wide), containing no N, at least 10 kb from any PCAn call and at least
2 kb from a contig end. Positions are drawn uniformly without replacement
among those that qualify, with a seed derived from the strain and centromere
(the first 8 hexadecimal digits of the SHA-256 of
`strain|contig|start|null`). If fewer than 20 qualify, all are used and the
shortfall is reported.

Null windows go through exactly the same flank extraction, alignment,
placement thresholds and status function as the centromeres. A null that
skipped the liftover would measure the liftover.

## Confirmatory outcomes

C1. Centromere count. Per strain, PCAn's call count on each assembly. Reported:
the fraction of strains with equal counts, with a Wilson interval over
strains.

C2. Status. The fraction of long-read centromeres in each of the five
statuses, pooled, with a percentile bootstrap over strains (10,000
resamples, seed 20261009).

C3. CDEII length. Among `intact_called` centromeres, the fraction whose CDEII
length is identical in the two assemblies, with a bootstrap over strains.

C4. Breaks. The break-rate ratio R = (broken centromeres / centromeres) /
(broken null windows / null windows), pooled over strains. Effect size: R with
a percentile bootstrap interval over strains (10,000 resamples, seed
20261009). Test: a stratified permutation test of log R, permuting the
centromere and null labels within each strain, 10,000 permutations, two-sided.
The difference in break proportions, centromere minus null, is reported with
its own bootstrap interval.

C5. Which arm 1 breakpoint model the real breaks resemble. Each short-read
contig is aligned to its strain's long-read assembly (minimap2 `-x asm5`,
primary alignments of mapping quality at least 20 and at least 1 kb). Each
short-read contig end that is not at a long-read contig end gives one
breakpoint position on the long-read assembly. The long-read assembly is
tiled in 500 bp windows. Under the uniform model a breakpoint falls in a
window with probability proportional to the window's length; under the
AT-weighted model with probability proportional to its length times
exp(gamma × (a − a_genome)). The log-likelihood of the observed breakpoints is
computed under the uniform model (gamma = 0), under arm 1's AT-weighted model
(gamma = ln(10) / 0.30), and at the maximum-likelihood gamma. Decision rule:
the real breaks resemble the AT-weighted model more if arm 1's AT-weighted
model has the higher log-likelihood of the two and the 95 per cent bootstrap
interval of the maximum-likelihood gamma (over strains, 10,000 resamples)
excludes zero; the uniform model more if the uniform model has the higher
log-likelihood and the interval includes zero; otherwise neither is favoured,
and the report says so. The likelihood-ratio test of gamma = 0 against free
gamma (chi-squared, one degree of freedom) gives the p-value.

The loss decomposition of arm 1 is applied to these real breaks: a long-read
call missing from the short-read assembly is a geometric loss if its status is
`split`, `n_run` or `absent`, and a pipeline loss if it is `intact_uncalled`.
Short-read calls with no long-read counterpart are gains.

## Statistics

Proportions carry Wilson intervals, or Clopper-Pearson where a count is below
10. Differences and ratios carry percentile bootstrap intervals over strains,
10,000 seeded resamples. No interval resamples centromeres. P-values appear
only beside effect sizes and are Holm-corrected across the family of
confirmatory tests: C4 and C5 in arm 2, and the tests that the arm 3 extension
of this plan adds. There are no pre-registered subgroups, so no subgroup
conclusion is drawn. Everything else is exploratory and labelled so.

## Sensitivity analyses

Each is re-run in full and reported in `results/sensitivity/`; the README says
for each confirmatory conclusion whether any alternative changes it.

- Flank length: 1, 2 (primary) and 5 kb.
- Placement identity: 0.90, 0.95 (primary) and 0.99.
- Null AT bin width: 2.5, 5 (primary) and 10 percentage points.
- Heterozygous strains excluded.

## Clarification of the pair rule, 9 October 2026

Added after registration and before any arm 2 data were examined: no arm 2
assembly had been downloaded, run through PCAn or aligned, and the pair list
had not been built. The text above is unchanged.

Candidate long-read assemblies are those deposited in INSDC under an assembly
accession and version that NCBI Datasets returns as current. Two kinds of
public long-read assembly are therefore not candidates:

- the haplotype-1 assemblies of the ScRAP panel (O'Donnell et al. 2023), which
  NCBI Datasets does not serve;
- the assemblies of Loegler et al. 2025, which exist only inside a single
  16.8 GB archive on Zenodo and were scaffolded against the reference genome
  with ragout. The archive alone would take half of the disk budget, and
  reference-guided scaffolding would make the long-read side partly a copy of
  the reference.

Strains are matched through the standardized names (three-letter codes) that
Peter et al. 2018 (Supplementary Table S1) and O'Donnell et al. 2023
(Supplementary Table 1) both use, and the isolate names in the two tables
must also agree. A ScRAP strain listed under another strain name, such as a
spore or colony derivative, is not matched. The short-read assembly of a
Peter et al. 2018 strain is taken from the single archive that the 1002 Yeast
Genomes project distributes; the archive member for each code is identified
by listing the archive, not by assuming a file name.
