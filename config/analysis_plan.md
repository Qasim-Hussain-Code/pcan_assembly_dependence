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

Only my own PCAn runs are compared with each other.

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
Peter et al. 2018 assembly if the strain is in that study, otherwise the
short-read assembly with the highest contig N50. Ties go to the lower accession.

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
4 lists two different centromere sequences for that strain
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

## Arm 3 (added 9 October 2026, before any arm 3 data were examined)

No read run had been listed, downloaded or subsampled, and no assembly built
from reads, when this section was committed.

### Question

When the same reads are assembled at known depths with two assemblers, at
what depth do PCAn's calls and CDEII lengths settle, and do the real breaks
resemble either breakpoint model of arm 1?

### Strains and reads

A strain qualifies for arm 3 when all of these hold:

1. it is in the arm 2 pair set, so its long-read assembly is fixed;
2. that long-read assembly is a haploid or collapsed ScRAP assembly;
3. Peter et al. 2018 (Table S1) list it as homozygous and euploid, because in a
   heterozygous strain the assembly from short reads and the long-read
   assembly could differ at a centromere for reasons that have nothing to do
   with depth;
4. it has at least one paired-end Illumina run in ENA study PRJEB13017, linked
   to its code through the reads-file name in Peter et al. 2018 Table S17.

If a strain has more than one qualifying run, the run with the largest base
count is used. Every run accession is verified in ENA and recorded with its
base count, read count, file checksums and the evidence that links it to the
strain.

Strains are processed one at a time in an order fixed now: ascending by the
SHA-256 of `arm3|<code>`. The time budget (48 hours by default) decides how
many are processed; strains are cut from the end of that order, never depths,
seeds or controls. After the first strain, its measured time and disk are
used to project the rest, and the run refuses to start a strain that will not
fit.

### Depths and subsampling

Target depths: 5x, 10x, 20x, 40x, and the full run capped at 80x. The genome
size G is the total length of the strain's long-read assembly. The fraction
for a target depth d is d × G / B, where B is the run's base count (both
mates) from the ENA file report; a fraction of 1 or more means the whole run
is used. Seed 11 is used at every depth, and seeds 22 and 33 are added at 10x
to measure the variance between subsamples. Because seqtk keeps a read when
its seeded random number falls below the fraction, the seed-11 subsamples are
nested, each lower depth inside the next.

Reads are streamed from ENA and subsampled during download with seqtk 1.5
(`seqtk sample -s <seed> - <fraction>`), the same seed for both mates, all
subsamples of one strain drawn in a single pass over each mate's file. After
subsampling, the read names of mate 1 and mate 2 must agree read for read.
Taking the first N reads is not used, because read order follows flow cell
position. The realised depth (bases in the subsample divided by G, before and
after trimming) is reported beside every target. A subsample that misses its
target by more than 10 per cent is reported as it is, never redrawn.

### Trimming and assembly

fastp 1.4.0 with `--detect_adapter_for_pe` and otherwise default settings; its
JSON report is kept.

Every subsample is assembled twice, so that the assembler is a measured
factor: SPAdes 4.3.0 in `--isolate` mode, with the memory cap from
`project.conf` and its temporary directory inside the data directory, and
MEGAHIT 1.2.9 with default k-mers. Contigs are the assembly (`contigs.fasta`,
`final.contigs.fa`). Scaffolds are not used: they insert runs of N whose
lengths are estimated, not observed. No contig is filtered by length. An
assembler run that fails, for example by exceeding the memory cap, is recorded
with its error and not repeated with other settings; that cell is missing,
not imputed.

### Calls and comparison

PCAn with the frozen arm 0 settings, exactly as in arm 2. Each assembly is
compared with the strain's long-read assembly by the arm 2 method: the same 2
kb flanks, placement thresholds, five statuses and null windows. The null
windows of a strain are the same windows as in arm 2, because they are drawn
from the same long-read assembly with the same seeds. The `n_run` status
cannot occur, because contigs carry no N.

### Confirmatory outcomes

D1. Recall relative to the long-read calls: the fraction of a strain's
long-read centromeres that are `intact_called`, per assembler and depth.
Reported against realised depth as the median and interquartile range over
strains, and as the mean over strains with a percentile bootstrap interval
over strains (10,000 resamples, seed 20261009).

D2. CDEII agreement: among `intact_called` centromeres, the fraction whose
CDEII length equals the long-read call's, per assembler and depth, with the
same summaries.

D3. Settling depth, per assembler, for D1 and for D2: the lowest target depth
at which the mean over strains is within 0.05 of its value at full depth and
stays within 0.05 at every higher target depth. The value 0.05 is arbitrary.
This is a summary of D1 and D2 and carries no test.

D4. Break-rate ratio against the AT-matched null at each target depth (seed
11) and for each assembler, as in arm 2 C4: effect size with a bootstrap
interval over strains, and a stratified permutation test.

D5. Which arm 1 breakpoint model the real breaks resemble, as in arm 2 C5,
for each assembler at full depth. The same analysis at the lower depths is
exploratory.

Descriptive, reported for every assembly: contig N50, L50 and contig count, so
that arm 3 can be placed on the arm 1 axis; and, at 10x, the range of D1 and D2
across the three seeds for each strain and assembler.

The Holm family of confirmatory tests becomes: arm 2 C4 and C5, arm 3 D4 (five
depths for each of two assemblers) and arm 3 D5 (two assemblers), fourteen
tests in all.

### Controls

Positive control for the liftover: in arms 2 and 3, each long-read assembly is
compared with itself; every centromere must come out `intact_called`.
Negative control: the null windows. Technical replicates: the three 10x seeds.

### Sensitivity

The arm 2 alternatives apply to arm 3: flank length 1, 2 and 5 kb, placement
identity 0.90, 0.95 and 0.99, null AT bin width 2.5, 5 and 10 percentage
points.

## Changes recorded on 9 October 2026, after the arm 2 results and before any arm 3 data

Written after the arm 2 confirmatory results had been read and before any arm
3 reads were downloaded. The text above is unchanged.

Arm 3 downloads overlap assembly. Strains are still assembled one at a time,
in the order fixed above, but while one strain assembles, the next strain's
reads download. Done in turn, the download leaves the processors idle and the
assembly leaves the link idle, and the time budget then covers fewer
strains. The depths, seeds, assemblers, controls and the rule that cuts
strains from the end of the order are unchanged. The projection before each
strain now uses the longest download and the longest assembly seen so far,
because a strain's elapsed time no longer adds up from the two.

Two exploratory additions to arm 2, outside the confirmatory outcomes and
labelled exploratory wherever they appear. First, the kept flanks are aligned
again and every full-length match is counted, secondary ones included, to
tell a flank that matches two places (as both haplotypes of a heterozygous
strain can in a short-read assembly) from one that matches nowhere; the
registered placement rule leaves both unplaced. Second, for each centromere
that is intact and called, whether the short-read contig extends 10 kb beyond
both ends of the call, the arm 1 definition of synteny-checkable.

Two points of implementation, which change no arm 2 result. In the
permutation test of C4 and D4, a permutation in which one group has no break
now counts as at least as extreme as the observed split, where it was
skipped before; no arm 2 permutation had that property. The Holm adjustment
always uses the fourteen registered tests, and a test not yet run enters as
p = 1.

## Change recorded on 9 October 2026, before any arm 3 read was assembled

The first arm 3 strain, CCN, could not be downloaded as a stream. Each of the
two read files broke after about 40 minutes, and one file broke twice more
after 16 and 7 minutes, the last time with a TLS record that failed its
integrity check; a stream that breaks must start again from its first byte,
and all three attempts for that file failed. No read reached an assembler.

Reads are therefore downloaded to a temporary file, resuming each broken
transfer from its last verified byte, and checked against ENA's MD5. Every
depth and seed is then drawn in one pass over the file, with the same seqtk
command and seeds as before, and the file is deleted. The subsamples are the
ones streaming would have drawn: on the test run used for the pipeline, the
fourteen assemblies built both ways have identical contig N50s. A strain
whose reads still cannot be downloaded intact is logged and left for a later
run instead of being counted as processed. CCN is processed again from the
start, first in the order as before.

## Change recorded on 9 October 2026, after the third arm 3 strain

The first three strains, CCN, AEG and BAM, were assembled with SPAdes on 16
threads under the 4 GB cap from `project.conf`. SPAdes failed at 40x and at
full depth for all three, and at 20x for CCN, with "mmap(2) failed: Cannot
allocate memory" while it held about 1.5 GB resident. SPAdes' `-m` limits its
address space, not its resident memory, and each thread reserves address
space of its own; on this machine, whose Linux side has 6 GB, the cap could
not be raised far. On BAM's 40x reads, SPAdes on 4 threads under the same 4
GB cap finished at 1.5 GB resident, and its contigs were identical, sequence
for sequence, to those of a 16-thread run with no cap.

From the fourth strain on, SPAdes runs with one thread per GB of its cap,
here 4. The cap itself, the assembler version, its mode and every other
setting are unchanged, and MEGAHIT keeps all 16 threads. As registered above,
the failed runs of the first three strains are recorded with their errors and
not repeated; those cells stay missing. The thread count of every assembler
run is recorded in each strain's manifest from the fourth strain on.

## Edit recorded on 10 October 2026

Three passages above were shortened to remove background remarks: one
sentence after the arm 2 question, a parenthesis in the pair rule, and three
words in the ploidy section. No rule, test, threshold or result changed.

## Change recorded on 10 October 2026, after the ninth complete arm 3 strain

Inside WSL, read files downloaded by curl kept arriving at full size with a
wrong MD5, and with a different wrong MD5 on every try: both files of CCN
once, one file of ADS once, one file of BCN on all three tries and one file of
ANE on its first two. A file that is wrong at the source would fail with the
same MD5 each time. The machine's disk gave the same MD5 on five reads of a
1.5 GB file, with the page cache dropped before each, and on a copy of it. The
second file of ANE, downloaded by Windows' own curl.exe on the same machine,
matched ENA's MD5 at the first try (1,062 s), and two copies of it into WSL
matched as well.

From BCN on, under WSL, `scripts/11_assemble_reads.sh` downloads each read
file with Windows' curl.exe into a directory under the Windows temporary
directory, copies it into the data directory and checks the copy against
ENA's MD5. Everything after that is unchanged. The stage was stopped at 16:28
UTC during ANE's third try and restarted with this change; BCN and ANE have no
done marker, so the restarted run takes them first, in the registered order.
Depths, seeds, assemblers and the 48 h budget are unchanged. Two periods of
download that no strain row recorded, AGK from 18:42 to 19:50 UTC on 9 October
(until the Linux VM shut down) and ANE from 14:18 to 16:28 UTC on 10 October,
were added to `logs/arm3_strains.tsv` with the status `stopped`, so that the
budget counts them.

## Correction recorded on 10 October 2026, after BCN's second run

The conclusion of the previous section was wrong. Downloaded by Windows'
curl.exe, BCN's files came back with a wrong MD5 in five of six transfers,
again a different one each time. Six downloads of the first 256 MiB of BCN's
first file, three by Windows' curl.exe and three by curl inside WSL, gave five
identical copies; the sixth differed from them in 201 consecutive bytes, which
held an error message from ENA's server written in place of the data:
`<Error><Code>ConnectionClosedException</Code><Message>Premature end of
Content-Length delimited message body (expected: 1,438,572,164; received:
27,444,504)</Message><ErrorMessage/><RequestId/></Error>`. Every other byte
matched. The damage comes from the server, not from this machine, and
downloading on Windows does not avoid it.

From the restart at 20:11 UTC, read files are downloaded by curl inside WSL
again, as before the previous change. After each complete download, every
such message in the file is located, and the bytes from 4 MiB before it to 4
MiB after it are fetched again as a byte range and written back in place; the
whole file is then checked against ENA's MD5 as before, and a mismatch still
starts the file again, at most twice. On the damaged 256 MiB copy this gave
the MD5 of the five intact copies in 19 s, and it found no message in an
intact copy. The stage was stopped at 20:10 UTC during ANE's first try; that
period is in `logs/arm3_strains.tsv` as `stopped`. BCN's second run is logged
as `download failed`. BCN and ANE have no done marker, so the restarted run
takes them first, in the registered order. Depths, seeds, assemblers and the
48 h budget are unchanged.

## Change recorded on 11 October 2026, after the fourteenth arm 3 strain

The arm 3 time budget is raised from 48 h to 96 h so that all 33 eligible
strains are assembled. When the change was made, at 02:19 UTC, 14 strains
were complete, about 29 h of the budget had been used, and no strain had been
cut. The running stage read its budget when it started and keeps 48 h: it
will refuse the first strain that would end past that, and its log will say
so. A second run, started automatically when the first ends, reads the new
budget and assembles the remaining strains in the registered order. Depths,
seeds and assemblers are unchanged.

