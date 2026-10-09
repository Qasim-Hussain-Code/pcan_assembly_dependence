# Arm 0: how the reproduction is judged

Written and committed before any arm 0 assembly was run through PCAn. Arm 0
compares my runs of PCAn v1.0 with the calls published in Supplementary Data 2
of Helsen et al. 2026. Agreement with that table is agreement with the
authors' output, not with biology. The only external truth in arm 0 is the
S288C reference against SGD.

## Species and assemblies

`config/species_arm0.tsv` assigns every row of Supplementary Data 5 a role,
with its reason.

- `primary`: the 130 Saccharomycetaceae species whose published calls carry a
  CDEI and a CDEII length, the output the released two-motif pipeline
  produces. These form the reproduction set.
- `cdeiii_only`: three *Yueomyces* species and *Grigorovia jiainica*, whose
  published calls have no CDEI. They are run and their CDEIII loci compared,
  outside the reproduction fraction.
- `no_published_calls`: *Arxiozyma slooffiae*. Run and reported.
- `excluded`: the three *Naumovozyma* species (published with an adapted
  single-motif PCAn that the repository does not provide) and the 28 rows that
  are not Saccharomycetaceae.

Assemblies are downloaded by the exact accession and version in Supplementary
Data 5. A version NCBI marks as replaced or suppressed is not used, nothing is
substituted for it, and its species leaves the reproduction set with that
reason.

PCAn runs with the genus (and, for *Kazachstania*, the species) named in
Supplementary Data 5, which selects the motifs and thresholds from PCAn's own
table. No threshold is changed for any species.

## Matching a published call to a reproduced call

Supplementary Data 2 gives, for each call, the contig and the full centromere
sequence from the start of CDEI to the end of the CDEIII motif, but no
coordinates. Each published call is located by searching for its sequence,
ignoring case, on the named contig and on both strands. Contig names are
translated between GenBank and RefSeq with the NCBI sequence report, because
20 species are published under RefSeq names. A published sequence that occurs
zero times, or more than once, is not forced onto a position; it is reported
as unlocated with that reason.

Each published call then gets one status, checked in this order:

1. `exact`: a reproduced call on the same contig with the same start, end and
   sequence (ignoring case), and therefore the same CDEII length.
2. `same_locus_other_boundaries`: a reproduced call overlapping the published
   interval, with the same CDEII length, that is not exact.
3. `different_cdeii_length`: a reproduced call overlapping the published
   interval with a different CDEII length.
4. `missing`: no reproduced call overlaps the published interval.
5. `unlocated`: the published sequence could not be placed on the assembly.

A reproduced call that overlaps no published call is `extra`.

For every call that is not `exact`, a reason is attached from the filter
trace (`scripts/lib/pcan_trace.py`): whether the CDEIII hit was among FIMO's
first-pass hits, whether a CDEI hit was found in its window, and, if a
candidate was formed, which of PCAn's filters removed it.

## The gate for arms 1 to 3

Primary measure: the fraction of published calls in the reproduction set that
are `exact`.

Interval: percentile bootstrap over species, 10,000 resamples, seed 20261009.
Species are the independent unit; the calls within one species share one
assembly, one motif pair and one set of thresholds, so they are not
independent of one another.

Threshold for proceeding: 0.95, on the point estimate. The value is arbitrary.
It is re-evaluated at 0.90 and 0.99 in `results/sensitivity/`, and the README
says whether either changes the decision. If the point estimate falls below
0.95, arms 1 to 3 do not start until the cause has been found and reported,
checking software version first and input differences second.

Secondary measures, all reported: per species, published and reproduced call
counts; the fraction of species whose calls are all exact (Wilson interval,
species as the unit); the CDEIII-only species compared on CDEIII loci.

## Controls in arm 0

Positive control: S288C (GCA_000146045.2) against the centromere features of
SGD release R64-5-1. SGD lists sixteen; the control passes if PCAn calls a
centromere overlapping each of them and nothing else. Offsets of both ends and
the CDEII length against SGD's CDEII feature are reported. Before the
comparison, the NCBI sequence is checked to be identical to SGD's reference
sequence, chromosome by chromosome.

Negative control: determinism. S288C is run five times; the five call tables
and the five pairs of FIMO outputs must be identical byte for byte.

FIMO: MEME suite 4.11.2, the version PCAn pins and the paper used, installed
from bioconda. PCAn does not set `--max-stored-scores`, so FIMO's default cap
of 100,000 stored matches applies (`fimo` 4.11.2 prints this default in its
usage text). The number of first-pass hits per assembly is recorded, with
whether FIMO reported reaching the cap. PCAn does not set `--bgfile` either,
so FIMO 4.11.2 uses its built-in background frequencies and a hit's p-value
depends only on the bases it covers, not on the rest of the assembly.
