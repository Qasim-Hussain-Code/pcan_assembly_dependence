# Arm 0 gate: the decision to proceed

Written on 2026-10-09 (UTC) to record a decision taken earlier the same day,
after the arm 0 results were committed (05:01 UTC) and before the first arm 1
stage started (06:32 UTC, `logs/resources.tsv`). `config/arm0_reproduction.md`
is unchanged. `run_all.sh` reads the last line of this file.

The gate was missed. PCAn v1.0 reproduced 1,356 of 1,465 published calls
exactly, a fraction of 0.926 (95 per cent interval 0.897 to 0.949, bootstrap
over 129 species), below the threshold of 0.95 (`results/arm0/gate.tsv`). The
protocol says that arms 1 to 3 do not start until the cause has been found
and reported, checking software version first and input differences second.

Software version. The pinned commit, a5fa46f, is the PCAn v1.0 release
archived on Zenodo. Its pipeline code is identical to the GitHub v1.0 tag and
to the repository's current head; only README files differ. FIMO is 4.11.2,
the version the paper used. Rescoring every candidate with the AT formula of
the published CDEII column, and replaying PCAn's filters, leaves 934 published
calls standing, against 1,356 with v1.0's own formula
(`results/arm0/diagnostic_at_formula.tsv`). The published calls were made with
v1.0's scoring.

Input differences. One species, *Kluyveromyces aestuarii*, was published on
the assembly version before the one Supplementary Data 5 lists; on that
version all 8 of its calls are reproduced exactly
(`results/arm0/diagnostic_replaced_versions.tsv`). Five more calls could not
be placed because of the table itself: a published sequence absent from the
assembly (2), repeated in it (2), or on a contig whose name differs (1).

What remains matches the authors' own expectations. PCAn ships a file,
`ExpectedFalsePosFalseNegs.txt`, that lists for each species how many
published calls v1.0 is expected to miss. For 113 of the 129 species, the
number of published calls this installation did not recover equals that
number (`results/arm0/authors_expected_vs_observed.tsv`). Summed over the
species with more misses than that, there are 27 extra misses, and 14 of them
are in one species, *Tetrapisispora taiwanensis*, where false candidates set
the length anchor and 13 published calls fell outside it. The published table
therefore holds calls that v1.0, by its authors' account, does not produce,
and every non-exact call has a recorded cause
(`results/arm0/non_exact_causes.tsv`).

Every later arm compares PCAn's runs with each other: on an assembly and its
altered copy in arm 1, and on two assemblies of one strain in arms 2 and 3.
None is compared with the published table. The shortfall measures how far
this installation reproduces the published output, which is reported as it
is; it does not enter the later comparisons.

decision: proceed with arms 1 to 3
