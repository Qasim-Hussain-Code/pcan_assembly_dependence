# Results tables

Every table here is tab-separated, with one header row, written by a script in
`scripts/`. Coordinates are 1-based and inclusive. Accessions carry their
version. "PCAn" always means PCAn v1.0 (commit a5fa46f) with the patch in
`patches/`, run through `scripts/04_run_pcan.sh`. A "call" is a predicted
centromere, not a validated one.

## Arm 0, reproduction (`results/arm0/`, written by `05_reproduce.py` unless noted)

### `gate.tsv`

The pre-registered gate (`config/arm0_reproduction.md`), one row.

- `measure`: what is estimated.
- `unit_of_replication`: the unit the interval resamples (species).
- `n_species`, `n_calls`, `n_exact`: species in the reproduction set, their published calls, and how many were reproduced exactly.
- `estimate`, `ci_low`, `ci_high`: fraction exact and its 95 per cent percentile bootstrap interval.
- `interval`: how the interval was computed (resamples, seed).
- `threshold`, `decision`: the pre-registered threshold and whether the point estimate reaches it (`proceed` or `stop`).
- `species_all_exact`, `species_all_exact_fraction`, `species_all_exact_wilson_low`, `species_all_exact_wilson_high`: species whose every published call is exact, with a Wilson interval over species.

### `per_centromere.tsv`

One row per published call of the reproduction set and of the CDEIII-only species.

- `species`, `accession`: as in Supplementary Data 5.
- `published_contig`: the contig name Supplementary Data 2 gives; `contig`: the same sequence's GenBank name in the downloaded assembly.
- `chr_no`: the centromere number Supplementary Data 2 gives.
- `published_cdeii_len`: CDEII length in Supplementary Data 2.
- `start`, `end`, `strand`: where the published sequence lies in the assembly (empty if it could not be placed).
- `status`: `exact`, `same_locus_other_boundaries`, `different_cdeii_length`, `missing` or `unlocated`; for CDEIII-only species `same_cdeiii`, `overlapping_other_cdeiii` or `missing`.
- `reproduced_start`, `reproduced_end`, `reproduced_cdeii_len`: the overlapping reproduced call, if any.
- `reason`: for any call that is not exact, what the filter trace shows.

### `non_exact_causes.tsv`, `non_exact_cause_counts.tsv`

Every published call of the reproduction set that is not exact, with one cause; and the count of calls per cause and status.

- `cause`: the category (for example, the CDEIII motif scores below the released threshold, or a named PCAn filter removed the candidate).
- `detail`: the evidence behind the category.
- `median_top5_cdeii_len`, `published_median_cdeii_len`: the CDEII length anchor of PCAn's run (median of the five best candidates) and the median published CDEII length of the species; the cause "length anchor set by false candidates" is assigned when they differ by more than 30.

### `authors_expected_vs_observed.tsv`

Per species of the reproduction set: published calls, reproduced calls and each status count, beside the false negatives (`authors_false_neg`) and false positives (`authors_false_pos`) that the authors' own `ExpectedFalsePosFalseNegs.txt` gives for PCAn v1.0, and whether the missing and extra counts equal them.

### `diagnostic_replaced_versions.tsv`

PCAn on an assembly version that Supplementary Data 5 does not list, run only to find the cause of a discrepancy, outside the reproduction fraction (`config/arm0_diagnostics.tsv` says which and why). `exact` counts published calls reproduced exactly on that version.

### `species_counts.tsv`

One row per Saccharomycetaceae species in Supplementary Data 5: its arm 0 role and reason, the PCAn genus used, published and reproduced call counts, PCAn's run status, and the count of each per-centromere status and of extra calls.

### `reproduced_calls.tsv`

Every call PCAn made on every arm 0 assembly.

- `contig`, `start`, `end`: the call, from the start of CDEI to the end of the CDEIII motif.
- `contig_hit`: PCAn's identifier of the CDEIII hit (the contig name and the hit's row in FIMO's first pass).
- `strand`: the strand of the CDEIII hit, from FIMO's first pass (PCAn does not report it).
- `cdei_start`, `cdei_end`, `cdeiii_start`, `cdeiii_end`: the two motifs, in forward coordinates.
- `cdeiii_pvalue`: FIMO's p-value for the CDEIII hit.
- `cdeii_len`, `cdeii_at`: CDEII length and AT percentage as PCAn computes them (the AT percentage divides by the length plus one).
- `fimo_i_score`, `fimo_ii_score`, `overall_score`: PCAn's CDEIII and CDEI motif scores and its combined score.
- `sequence`: the call's sequence, in the orientation PCAn reports, case as in the assembly.

### `extra_calls.tsv`

Reproduced calls that overlap no published call, with the columns of `reproduced_calls.tsv`.

### `s288c_vs_sgd.tsv`

The positive control: each SGD centromere of S288C (release R64-5-1) beside PCAn's call.

- `sgd_start`, `sgd_end`, `sgd_strand`, `sgd_cdei_len`, `sgd_cdeii_len`, `sgd_cdeiii_len`: SGD's feature and its three elements.
- `called`: whether a PCAn call overlaps it; `pcan_start`, `pcan_end`, `pcan_strand`, `pcan_cdeii_len`: that call.
- `start_offset`, `end_offset`: PCAn's start and end minus SGD's.
- `cdeii_len_difference`: PCAn's CDEII length minus SGD's.
- `cdeii_cdeiii_boundary_offset`: where PCAn places the start of the CDEIII motif minus where SGD places the start of CDEIII (forward strand), or the same at the other end (reverse strand).

### `s288c_sequence_check.tsv`

Each S288C chromosome in NCBI's GCA_000146045.2 against SGD's R64-5-1 reference sequence: lengths and whether the sequences are identical.

### `determinism.tsv`

The negative control: five PCAn runs on S288C. SHA-256 of each run's call table, FIMO passes and candidate table, the run time, and whether all five are identical.

### `fimo_cap.tsv`

Per assembly, FIMO's first-pass (CDEIII) and second-pass (CDEI) hit counts, whether FIMO reported reaching its stored-score cap, and the cap (100,000).

### `diagnostic_at_formula.tsv` (written by `scripts/lib/at_formula_check.py`, run from `05_reproduce.py`)

Whether a different CDEII AT formula would recover published calls that PCAn v1.0 misses. Each assembly's candidate table is rescored three ways (`formula`, described in `formula_description`) and PCAn's filters are replayed inside PCAn's own environment. Per assembly and formula, and summed over assemblies (`accession` = `all`): the published calls, how many of them are still standing after the filters, and how many calls stand in all.

### `assembly_stats.tsv` (written by `03_fetch_assemblies.sh`)

Per cached assembly: NCBI assembly level, sequence count, total length, sequence N50, contig N50 and L50 and contig count (sequences split at every run of N), N runs and N bases, nuclear chromosome and mitochondrial sequences (from the NCBI sequence report), NCBI's own contig and scaffold N50, and the genome's AT fraction.

### `excluded_assemblies.tsv` (written by `03_fetch_assemblies.sh`)

Accession versions not used, with NCBI's status, the current version, and the reason.

### `source_concordance.tsv` (written by `02_fetch_tables.sh`)

The publisher's and Figshare's copies of Supplementary Data 2, 5 and 6 compared cell by cell.

## Arm 1, controlled perturbation, a simulation (`results/arm1/`)

Truth in arm 1 is PCAn's call on the unaltered assembly. Every measure here is
relative to that call: self-consistency, not accuracy.

### `eligibility.tsv`, `genomes.tsv` (written by `06_fragment.py`)

Every chromosome-level or complete assembly of a reproduction-set species, with the eligibility rule of `config/arm1_design.md` applied (`eligible`, `reason`), and the one genome per genus chosen from them (`chosen`; `genomes.tsv` holds only those). Columns: accession, species, PCAn genus (and *Kazachstania* species), NCBI assembly level, contig N50, total length, genome AT fraction, nuclear chromosome sequences, and PCAn's call count.

### `calibration.tsv` (written by `06_fragment.py`)

Per genome, breakpoint model and target contig N50: the cut rate (cuts per base) found by bisection, the model's gamma (empty for the uniform model) and the mean contig N50 of the 12 calibration draws at that rate. A rate of 0 means the assembly is already at or below the target.

### `cuts.tsv` (written by `06_fragment.py`)

Replicates cut from the plan to fit the time budget, if any, with the projected hours before the cut.

### `fragmentation_runs.tsv` (written by `06_fragment.py`)

One row per replicate: `replicate_id`, genome, `model`, `gamma`, `target_n50`, `replicate`, `seed`, `rate_per_bp`, `n_cuts`, the realised contig N50, L50 and contig count, PCAn's `status` (`ok`, `failed`, or `unperturbed` when no cut was needed), `n_calls`, run time, peak memory, and whether the filter trace reproduced PCAn's calls.

### `fragmentation_calls.tsv` (written by `08_score_perturbations.py`)

One row per unfragmented call per replicate.

- `contig`, `start`, `end`, `strand`: the call on the unfragmented assembly.
- `status`: `retained_exact`, `retained_altered`, `lost_geometric` or `lost_pipeline`.
- `window_cut`: a cut falls inside PCAn's extraction window (the CDEIII motif and 249 bp upstream on the forward strand, 250 bp on the reverse).
- `geometric_loss_predicted`: the cuts alone remove the call (forward strand: a cut in the window; reverse strand: a cut between CDEI and the end of the CDEIII motif).
- `synteny_checkable`: the fragment holding the call extends 10 kb beyond both of its ends.
- `fragment_length`: length of that fragment.
- `trace_step`: for a call lost with no cut to explain it, the PCAn filter that removed its candidate, or `no candidate formed`.
- `reported_cdeii_len`, `truth_cdeii_len`: CDEII length in the replicate and on the unfragmented assembly.

### `fragmentation_gained.tsv` (written by `08_score_perturbations.py`)

Calls in a replicate that overlap no unfragmented call: where they lie on the unfragmented assembly (`source`, `source_start`, `source_end`), their CDEII length and AT, the length of the fragment they are on, and what happened to the same candidate in the unfragmented run (`unfragmented_fate`, the filter step that removed it there, and its rank).

### `fragmentation_replicates.tsv` (written by `08_score_perturbations.py`)

Per replicate: unfragmented calls (`n_truth`), counts retained, retained exactly, lost by geometry, lost with the window intact, predicted lost, with a cut window, and synteny-checkable; `recall` (retained over unfragmented calls), `exact_recall`, `geometric_expectation` (one minus the fraction predicted lost), `synteny_checkable_fraction`, `gained`, realised L50 and contig count.

### `recall_by_genome.tsv`, `recall_pooled.tsv` (written by `08_score_perturbations.py`)

Per genome, model and level, and pooled over genomes per model and level: median and interquartile range of realised N50 and of recall over replicates, exact recall, geometric expectation, synteny-checkable fraction, and loss and gain totals. The pooled table adds the mean over genomes of each genome's mean recall and of its geometric expectation, each with a percentile bootstrap interval over genomes (10,000 resamples, seed 20261009).

### `pipeline_losses_by_filter.tsv`, `gained_calls_by_origin.tsv` (written by `08_score_perturbations.py`)

Calls lost with no cut in their window, counted by the filter step that removed them; and gained calls counted by what happened to them in the unfragmented run.

### `variant_genomes.tsv`, `variant_edits.tsv`, `variant_runs.tsv` (written by `07_plant_variants.py`)

The three edited genomes and the rule that chose each; every edit (genome, centromere, strand, size, position, the bases duplicated or deleted, the expected call and CDEII length after the edit; for the progressive design the number of edited centromeres `k` and the replicate); and every PCAn run on an edited genome.

### `planted_single.tsv`, `detection_envelope.tsv`, `planted_progressive.tsv` (written by `08_score_perturbations.py`)

For each single edit: whether PCAn still calls the centromere, the CDEII length it reports, whether that is the planted length or the unedited one, and for a centromere not called, the filter step. The detection envelope gives, per genome and pooled, the fraction called and the fraction called at the planted length for each size: Clopper-Pearson intervals over centromeres within a genome, and a bootstrap over the three genomes for the pooled rows. For the progressive design: per replicate, the calls in the genome, the edited centromeres called and called at +10 bp, and the unedited centromeres left unchanged.

## Arm 2, published assembly pairs (`results/arm2/`, written by `10_compare_pairs.py` unless noted)

Each strain has two assemblies: its long-read assembly from ScRAP (O'Donnell
et al. 2023), which is the reference, and the short-read assembly that Peter et
al. 2018 published. PCAn's calls on the long-read assembly are the elements.
Each element is looked for in the short-read assembly through its two flanks.
The method and the thresholds are registered in `config/analysis_plan.md`.

### `pair_candidates.tsv` (written by `09_build_pairs.py`)

Every strain in both collections: the ScRAP and Peter et al. names, the long-read assemblies available for it (`long_read_options`), whether it entered the pair set (`decision`) and why (`reason`). The pairs used are in `config/arm2_pairs.tsv`.

### `assembly_stats.tsv`, `excluded_assemblies.tsv` (written by `03_fetch_assemblies.sh`)

The long-read assemblies, with the columns of the arm 0 tables of the same names.

### `short_read_assembly_stats.tsv`

Per strain, the published short-read assembly: sequence count, total length, sequence N50, contig N50, L50 and count (sequences split at every run of N, the rule used in arms 0, 1 and 3), N runs and N bases.

### `pcan_runs.tsv`

One row per PCAn run, two per strain (`assembly`: `long` or `short`): call count, run time, peak memory, whether the filter trace reproduced PCAn's calls, and FIMO's first-pass hit count.

### `call_counts.tsv`

Per strain: Peter et al.'s isolate name, the long-read accession, PCAn's call count on each assembly, and the zygosity and ploidy that Peter et al. 2018 (Table S1) give. Outcome C1.

### `centromere_status.tsv`

One row per long-read call (an element). Outcomes C2 and C3.

- `element_id`: `cen|<contig>|<start>`; `kind`: `centromere`.
- `contig`, `start`, `end`, `length`, `at`: the call on the long-read assembly (start of CDEI to end of the CDEIII motif) and its AT fraction.
- `liftable`: both 2 kb flanks fit on the long-read contig; `cdeii_len`: the long-read call's CDEII length.
- `status`: `intact_called`, `intact_uncalled`, `split`, `n_run` or `absent` (defined in `scripts/lib/status.py` and the analysis plan), or `not_liftable`. `split`, `n_run` and `absent` count as broken.
- `detail`: the evidence for the status: where the flanks were placed, or why they were not.
- `target`, `region_start`, `region_end`: the short-read sequence holding the element and the region between its two placed flanks.
- `call_start`, `call_end`, `query_cdeii_len`: the short-read call overlapping that region, for `intact_called`, and its CDEII length.
- `cdeii_difference`: that CDEII length minus the long-read call's.
- `zygosity`: as in `call_counts.tsv`.
- `possibly_allelic`: `yes` when Supplementary Data 4 of the PCAn paper lists two different centromere sequences for this strain at the chromosome whose listed sequence matches either assembly's call; `no` when it lists one, or the same one twice; `chromosome not identified` when no listed sequence matches either call; `no Supplementary Data 4 row` when the strain is not in that table.

### `null_status.tsv`

The AT-matched null windows: 20 per liftable centromere, of the same length, with AT within 2.5 percentage points, at least 10 kb from any call and 2 kb from a contig end, drawn with a seed derived from the strain and centromere. Each goes through the same flanks, placement and status function as the centromeres. `for_centromere` is the `element_id` of the centromere the window was matched to; the other columns are as in `centromere_status.tsv`. Outcome C4.

### `liftover_control.tsv`

The positive control: each long-read assembly compared with itself by the same method, with the columns of `centromere_status.tsv`. Every centromere should come out `intact_called`.

### `short_calls_without_long_counterpart.tsv`

Per strain: the short-read calls, and how many of them are not the `intact_called` match of any long-read centromere.

### `breakpoint_summary.tsv`

Per strain: the short-read contig ends located on the long-read assembly (`breakpoints`), from the short-read contigs' alignments to it (minimap2 `-x asm5`, primary, mapping quality at least 20, at least 1 kb; an end is located when its outermost alignment reaches within 100 bp of it, and ends within 1 kb of a long-read contig end are dropped); the ends that could not be located; the mean AT fraction of the 500 bp long-read windows holding a breakpoint; and the genome's length-weighted mean window AT. The input to outcome C5.

### `exploratory_centromere_flanks.tsv`, `exploratory_flank_summary.tsv` (exploratory, not in the analysis plan)

Why a broken element is broken. The flanks of the primary analysis are aligned again exactly as in the comparison, and every alignment that passes the identity and coverage thresholds is counted, secondary ones included. A flank whose sequence occurs twice in the short-read assembly, as both haplotypes of a heterozygous strain can, has a low mapping quality and is left unplaced by the registered rule, so its element counts as broken although nothing is broken there.

- `left_flank`, `right_flank`: `placed`, `no full-length match`, `two or more full-length matches`, or `one full-length match, not placed` (mapping quality below 20, or only a secondary alignment).
- `element_reason`: one reason per element, a second copy taking precedence over a missing flank.
- `short_read_synteny_checkable` (centromeres `intact_called` only): the short-read contig extends 10 kb beyond both ends of the call, the arm 1 definition.

The summary counts elements, centromeres and null windows, by zygosity, kind, status and reason, with the fraction of that zygosity's elements of that kind.

### `confirmatory.tsv`

The registered outcomes, one row per measure.

- `outcome`: C1 to C5; `measure`: what is estimated.
- `n_strains`, `n_centromeres`, `k`: strains, centromeres and the count in the numerator; `n_null`, `k_null`: null windows and broken null windows (C4).
- `estimate`, `ci_low`, `ci_high`, `interval`: the estimate and its 95 per cent interval, which resamples strains (10,000 resamples, seed 20261009), except for C1 (Wilson or Clopper-Pearson over strains).
- `p_value`, `test`, `permutations_used`: the unadjusted p-value of a registered test (C4: stratified permutation of the log break-rate ratio within strains; C5: likelihood ratio). Holm-adjusted values are in `results/reporting_summary.tsv`.
- `n_breakpoints`, `loglik_uniform`, `loglik_arm1_at_weighted`, `loglik_ml`, `decision` (C5): breakpoints used, the log-likelihood of the breakpoints under the uniform model, under arm 1's AT-weighted model and at the maximum-likelihood gamma, and the decision by the registered rule.

## Arm 3, read depth (`results/arm3/`, written by `10_compare_pairs.py`)

The reads, subsamples and assemblies are made by `11_assemble_reads.sh`, which logs each strain in `logs/arm3_strains.tsv`. Each assembly is compared with the strain's long-read assembly by the arm 2 method.

### `assemblies.tsv`

One row per assembly attempted: strain, `assembler` (`spades` or `megahit`), `target_depth` (5, 10, 20, 40 or `full`, the whole run capped at 80x), `seed`, the subsampling `fraction`, read pairs, bases before and after fastp and the realised depths (bases over the long-read assembly's length), the check that mate names agree (`names_checked`, `name_mismatches`), `status` and `error`, the assembler's run time and peak memory, contig N50, L50, count and total length, the long-read call count (`long_calls`) and PCAn's call count on the assembly (`calls`).

### `centromere_status.tsv`

As in arm 2, one row per long-read centromere per assembly, for the primary analysis, with `assembler`, `target_depth`, `seed` and `variant` added.

### `confirmatory.tsv`

The registered outcomes D1, D2, D4 and D5 per assembler, target depth and seed, with the columns of the arm 2 table; D1 and D2 rows add the median and interquartile range over strains. D3, the settling depth, is computed from D1 and D2 in `results/reporting_summary.tsv`.

## Sensitivity (`results/sensitivity/`)

### `arm0_gate_threshold.tsv`

The gate decision at the pre-registered threshold (0.95) and at 0.90 and 0.99, and whether each alternative gives the same decision.

### `arm1_at_weight_gamma.tsv` (written by `08_score_perturbations.py`)

Arm 1's AT-weighted fragmentation repeated with other values of its weight gamma, with the columns of `results/arm1/recall_pooled.tsv` and the gamma used.

### `arm2_variants.tsv`

Outcomes C1 to C4 under each registered alternative (`variant`): flanks of 1 kb or 5 kb, null AT bins half or twice as wide, placement identity 0.90 or 0.99, and the heterozygous strains excluded. Columns as in `results/arm2/confirmatory.tsv`.

### `arm3_variants.tsv`

The same alternatives for arm 3 at seed 11: D1 and D4 per variant, assembler and target depth.

### `summary.tsv` (written by `12_analyse.py`)

Every alternative beside the primary result, and whether it changes the conclusion by the rule set before any sensitivity result was read: for a registered test, a different Holm-adjusted significance at 0.05 or a different direction of effect; for an estimate without a test, an alternative estimate outside the primary 95 per cent interval.

## Across arms (written by `12_analyse.py`)

### `reporting_summary.tsv`, `reporting_summary.md`

Every claim the top-level README makes from a registered or pre-specified analysis: arm, status (confirmatory or pre-specified), claim, n, unit of replication, test, estimate, 95 per cent interval, unadjusted (`p_raw`) and Holm-adjusted (`p_holm`) p-values, and the source table. The Holm family is the fourteen tests registered in the analysis plan; `holm_tests_available` says how many of them had been run when the table was written.

### `timing_by_stage.tsv`

A copy of `logs/resources.tsv`: per stage run, its start and end (UTC), elapsed seconds, peak memory, the data directory's footprint at the start and at its peak, the lowest free space on the drive, threads, jobs and exit status.
