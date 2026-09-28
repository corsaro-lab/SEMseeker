# semseeker NEWS

## semseeker 0.99.5 (development)

### Breaking changes

- **`association_analysis()` runs one aggregation per request, and it has to be
  named.** `SCOPE` has two values: `SAMPLE` reduces the positions of a region
  class to one number per sample, `INSTANCE` reduces them to one number per
  instance of that class: one row per gene, per island, per probe. Until now
  every request produced **both** and wrote them into the same file, so each
  result was the union of two different questions, with no way to have asked for
  only one.

  `inference_details` therefore gains a required `scope` column, `"SAMPLE"` or
  `"INSTANCE"`. It has no default: any default would answer half of what a
  previous script asked and say nothing about the other half.

  ```r
  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    aggregation          = "SUM",
    scope                = "SAMPLE"      # or "INSTANCE"
  )
  ```

  **Migration.** A request that wants both branches writes two rows, one per
  scope; both still land in the same result file, told apart by the `SCOPE`
  column. `limma_`/`voom_` families are refused at `scope = "SAMPLE"`: they
  estimate a prior variance across instances, and a collapsed artefact holds one
  row, whereas before the collapsed keys were skipped and the run continued on
  the per-instance ones.

- **`inference_details$scopes` is removed.** It named the region classes a
  second time. They are the `(AREA, SUBAREA)` pairs of the run, declared with
  the `areas` and `subareas` arguments of `association_analysis()` and built at
  runtime; both scopes range over the same pairs. A request still carrying
  `scopes` stops with a message naming the column and what replaced it.

- **A p-value adjustment now names the family it was controlled over, and there
  are three.** `PVALUE_ADJ` said nothing about its family, and at
  `SCOPE = SAMPLE` that family holds a single row — so the column equalled the
  raw p-value by arithmetic while its name promised a correction. It is replaced
  by three columns, nested, each named after its own family:

  ```
  PVALUE_ADJ_KEY_<m>    the identity key: MARKER, FIGURE, SCOPE, AREA,
                        SUBAREA, AGGREGATION — members are the instances
  PVALUE_ADJ_SCOPE_<m>  every row of the same SCOPE — members are its region
                        classes, figures and aggregations
  PVALUE_ADJ_ALL_<m>    the whole file — members are everything tested on this
                        marker under this model (unchanged)
  ```

  The middle level is the one that was missing: a collapsed artefact had only a
  correction over one row, or one shared with the tens of thousands of
  per-instance rows in the same file, and neither answers how many things were
  tested on that sample.

  `<m>` is the method the run declared. All three now honour
  `multiple_test_adj`; the narrow level had `BH` written into it, so a run
  asking for another estimator got a column computed as BH and named for it.
  Where an estimator cannot work on a family — `qvalue` needs enough p-values to
  estimate pi0 — the answer is `NA` and a line in the log, rather than a
  silently substituted estimator that would make the column name wrong.

  **Migration.** Read `PVALUE_ADJ_KEY_<m>` wherever you read `PVALUE_ADJ`. Files
  written by an earlier release keep their `PVALUE_ADJ` column until the folder
  is recomputed, and it is dropped when they are. The three answer three
  different questions and are **not** a severity ladder: Benjamini-Hochberg is
  adaptive, so a wider family is usually but not necessarily more conservative.


- **The extent a number is valid over is a coordinate of its own, `SCOPE`.**
  A burden over the whole sample and a burden per gene are the same
  marker reduced over different extents. They used to live in artefacts of
  different *shape* — a sibling CSV with samples down the rows, and a pivot with
  areas down the rows — and that difference in shape is what hid the difference
  in meaning. Both are now pivots, and the name says which is which:

  ```
  <MARKER>_<FIGURE>_<SCOPE>_<AREA>_<SUBAREA>_<AGGREGATION>_<GENOME_BUILD>

  MUTATIONS_HYPER_SAMPLE_PROBE_WHOLE_SUM_HG19     burden of the whole sample
  MUTATIONS_HYPER_SAMPLE_GENE_WHOLE_SUM_HG19      same, over gene probes only
  DELTAS_HYPO_INSTANCE_GENE_TSS200_MEDIAN_HG19    median per gene, TSS200 window
  SIGNAL_BETA_INSTANCE_PROBE_WHOLE_VALUE_HG19     the beta value per probe
  ```

  `SCOPE = SAMPLE` collapses the region class to one number per sample — such an
  artefact is one row tall — while `SCOPE = INSTANCE` keeps one row per gene,
  island or probe.

- **The aggregation requested at `SCOPE = INSTANCE` is now the one you get.** It
  was validated and then dropped on the way to the read, so asking for the
  median of a region returned its mean: the file existed, its name said nothing
  about which operator had produced it, and nothing complained. Result folders
  from earlier versions do not match the new names and recompute once.

- **`VALUE` names the identity.** A `PROBE` or `POSITION` row already holds a
  single position, so there is nothing to reduce; `VALUE` says so, where before
  the aggregation segment was simply absent. Its absence is again an error
  everywhere else.

- **`SAMPLE_STATS_RESULT.csv` is gone.** Its content is the `SCOPE = SAMPLE`
  artefacts, and the readable per-sample table is composed on read by
  `sem_study_summary_get(regions = ...)`, which joins them onto the sample
  sheet. `sem_study_summary_get()` with no arguments behaves as before.

- **`semseeker(sample_stats_scopes = ...)` was removed.** Which region classes
  you want is no longer a decision to be made before the run: the artefacts are
  built when they are asked for. The error that used to say *"produce it with
  `semseeker(sample_stats_scopes = ...)` and rerun the analysis"* is gone —
  changing your mind now costs one scan of the position pivot instead of a whole
  SEM run. Ask for a region class at analysis time with
  `association_analysis(inference_details$scopes)` or
  `sem_study_summary_get(regions = ...)`.

  The cost moved rather than vanished: the call that asks for a region class has
  to know it, so `association_analysis(inference_details$scopes = "GENE_TSS1500")`
  needs `areas` and `subareas` to cover `GENE`/`TSS1500` — a run that does not
  register the pair refuses the request instead of silently testing nothing. The
  registry is also what resolves `GENE_TSS1500` into its two coordinates without
  splitting the string, which cannot be done safely: `N_SHORE` carries an
  underscore of its own.

- **`MODE_LOW` / `MODE_HIGH` are spelled `MODELOW` / `MODEHIGH`**, and they are
  admissible only at `SCOPE = SAMPLE`. Per instance a gene holds about nineteen
  probes and a TSS200 window two to five: a two-peak density estimate on those
  is shaped by the bandwidth rather than by the data, and the old code returned
  that number instead of refusing. A minimum-numerosity guard now returns `NA`
  rather than a plausible-looking value.

- **`depth_analysis` is gone, and so is the `DEPTH` column.** Not derived, not
  deprecated: removed. Once every artefact has the same shape — a key column and
  one column per sample — a model handed a row does not know, and has no reason
  to ask, whether the key of that row is a gene symbol or `PROBE_WHOLE`. It
  fits. The two consumers merged into one (`assoc_run_marker()`, which replaces
  `sem_run_depth1_marker()` and `sem_run_depth_n_marker()`), and the granularity
  is said by `SCOPE`, `AREA` and `SUBAREA` — which say more, because those pairs
  are only *partially* ordered and an integer scale projected a lattice onto a
  line.

- **The `TOTAL` column is gone.** It was `sum` of the per-area aggregates, and
  `depth_analysis == 2` meant "test only that" — but the partition into areas is
  not disjoint, so a probe the annotation maps onto three genes entered that
  total three times. It was the composition of aggregates this release forbids,
  in the one place nobody looked for it. What it meant to compute — the whole
  region class, one number per sample — is now an artefact of its own:
  `SCOPE = SAMPLE` on that `(AREA, SUBAREA)`, derived from the masked position
  pivot where every position counts once.

- **`AREA_OF_TEST` names the aggregate for collapsed artefacts.** At
  `SCOPE = SAMPLE` a row is the whole region class reduced to one number, and
  which class that is is already said by `AREA` and `SUBAREA`; what the row
  needs to say is *which aggregate*. So it carries `MEAN`, `MEDIAN`, `SUM`, …
  Previously it repeated the region class, leaving the mean and the median of
  the same class with the same `AREA_OF_TEST`.

- **The taxonomy key leads the result file**, in the order `MARKER`, `FIGURE`,
  `SCOPE`, `AREA`, `SUBAREA`, `AGGREGATION`, `AREA_OF_TEST`, and a missing value
  in any of them stops the write. The key is the identity of a row: two rows with
  a hole in the same place cannot be told apart. The old code filled such holes
  instead — `SUBAREA` became the literal `"TOTAL"`, a value outside its own
  vocabulary — which is how a defect hides inside a key.

- **`N_PROBES` moved to `SAMPLE_SHEET_RESULT.csv`.** It is a property of the
  imputation — how many positions of that sample survived the treatment of
  missing values — not an aggregation of a marker, so it sits with the
  descriptive properties of the sample. Nothing is lost for the density: the
  `MEAN` of a binary marker *is* the density, denominator included.

- **Every computed quantity is now named by four coordinates.** A
  number produced by SEMseeker is the reduction of a set of genomic positions,
  and until now the operator that reduced them was implicit — one per marker, so
  the marker and its direction identified the value. A scope can carry several
  reductions, so the operator became an axis of its own and the names changed
  accordingly:

  ```
  <SCOPE>_<MARKER>_<FIGURE>_<AGGREGATION>
  ```

  | before | after |
  |---|---|
  | `SAMPLE_MUTATIONS_HYPER` | `SAMPLE_MUTATIONS_HYPER_SUM` |
  | `SAMPLE_DELTAS_HYPO` | `SAMPLE_DELTAS_HYPO_MEAN` |
  | `SAMPLE_MEDIAN` | `SAMPLE_SIGNAL_BETA_MEDIAN` |
  | `SAMPLE_MODE_LOW` | `SAMPLE_SIGNAL_BETA_MODELOW` |
  | `SAMPLE_N_PROBES` | moved to `SAMPLE_SHEET_RESULT.csv` as `N_PROBES` — see the entry above |

- **The figure of `SIGNAL` is the scale of the values, `BETA` or `MVALUE`.** It
  used to be `MEAN`, a placeholder that made the key unique while describing a
  way of aggregating — which now has an axis of its own. Pivot files carry it,
  so a beta run and an M-value run no longer write the same file name and
  overwrite each other in one folder. WGBS, ONT and PacBio stay `BETA`: a
  methylation fraction from reads is the same bounded scale as an array beta
  value, and the package compares the two on purpose.

- **Pivot file names carry the aggregation**, e.g.
  `MUTATIONS_HYPER_CHR_CYTOBAND_SUM_HG19.parquet`. Result folders produced by
  earlier versions do not match the new names and recompute once. This buys a
  guarantee that did not exist: the name used to say nothing about which
  operator had produced the file, so an existing pivot was reused on trust.

- **`inference_details$aggregation` is required.** A request that does not name
  the reduction does not identify what it wants tested. A malformed request (no
  aggregation, or a name outside the taxonomy) stops the run; a request that no
  marker of the run admits — the median of a 0/1 marker, say — drops that row
  with a warning naming it and lets the other rows proceed.

- **Results carry an `AGGREGATION` column**, part of the row identity along with
  marker, figure, area and subarea, in the deduplication, in the resume match
  and in the cross-study overlaps.

### Bug fixes

- **A request for `multicore` on a platform without `fork()` was granted with a
  single worker.** `future::plan(multicore)` does not refuse such a platform: it
  accepts the request and degrades to one worker. Windows has no `fork()`, so a
  run taking the default `parallel_strategy = "multicore"` executed the whole
  pipeline in series while believing it was parallel, and nothing in the log
  said so. macOS already had an explicit conversion, for an unrelated reason
  (fork is unsafe there alongside Polars); a platform that simply lacks `fork()`
  had none. `core_parallel_session()` now converts `multicore` to `multisession`
  wherever `parallelly::supportsMulticore()` is false, and logs the substitution
  the way the macOS path does. Measured on CI, same commit and same suite: 285
  minutes on Windows against 134 on Linux. The ratio was the number of workers,
  not the speed of the platform.

- **A per-sample burden restricted to a region class was computed over the whole
  sample.** At `SCOPE = SAMPLE` the mask that selects the positions of a class
  was built from the full probe annotation table, including the probes with no
  annotation for that class. On a 450k run the mask covered all 485,512
  positions whatever the class was, against 365,860 for `GENE_WHOLE`, 84,342
  for `GENE_TSS1500` and 62,870 for `ISLAND_N_SHORE`, so every restricted
  burden equalled the burden of the whole sample, and all classes returned the
  same number. Per-instance results were never affected, and neither were
  WGBS/long-read runs, where the annotation is resolved by coordinate overlap
  and carries only the positions of the class.

  **Migration.** Recompute the collapsed artefacts and any inference run on
  them: delete the `SCOPE = SAMPLE` pivots of the result folder, or rerun the
  analysis with `start_fresh = TRUE`. Per-instance results and the SEM layer
  need no action.

- **A consumer wrote over the file it was consuming.**
  `assoc_data_extractor()` read the canonical inference CSV, applied three
  filters — the figures of the marker, the `areas_sql_condition` of the request,
  the dropped `SAMPLES_SQL_CONDITION` columns — and wrote the filtered frame
  back over the same path. Extracting was therefore destructive: a narrow
  `areas_sql_condition` permanently reduced the result of the run, and the next
  iteration of its own loop read a file a previous iteration had truncated. It
  writes only to `destination_folder` now. `sem_metrics_name_collect()` left the
  read path with it — its body is commented out in full, but what it used to do
  is why it does not belong there: it wrote a registry to disk while a consumer
  was reading.

- **`assoc_results_get()` guessed which region class and which extent to read.**
  `area` defaulted to `"GENE"` and `scope` to `NULL`, which meant no
  filter at all. Both are required now, and removing the defaults immediately
  showed what they had been hiding: two enrichment backends declared neither, so
  they were reading every extent — the collapsed `SCOPE = SAMPLE` rows
  included — into a gene set, and four cross-study overlap call sites mixed the
  collapsed row of a region class with its per-instance rows. All seven now
  declare what they read.

- **Two dead adjustment flags removed.** `adjust_per_area` and
  `adjust_globally` re-ran `p.adjust()` at read time on top of an already
  adjusted column — `pvalue_column` defaults to `PVALUE_ADJ_ALL_BH`, so either
  one meant BH over BH. No call site passed `TRUE`. `adjust_per_area` also
  reused the name `area` as its loop counter, shadowing the parameter, so the
  filter that followed selected the last area of the loop rather than the
  requested one.

- **The significance flags matched the method string, not the level.**
  `SIGNIFICATIVE_ADJ_ALL` and `SIGNIFICATIVE_ADJ` were computed over *every*
  column whose name contained the method — `grepl("BH", colnames)` — so a second
  adjusted level would have turned them into "significant at every level at
  once" without a visible change. They now select the `ALL` level by name. An
  empty selection yields `NA` instead of `TRUE`, which is what
  `all(logical(0))` used to return.

- **`filter_p_value` filtered on a column that is never created.**
  `assoc_analysis_save_results()` subset on `SIGNIFICATIVE_ADJ`, which is built
  by `assoc_results_get()` on the way out and never exists on the way in. The
  parameter defaults to `TRUE`, so the default path could not complete; every
  fixture in the suite sets it to `FALSE`, which is why nothing caught it. It
  filters on `SIGNIFICATIVE_ADJ_ALL`, the flag that exists at that point.

- **Two models adjusted p-values over a batch.**
  `assoc_apply_stat_model()` and `assoc_glm_model_bulk()` each wrote a
  `PVALUE_ADJ` over whatever chunk they had been handed, making the family a
  memory parameter; the value was then recomputed downstream anyway. Both
  removed. The one in `assoc_apply_stat_model()` also split its rows on
  `grepl("TOTAL", AREA_OF_TEST)`, and `TOTAL` was removed in this release — the
  predicate had been false on every row since, so its two branches had quietly
  become one.

### New features

- **Descriptors on any scope.** The signal descriptors are no longer a separate
  path: `SIGNAL` is a marker like the others, so `sample_stats_scopes` now
  produces its median, mean, variance, IQR and — on the beta scale — its two
  modes restricted to a region class, e.g. `GENE_TSS1500_SIGNAL_BETA_MEDIAN`.

- **Density.** The mean of a binary marker over a scope is the fraction of its
  positions classified as epimutations. Unlike the raw burden it is comparable
  between region classes of different size, which the burden is not.

- Which reductions a marker admits is declared in one place. Continuous markers
  (`SIGNAL`, `DELTAS`, `DELTAR`) take the full descriptor set; counts take the
  sum and the mean, because their median and IQR are degenerate on a vector of
  zeros and ones. The two modes assume a bounded bimodal scale and are therefore
  restricted to the signal on the beta scale.


- **Per-sample burden restricted to a region class.** The
  statistics sibling can now carry a scope other than the whole sample:
  `semseeker(sample_stats_scopes = c("SAMPLE", "GENE_TSS1500"))` adds
  `GENE_TSS1500_<MARKER>_<FIGURE>`, the burden computed over the probes of that
  region class only. Any registered `(AREA, SUBAREA)` pair of the run is a
  legal scope; `SAMPLE` is always produced. A scope that does not resolve stops
  the run instead of being ignored.

  Each position is counted **once**, even when the annotation maps it to
  several genes. The burden is computed from the POSITION pivot restricted to
  the probes of the region, not by summing the per-gene pivot — summing the
  latter would count a multi-gene probe once per gene.

- **Region scopes as dependent variable at `depth = 1`.** The new
  `inference_details$scopes` column (several separated by `"+"`, default
  `"SAMPLE"`) selects which scopes `association_analysis()` tests at depth 1.
  A region scope is still one value per sample, so its rows keep `DEPTH = 1`
  and `SUBAREA = "SAMPLE"`, with the scope in `AREA` (e.g. `GENE_TSS1500`) and
  `AREA_OF_TEST` unchanged. Requesting a scope that was never produced is an
  error naming the scope, not a silently empty result.

### Documentation

- The vignettes now show the names this version writes: `SIGNAL` pivots carry
  the scale of the run (`SIGNAL_BETA_*`, or `SIGNAL_MVALUE_*`), file names are
  uppercase, and the optional `AGGREGATION` segment is part of the documented
  pattern. `sem_coverage_analysis_report()` describes its `signal_data` path
  the same way.

## semseeker 0.99.4

### Breaking changes

- **The per-sample burden moved out of `SAMPLE_SHEET_RESULT.csv` into a new
  sibling file, `SAMPLE_STATS_RESULT.csv`.** The columns that used to
  be appended to the sample sheet (`MUTATIONS_HYPER`, `DELTAS_HYPO`, …, and
  `PROBES_COUNT`) were aggregated over every probe with no genomic filter —
  that is the *sample* scope of the new artefact — and are now written there as
  `SAMPLE_<MARKER>_<FIGURE>` and `SAMPLE_N_PROBES`. The sample sheet stays the
  description of the study; the two files join on `Sample_ID`, and
  `sem_study_summary_get()` performs that join for you, so analyses inside the
  package are unaffected. Code that read the burden straight from
  `SAMPLE_SHEET_RESULT.csv` must read the sibling instead. *(Superseded in
  0.99.5: the sibling file no longer exists, the join is composed from the
  `SCOPE = SAMPLE` artefacts, and `N_PROBES` returned to the sample sheet.)*

### New features

- **Per-sample signal descriptors.** Alongside the burden, the new
  sibling carries `SAMPLE_MEDIAN`, `SAMPLE_MEAN`, `SAMPLE_VARIANCE`,
  `SAMPLE_IQR` and `SAMPLE_N_PROBES` for every sample. On the beta scale it
  also carries the two modes of the bimodal distribution, `SAMPLE_MODE_LOW` and
  `SAMPLE_MODE_HIGH`, estimated as the highest density peak on each side of
  0.5. On M-values the distribution is unimodal and the split carries no
  meaning, so the two columns are omitted rather than filled with a
  meaningless number. Column names are composed by a single shared helper, so
  producer and consumer cannot drift apart. Region scopes
  (`<AREA>[_<SUBAREA>]_*`) are the next step; the naming already accommodates
  them.

- **Coverage is now a mandatory pre-step of every SEM analysis.**
  Coverage used to be an opt-in report, so a run could go straight to SEM
  detection on input that barely overlaps the reference annotation (long-read
  positions against an Illumina manifest, a batch with a probe mismatch, the
  wrong genome build) and produce thresholds computed on a handful of probes.
  Every run now: writes the coverage charts, logs a coverage banner
  (`input_positions` / `covered_by_reference` / percentage), writes a
  `COVERAGE_GATE.json` sidecar for audit, and **stops** when coverage falls
  below the new `coverage_minimum` parameter (default 80%). The charts are
  produced even when the gate rejects the run — they are what shows which
  regions were missing. Coordinate-based technologies (WGBS, long reads)
  derive their annotation from the positions themselves, so the threshold is
  not applied there. Lower `coverage_minimum` explicitly for deliberate
  cross-technology runs.
- The coverage charts are now computed over the full annotation set instead of
  the areas selected for the run. `areas` defaults to `POSITION` only, which
  the coverage report skips, so on a default run the charts were previously
  never produced.

### Documentation

- `inst/CITATION` now carries the Zenodo DOI `10.5281/zenodo.5095416` instead
  of a Bioconductor DOI that does not resolve (the package is not on
  Bioconductor yet), and the README no longer says a Zenodo DOI "will be
  registered" — the archive has existed since 2021.

### Bug fixes

- **A sample sheet that shares no identifier with the computed pivots no
  longer produces a silently empty result.**
  `sem_study_summary_total()` merged the per-sample burden onto the sample
  sheet by name with `all.x = TRUE`: when the two sides had no identifier in
  common, every burden column landed as `NA` while `PROBES_COUNT` stayed
  populated, and each `depth = 1` inference downstream then failed with
  "data are not the same size". The merge now stops with a message naming
  both sides, and reports partially missing samples as a warning. The
  accumulator inside the function is also a proper local binding now: it was
  resolved with `exists()`, which reaches the global environment.
- The same `exists()` pattern was removed from `sem_coverage_analysis()`,
  which could also raise "object 'tot_result' not found" when no area
  produced a usable table.

- **Sample identifiers containing `-`, `.` or spaces are now handled correctly.**
  Sample-column names of the signal matrix are normalised
  (uppercase, non-alphanumeric characters replaced by `_`), but the sample
  sheet identifiers were used as-is when selecting those columns. With
  identifiers such as `C3L-00001-06` the two sides never matched:
  `util_exploratory_analysis()` wrote a cleaned signal parquet holding the
  `PROBE` column only, and the SEM pipeline stopped with
  `undefined columns selected`. Datasets whose identifiers are already
  alphanumeric (e.g. `GSM123456`) were unaffected. Both sides are now
  normalised through a single, idempotent internal helper before any
  name-based column selection; distinct identifiers that would collapse onto
  the same normalised name, and identifiers without a matching signal column,
  are reported explicitly instead of failing later.

## semseeker 0.99.3

### Documentation

- Aligned the delta-metric documentation in the vignettes and README with the
  implementation: the `DELTAS`/`DELTAR`/`DELTAP`/`DELTARP`/`DELTAQ`/`DELTARQ`
  descriptions now reflect the relative-delta ratio and its equal-width and
  quantile ranked variants. Corrected the pivot file-name examples to
  `Data/Pivots/<MARKER>/<MARKER>_<FIGURE>_<AREA>_<SUBAREA>_<build>.parquet` and
  updated a stale `enrichment_analysis()` reference.

### Bug fixes

- **`plot_box_plot()`: fixed R CMD check ERROR under ggplot2 >= 4.0.**
  `ggpubr::stat_compare_means(label = "p.format")` builds an internal
  `after_stat(create_p_label(...))` expression that ggplot2 >= 4.0 cannot
  resolve when ggpubr is loaded via `::` but not attached, causing
  "could not find function 'create_p_label'" (test-2-box-plot.R). The
  comparison p-value is now computed directly from `stats`
  (`t.test`/`wilcox.test`/`kruskal.test`/`anova`) and drawn with
  `ggplot2::annotate()`, removing the ggpubr rendering-path dependency.
  Behaviour change: the `kruskal.test` branch now shows the overall test
  p-value only; the per-pair `ggpubr` brackets have been dropped. Also fixed
  the `unit=`/`units=` partial-argument-match warning in `ggsave()`.

### Dependencies

- Bumped pinned GitHub Actions: `actions/cache` 4→6, `actions/upload-artifact`
  4→7, `actions/dependency-review-action` 4→5, `codecov/codecov-action` 4→7.

## semseeker 0.99.2

### Breaking changes

- **LESIONS detection now uses genomic distance, not probe count.**
  The legacy `sliding_window_size` parameter (probe-count based, default 11)
  has been REMOVED. It is replaced by `LESIONS_BP` (default 2000), the
  maximum bp distance for two probes to be considered part of the same
  enrichment window. Same registration pattern as `DELTAQ_Q` / `DELTAP_B`:
  registered in `core_init_env()`, surfaces through the `Q` column of
  `keys_markers_default_discrete`, exposed via `semseeker(LESIONS_BP = ...)`.
  Rationale: probe-density varies dramatically across the genome and across
  array platforms (450K ~485k probes, EPIC ~865k, LONGREAD bedmethyl
  variable). A fixed 11-probe window covers ~70kb on average on 450K but
  ~500bp inside a CpG island — biologically incomparable. `LESIONS_BP=2000`
  matches the typical span of a CpG island or DMR and yields a single
  semantics across platforms. Migration: any caller passing
  `sliding_window_size=N` must replace it with `LESIONS_BP=M` (no equivalence
  formula; the metric has changed). Multi-window sensitivity will be tackled
  in a later release (vector-valued `LESIONS_BP`).

### New features

- **CpG-island subareas aligned to Illumina `Relation_to_Island`.**
  The `ISLAND` area now exposes all six Illumina contexts plus the whole
  neighbourhood: `WHOLE`, `ISLAND`, `N_SHORE`, `S_SHORE`, `N_SHELF`,
  `S_SHELF`, `OPENSEA`. Previously the island core (`Island`) and the
  open-sea compartment (`OpenSea`) were lost. `ISLAND_WHOLE` is **redefined**
  to mean the whole island neighbourhood (core + shores + shelves, ±4 kb),
  mirroring `GENE_WHOLE`; use the new `ISLAND` subarea for the core alone.
  `OPENSEA` groups each open-sea CpG by the inter-neighbourhood genomic gap
  that contains it, labelled `chr:start-end` (never spanning a chromosome).
  Semantics are centralised in `island_opensea.R` and shared by both the
  Illumina (`anno_probe_annotation_build`) and coordinate/AnnotationHub
  (`anno_area_granges_build`) backends. See the getting-started vignette.

- **`binomial_bulk` family + goodness-of-fit metrics extension.**
  New `family_test = "binomial_bulk"` dispatches to `assoc_glm_model_bulk()` for
  bulk per-probe logistic regression via `Rfast::glm_logistic` (parallelised
  with `foreach %dorng%`), ~10-20× faster than the per-probe `stats::glm`
  path. Drop-in replacement for `family_test = "binomial"` at PROBE-level
  inference details (LESIONS, MUTATIONS). Same legacy schema: per-coef
  PVALUE/ESTIMATE columns plus top-level PVALUE/PVALUE_ADJ.

  Standard errors are derived from the Fisher information at the MLE
  (`Var(β̂) = (X' diag(p(1-p)) X)^{-1}`), matching `stats::glm` Wald output.

  `io_data_preparation()` gains a universal degenerate-burden filter: columns
  with `var(Y) == 0` are dropped before reaching the model. Critical for
  LESIONS @ PROBE where ~92% of probes are all-zero across samples
  (manifest-aligned pivot, retained for positional join with annotations).
  The filter applies to every family — binomial GLM no longer diverges on
  constant Y, limma/voom no longer produces NaN t-stats, polynomial no
  longer fits rank-deficient designs.

  Ten new metrics registered in `metrics_properties.rda`:

  | Metric | Engine | Direction | Notes |
  |---|---|---|---|
  | `T_STAT_MODERATED` | limma_2, voom_2 | Higher better | eBayes-moderated t-stat per coef |
  | `B_STATISTIC` | limma_2, voom_2 | Higher better | log-odds of differential expression (lods) |
  | `F_STAT_MODERATED` | limma_2, voom_2 | Higher better | joint F across non-intercept coefs (parabolic test) |
  | `POSTERIOR_RESIDUAL_VAR` | limma_2, voom_2 | Lower better | s2.post diagnostic |
  | `MCFADDEN_R2` | binomial_bulk | Higher better | 1 − devi/null_devi pseudo-R² |
  | `NAGELKERKE_R2` | binomial_bulk | Higher better | scaled Cox-Snell pseudo-R² |
  | `C_STATISTIC_AUC` | binomial_bulk | Higher better | discrimination (= AUC) |
  | `DEVIANCE_RATIO` | binomial_bulk | Lower better | devi / null_devi |
  | `AIC_VALUE` | all GLM | Lower better | uppercase canonical of legacy `aic_value` |
  | `BIC_VALUE` | all GLM | Lower better | new |

  The metrics replace R²/R²_adj for limma/voom (lmFit doesn't return R²
  natively, and the eBayes-moderated diagnostics are more informative)
  and provide the analogous goodness-of-fit signal for binomial logistic
  regression (R² doesn't apply to {0,1} outcomes).

- **Session provenance metadata** (C-06).
  Every `semseeker()` run now writes `session_metadata.json` to the result
  folder at the start of analysis:
  `{"genome_build":"hg19","tech":"K850","semseeker_version":"0.99.0","created":"...","sample_n":120}`
  The file is human-readable and parseable without R, suitable for pipeline
  auditing and automated compatibility checks.
  Pivot files (parquet) additionally receive a sidecar `*_meta.json` with the
  same build/tech stamp; pivot file names now include the genome build as a
  suffix before the extension (e.g. `MUTATIONS_HYPER_GENE_TSS1500_hg19.parquet`)
  as belt-and-suspenders provenance.
  Inference CSV outputs gain two constant columns — `GENOME_BUILD` and `TECH`
  — that survive any downstream merge or stack operation.
  New exported function `core_check_session_compatibility(session_list)`:
  **stops** if `genome_build` differs across sessions (physically incomparable
  coordinates); **warns** if `tech` differs (cross-array meta-analysis is valid
  on the probe intersection but must be explicit).
  `assoc_intra_study_association_replication()` now calls this guard internally: stops
  with a clear message when origin results carry a different `GENOME_BUILD` than
  the current session.

### Breaking changes

- **`semseeker()` no longer auto-converts M-values to beta.**
  The `auto_convert_mvalues` parameter has been removed and the
  `.looks_like_mvalues()` helper deleted. Input values are now passed
  through unchanged; if you need conversion, call `sem_mvalue_to_beta()`
  explicitly before `semseeker()`. The beta-vs-M-value flag is still
  detected and stored in `ssEnv$beta` by `core_get_meth_tech()`, so downstream
  code that consults that flag continues to work.

- **`start_fresh` defaults to `FALSE`.**
  `core_init_env()` no longer deletes the result folder before starting. Existing
  results are preserved unless the caller explicitly passes `start_fresh = TRUE`.
  The Shiny UI exposes this as a checkbox ("Delete result folder before running",
  unchecked by default).

### Bug fixes

- **macOS: tests default to `multisession` instead of `multicore`** (E-14).
  `multicore` (fork) is unsafe on macOS with Polars' C++ thread pool — forked
  children can be killed by Mach exceptions. `setup.R` now selects `multisession`
  on Darwin, `multicore` on Linux. All `%dorng%` foreach bodies now call
  `core_update_session_info(ssEnv)` as their first statement to populate the worker's
  `.pkgglobalenv` — required because `multisession` workers are fresh R processes
  where the session singleton is empty.

- **`chr` prefix mismatch in Polars join** (E-13).
  `io_dump_sample_as_bed_file()` prepends `chr` to chromosome names, but
  `signal_thresholds` retains bare numbers from probe annotations. The inner
  join in `util_join_values_to_thresholds()` returned 0 rows → no mutations detected.
  Fix: `util_join_values_to_thresholds()` now strips the `chr` prefix from both sides
  before joining.

- **`exists("signal_data")` scoped to local environment in `sem_analyze_batch()` and
  `sem_analyze_population()`** (E-01).
  The previous `exists("signal_data")` used `inherits = TRUE` (R default), which walked
  all the way up to `.GlobalEnv`. If the user had a `signal_data` object in their session
  (the normal case — they load data before calling `semseeker()`), the guard fired
  even though the local copy had already been freed, and the subsequent `rm(signal_data)`
  produced a spurious `"object not found"` warning. In the worst case (interactive session
  where local and global scopes coincide) it could silently delete the user's data.
  Fix: `exists("signal_data", envir = environment(), inherits = FALSE)` +
  `rm("signal_data", envir = environment())` in both files.

- **Polars inner join replaces positional zip in `sem_mutations_get()`, `sem_delta_single_sample()`,
  `sem_deltar_single_sample()`; coverage banner in `sem_analyze_population()`** (A-10).
  The previous implementation sorted both `values` and `thresholds` by CHR/START/END and
  then compared them row-by-row (positional zip). When `signal_thresholds` came from a
  different run (cross-run analysis, e.g. Nanopore sample vs Illumina reference batch passed
  via `populationControlRangeBetaValues`), the two position sets could differ in size or
  overlap, producing silently wrong mutation/delta calls or an out-of-bounds crash.
  Changes:
  - Added `util_join_values_to_thresholds()` — a private helper that performs a Polars lazy inner
    join on `(CHR, START, END)`. Shared by all three per-sample functions to ensure consistent
    intersection logic. Required for Nanopore bedmethyl files (28M+ rows where base-R
    `merge()`/`match()` is too slow).
  - `sem_mutations_get()`, `sem_delta_single_sample()`, `sem_deltar_single_sample()` now use the helper;
    zero-overlap guard returns an empty result rather than crashing.
  - Coverage banner moved from per-sample (inside `sem_mutations_get()`) to per-batch (at the
    start of `sem_analyze_population()` before the per-sample loop). Emitted once per batch:
    `input_positions | beta_range_positions | covered_by_inner_join`.

### Bug fixes (A-09: assoc_bayes_analysis rewrite)

- **`assoc_bayes_analysis()`: 9 bugs fixed** (A-09).
  - **Loop off-by-one** (bug 1): `for (a in length(markers))` iterated only once
    (`a = 2`, only "LESIONS"), silently skipping "MUTATIONS". Fixed with `seq_along()`.
  - **Wrong file path** (bug 2): `read_delim(io_pivot_file_name)` passed the function
    object instead of the local variable `pivot_filename` → runtime connection error.
  - **Missing assignment** (bug 3): `tempDataFrame` used before being assigned from
    `pivot` → "object not found" crash. Added `tempDataFrame <- pivot`.
  - **`subset()` with string literal** (bug 4): `subset(df, "Sample_Group" != "Reference")`
    is always `TRUE` → Reference samples never filtered → contaminated Bayes estimates.
    Fixed to bare column reference `subset(df, Sample_Group != "Reference")`.
  - **Column drop off-by-one** (bug 3 cont.): `[, c(-1,-3)]` dropped the first data
    column or the `independent_variable` column depending on session config. Replaced
    with `colnames != "Sample_ID"` (remove only the merge key).
  - **`exists()` without scope** (bug 5): same pattern as E-01; fixed with
    `inherits = FALSE, envir = environment()`.
  - **Duplicate `max()` computation** (bug 6): `max_P_case` was computed twice,
    `max_P_control` never computed → filtered output missed Control criterion.
  - **Output filename typo** (bug 7): `"bayes_analisys"` → `"assoc_bayes_analysis"`.
  - **`c` as loop variable** (bug 8): foreach variable named `c` shadowed the base
    `c()` function. Renamed to `col_idx`.
  - **Hardcoded thresholds** (bug 9): 0.9 / 0.1 are now `bayes_case_threshold` and
    `bayes_control_threshold` parameters with the original values as defaults.

### Statistical model changes

- **`sem_lesions_get()`: replaced hypergeometric with binomial test** (A-01).
  The sliding window advances one probe at a time, so each probe participates
  in up to `sliding_window_size` consecutive windows — sampling with
  replacement. The hypergeometric distribution assumes sampling without
  replacement and produced inflated (too-small) p-values. The new test is:
  `P(X ≥ ENRICHMENT | Binomial(sliding_window_size, p0))` where
  `p0 = MUTATIONS_COUNT / PROBES_COUNT` is the empirical background
  mutation rate for the grouping unit (gene, chromosome, …). Expected
  impact: more conservative lesion calls, better calibration for samples
  with low or high global methylation variation.

## semseeker 0.11.0

### New features

- Added three pkgdown vignettes: Getting started, Association analysis (all 15+
  model families), Pathway and enrichment analysis (all 6 backends).
- Added GitHub Actions CI matrix: macOS, Ubuntu, Windows (`R-CMD-check.yml`).
- Added test coverage workflow (`test-coverage.yml`) with Codecov upload.
- Added `./ci-local.sh` for local Docker-based CI reproduction (`check` and
  `coverage` modes).

### Bug fixes

- Fixed `future::plan(multicore, workers = 0)` crash when `availableCores()`
  returns 1 (e.g. in covr subprocess): added `max(1L, nCore)` guard in
  `parallel_session.R`.

### Documentation

- Corrected SEM citations: replaced Teschendorff with correct attribution to
  Gentilini et al. 2015 (doi:10.18632/aging.100792) and Corsaro et al. 2023
  (doi:10.3390/cancers15164109).
- Added differential signal analysis section (SIGNAL_MEAN) to
  getting-started vignette.

## semseeker 0.10.0 and earlier

See git log for earlier changes.
