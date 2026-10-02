# Package index

## All functions

- [`anno_annotate_position_pivots()`](https://corsaro-lab.github.io/SEMseeker/reference/anno_annotate_position_pivots.md)
  :

  Materialise the area-level artefacts of the run Rewritten as a thin
  batch pass over \[io_pivot_build()\], which is now the single
  derivation path of the taxonomy. What used to live here - the join
  with the annotation, the multi-gene explode, the group-by - moved into
  \`io_pivot_build()\` and \[.anno_area_explode()\], so the batch
  pre-pass and the on-demand construction cannot diverge. What
  disappeared with it is the rule that chose the operator:

       if (localKeys[i, "DISCRETE"]) pivot$group_by("AREA")$sum() else pivot$group_by("AREA")$mean() 

  One file per key, the operator picked by a flag and \*\*absent from
  the name\*\*. Downstream, \`inference_details\$aggregation\` was
  validated and then ignored: asking for \`MEDIAN\` on \`GENE_TSS1500\`
  returned the mean, silently. The keys now carry the aggregation, so
  what is written says which operator produced it. This pass is a
  convenience, not a requirement: anything it does not materialise is
  built on first read by \[io_read_pivot()\].

- [`anno_area_granges_build()`](https://corsaro-lab.github.io/SEMseeker/reference/anno_area_granges_build.md)
  : Build a GRanges object for a given genomic area/subarea

- [`anno_manhattan_plot_marker_per_probe()`](https://corsaro-lab.github.io/SEMseeker/reference/anno_manhattan_plot_marker_per_probe.md)
  : Title anno_manhattan_plot_marker_per_probe

- [`anno_position_pivot_to_probe()`](https://corsaro-lab.github.io/SEMseeker/reference/anno_position_pivot_to_probe.md)
  : Convert a POSITION-keyed signal pivot into a PROBE-keyed lazy frame.

- [`anno_probe_features_get()`](https://corsaro-lab.github.io/SEMseeker/reference/anno_probe_features_get.md)
  : Retrieve probe feature annotations for a given genomic area

- [`anno_sort_by_chr_and_start()`](https://corsaro-lab.github.io/SEMseeker/reference/anno_sort_by_chr_and_start.md)
  : sort the dataframe using CHR and START sorting column first for CHR
  and after for START

- [`assoc_apply_stat_model()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_apply_stat_model.md)
  : Title

- [`assoc_bayes_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_bayes_analysis.md)
  : Bayesian posterior probability analysis of SEMseeker mutations and
  lesions

- [`assoc_compute_mean_delta_permutation_cpu()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_compute_mean_delta_permutation_cpu.md)
  : Compute mean delta permutation (CPU)

- [`assoc_compute_quantile_delta_permutation()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_compute_quantile_delta_permutation.md)
  : Compute quantile delta permutation (CPU)

- [`assoc_compute_quantreg_permutation()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_compute_quantreg_permutation.md)
  : Title

- [`assoc_compute_spearman_permutation()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_compute_spearman_permutation.md)
  : Compute Spearman permutation statistic (CPU)

- [`assoc_exact_pvalue()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_exact_pvalue.md)
  : Quantile regression result value, confidence interval and p.value

- [`assoc_glm_model()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_glm_model.md)
  : GLM association model

- [`assoc_mean_permutation()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_mean_permutation.md)
  : Mean permutation model

- [`assoc_mediation_linear_model()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_mediation_linear_model.md)
  : Mediation analysis using linear models

- [`assoc_model_performance()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_model_performance.md)
  : Compute model performance metrics

- [`assoc_polynomial_formula_build()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_polynomial_formula_build.md)
  : Build a polynomial regression formula with covariate interactions

- [`assoc_quantile_permutation_model()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_quantile_permutation_model.md)
  : assoc_quantile_permutation_model calculate differences between the
  same quantile of two distribution

- [`assoc_quantreg_model()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_quantreg_model.md)
  : Quantile regression model (lqm)

- [`assoc_quantreg_permutation_model()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_quantreg_permutation_model.md)
  : Quantile regression permutation model

- [`assoc_spearman_permutation()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_spearman_permutation.md)
  : Spearman permutation model

- [`assoc_test_model()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_test_model.md)
  : Statistical test model dispatcher

- [`assoc_volcano_plot_inference()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_volcano_plot_inference.md)
  : Volcano plot of association results for one inference_detail row.

- [`association_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/association_analysis.md)
  : Association analysis of SEMseeker results

- [`core_check_session_compatibility()`](https://corsaro-lab.github.io/SEMseeker/reference/core_check_session_compatibility.md)
  : Check compatibility of SEMseeker sessions before meta-analysis

- [`core_get_meth_tech()`](https://corsaro-lab.github.io/SEMseeker/reference/core_get_meth_tech.md)
  : Detect the Illumina methylation array technology from a signal
  matrix

- [`core_init_env()`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md)
  : init ssEnvonment

- [`cytoband_hg19`](https://corsaro-lab.github.io/SEMseeker/reference/cytoband_hg19.md)
  : cytoband_hg19

- [`diagnostic_performance()`](https://corsaro-lab.github.io/SEMseeker/reference/diagnostic_performance.md)
  : Diagnostic performance (sensitivity & specificity) of SEMseeker
  mutations and lesions

- [`dmr_annotation`](https://corsaro-lab.github.io/SEMseeker/reference/dmr_annotation.md)
  : dmr_annotation

- [`enrich_find_unique_gene_sets()`](https://corsaro-lab.github.io/SEMseeker/reference/enrich_find_unique_gene_sets.md)
  : Find gene sets unique to each group

- [`enrich_wrap_it()`](https://corsaro-lab.github.io/SEMseeker/reference/enrich_wrap_it.md)
  : Wrap long strings to a fixed width

- [`enrichment_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/enrichment_analysis.md)
  : Enrichment analysis (pathway + phenotype)

- [`io_build_data_set_from_geo()`](https://corsaro-lab.github.io/SEMseeker/reference/io_build_data_set_from_geo.md)
  : io_build_data_set_from_geo

- [`io_coord_probe_features()`](https://corsaro-lab.github.io/SEMseeker/reference/io_coord_probe_features.md)
  : Build a minimal probe_features data frame from synthetic probe IDs.

- [`io_coord_to_semseeker()`](https://corsaro-lab.github.io/SEMseeker/reference/io_coord_to_semseeker.md)
  : Convert a coordinate-based methylation data frame to SEMseeker
  internal format.

- [`io_data_preparation()`](https://corsaro-lab.github.io/SEMseeker/reference/io_data_preparation.md)
  : Title

- [`io_dir_check_and_create()`](https://corsaro-lab.github.io/SEMseeker/reference/io_dir_check_and_create.md)
  : Create a directory path, building any missing intermediate
  directories

- [`io_dump_sample_as_bed_file()`](https://corsaro-lab.github.io/SEMseeker/reference/io_dump_sample_as_bed_file.md)
  : given data and colnames dump as bed file

- [`io_is_coord_format()`](https://corsaro-lab.github.io/SEMseeker/reference/io_is_coord_format.md)
  : Detect whether a data frame is in coordinate format (CHR + START
  columns).

- [`io_normalize_signal_input()`](https://corsaro-lab.github.io/SEMseeker/reference/io_normalize_signal_input.md)
  : Transparently normalise signal input for any technology.

- [`io_pivot_to_long_format()`](https://corsaro-lab.github.io/SEMseeker/reference/io_pivot_to_long_format.md)
  : Get the pivot in long format instead of wide format

- [`io_probe_id_to_coord()`](https://corsaro-lab.github.io/SEMseeker/reference/io_probe_id_to_coord.md)
  : Parse synthetic probe IDs back to a CHR / START / END data frame.

- [`meta_association_across_studies()`](https://corsaro-lab.github.io/SEMseeker/reference/meta_association_across_studies.md)
  : Cross-study meta-analysis of association results

- [`meta_association_overlaps_studies()`](https://corsaro-lab.github.io/SEMseeker/reference/meta_association_overlaps_studies.md)
  : Compare association results across studies

- [`meta_association_overlaps_subsamples()`](https://corsaro-lab.github.io/SEMseeker/reference/meta_association_overlaps_subsamples.md)
  : Compare association results across subsamples of one study

- [`meta_association_replication()`](https://corsaro-lab.github.io/SEMseeker/reference/meta_association_replication.md)
  : Intra-study association (focused replication across models)

- [`meta_enrichment_compare_studies()`](https://corsaro-lab.github.io/SEMseeker/reference/meta_enrichment_compare_studies.md)
  : Compare enrichment results across studies

- [`meta_enrichment_overlaps_subsamples()`](https://corsaro-lab.github.io/SEMseeker/reference/meta_enrichment_overlaps_subsamples.md)
  : Enrichment stability across random subsamples of a study

- [`metrics_properties`](https://corsaro-lab.github.io/SEMseeker/reference/metrics_properties.md)
  : metrics_properties

- [`plot_create_heatmap()`](https://corsaro-lab.github.io/SEMseeker/reference/plot_create_heatmap.md)
  : plot_create_heatmap load the multiple bed resulting from analysis
  organized into files and folders per marker and produce a pivot

- [`plot_manhattan_plot_per_area()`](https://corsaro-lab.github.io/SEMseeker/reference/plot_manhattan_plot_per_area.md)
  : Manhattan plot of association results per genomic area

- [`sem_analyze_population()`](https://corsaro-lab.github.io/SEMseeker/reference/sem_analyze_population.md)
  : Calculate stochastic epi mutations from a methylation dataset as
  outcome report of pivot

- [`sem_analyze_single_sample()`](https://corsaro-lab.github.io/SEMseeker/reference/sem_analyze_single_sample.md)
  : sem_analyze_single_sample

- [`sem_coverage_analysis_report()`](https://corsaro-lab.github.io/SEMseeker/reference/sem_coverage_analysis_report.md)
  : Probe coverage analysis report

- [`sem_delta_single_sample()`](https://corsaro-lab.github.io/SEMseeker/reference/sem_delta_single_sample.md)
  : sem_delta_single_sample

- [`sem_deltar_single_sample()`](https://corsaro-lab.github.io/SEMseeker/reference/sem_deltar_single_sample.md)
  : sem_deltar_single_sample

- [`sem_manhattan_plot_marker_per_sample()`](https://corsaro-lab.github.io/SEMseeker/reference/sem_manhattan_plot_marker_per_sample.md)
  : Manhattan plot of SEM markers for a single sample

- [`sem_metrics_filter()`](https://corsaro-lab.github.io/SEMseeker/reference/sem_metrics_filter.md)
  : Filter metrics by transformation type

- [`sem_mutations_get()`](https://corsaro-lab.github.io/SEMseeker/reference/sem_mutations_get.md)
  : sem_mutations_get

- [`sem_signal_range_values()`](https://corsaro-lab.github.io/SEMseeker/reference/sem_signal_range_values.md)
  : calculate the range of signal values to define the outlier

- [`sem_signal_single_sample()`](https://corsaro-lab.github.io/SEMseeker/reference/sem_signal_single_sample.md)
  : sem_signal_single_sample

- [`semseeker()`](https://corsaro-lab.github.io/SEMseeker/reference/semseeker.md)
  : Run SEMseeker on methylation data from any supported source

- [`ssEnv`](https://corsaro-lab.github.io/SEMseeker/reference/ssEnv.md)
  : ssEnv

- [`test_master_features`](https://corsaro-lab.github.io/SEMseeker/reference/test_master_features.md)
  : test_master_features

- [`test_samplesheet_gse133774`](https://corsaro-lab.github.io/SEMseeker/reference/test_samplesheet_gse133774.md)
  : test_samplesheet_gse133774

- [`test_signal_gse133774`](https://corsaro-lab.github.io/SEMseeker/reference/test_signal_gse133774.md)
  : test_signal_gse133774

- [`util_boolean_check()`](https://corsaro-lab.github.io/SEMseeker/reference/util_boolean_check.md)
  : Coerce a value to logical

- [`util_data_frame_add_column()`](https://corsaro-lab.github.io/SEMseeker/reference/util_data_frame_add_column.md)
  : Add or update a column in a data frame

- [`util_describe_dataframe()`](https://corsaro-lab.github.io/SEMseeker/reference/util_describe_dataframe.md)
  : Describe a data frame with summary statistics per column

- [`util_exploratory_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/util_exploratory_analysis.md)
  : Preliminary exploratory analysys of the data

- [`util_test_match_order()`](https://corsaro-lab.github.io/SEMseeker/reference/util_test_match_order.md)
  : Title
