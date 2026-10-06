# Names R CMD check cannot bind, declared so that its "no visible binding"
# NOTE lists only what is really undefined.
#
# dplyr and ggplot2 column names are not declared here: they use the .data
# pronoun (.data$col, ggplot2::aes(.data$col)). What is declared is of two
# kinds, and each name says where it comes from.
#
# 1. Forward references guarded by exists(): the variable may or may not exist
#    at runtime.
# 2. Column names read through non-standard evaluation where the .data pronoun
#    does not apply: base subset(), the iterator of a foreach() call, and the
#    count column that dplyr::count() creates and tidyr::pivot_wider() reads.
#
# A name goes here only after its every use has been read and found to be one
# of these. A name that is simply never assigned is a defect, and declaring it
# would hide it from the one check that reports it.
utils::globalVariables(c(
  # 1. forward references
  "pathway_report",                         # pathway helpers, read from disk

  # 2. columns in subset()
  "AREA",                                   # assoc_results_get
  "FAMILY_TEST",                            # assoc_analysis_save_results
  "FIGURE",                                 # assoc_data_extractor
  "KEY",                                    # meta_association_overlaps_subsamples
  "MARKER",                                 # anno_create_position_pivots, assoc_run_marker, sem_deltaX_get
  "P_to_be_Case_cond_to_be_Epimutated",     # assoc_bayes_analysis
  "P_to_be_Control_cond_to_be_Epimutated",  # assoc_bayes_analysis
  "Q",                                      # sem_deltaX_get
  "Sample_Group",                           # assoc_bayes_analysis
  "SCOPE",                                  # assoc_results_get, meta_association_overlaps_subsamples
  "SIGNIFICATIVE_ADJ_ALL",                  # assoc_analysis_save_results
  "SOURCE",                                 # sem_deltaX_get
  "STUDY",                                  # meta_enrichment_compare_studies

  #    the iterator of foreach()
  "col_idx",                                # assoc_bayes_analysis

  #    dplyr::count() and tidyr::pivot_wider()
  "KEY_SELECTOR",                           # meta_enrichment_overlaps_subsamples
  "n"                                       # meta_enrichment_overlaps_subsamples
))
