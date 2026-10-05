util_keys_create <- function(ssEnv, arguments)
{

  # remove BOTH and BOTHSUM
  ssEnv$default <- list()
  ssEnv$keys_sample_groups <-  data.frame("SAMPLE_GROUP"=c("Reference","Control","Case"))
  ##### MANAGE MARKERS AND FIGURES #####
  # set default values
  keys_figures_default_discrete <-  data.frame("FIGURE"=c("HYPO", "HYPER", "BOTH","BOTHSUM")) #
  keys_markers_default_discrete <-  data.frame("MARKER"=c("MUTATIONS","LESIONS","DELTAQ","DELTARQ","DELTAP","DELTARP"))


  # AI-092: LESIONS uses LESIONS_BP (max bp distance for two probes to be in
  # the same enrichment window, default 2000 = typical CpG-island / DMR span)
  # instead of the legacy sliding_window_size (probe-count based, no longer
  # supported). Same registration pattern as DELTAQ_Q / DELTAP_B.
  keys_markers_default_discrete$Q <- c(1,ssEnv$LESIONS_BP,ssEnv$DELTAQ_Q,ssEnv$DELTARQ_Q,ssEnv$DELTAP_B,ssEnv$DELTARP_B)
  keys_markers_default_discrete$SOURCE <- c("DELTAS","MUTATIONS","DELTAS","DELTAR","DELTAS","DELTAR")
  keys_markers_default_discrete$SUFFIX <-  rep("",6)
  # keys_markers_default$SUFFIX <-  c("","","","", ssEnv$epiquantile ,"",ssEnv$epiquantile)
  keys_markers_default_discrete$EXT <-  c("bed","bed","bed","bed","bed","bed")
  keys_markers_figures_default_discrete <- merge(keys_markers_default_discrete, keys_figures_default_discrete, by = NULL)
  rm(keys_figures_default_discrete, keys_markers_default_discrete)


  keys_figures_default_continuos <-  data.frame("FIGURE"=c("HYPO", "HYPER", "BOTH","BOTHSUM")) # , "BOTH","BOTHSUM"
  keys_markers_default_continuos <-  data.frame("MARKER"=c("DELTAS","DELTAR"))
  keys_markers_default_continuos$SUFFIX <-  c("","")
  # keys_markers_default$SUFFIX <-  c("","","","", ssEnv$epiquantile ,"",ssEnv$epiquantile)
  keys_markers_default_continuos$EXT <-  c("bedgraph","bedgraph")
  keys_markers_figures_default_continuos <- merge(keys_markers_default_continuos, keys_figures_default_continuos, by = NULL)
  rm(keys_figures_default_continuos, keys_markers_default_continuos)

  # AI-248: the figure of SIGNAL is the SCALE of the value (BETA / MVALUE), not
  # a way of aggregating it - "MEAN" used to sit here as a placeholder that made
  # the key unique. The scale is unknown at init (no data read yet), so this is
  # the bounded default; core_get_meth_tech() realigns it via
  # util_keys_signal_figure_refresh() as soon as the first batch is seen.
  keys_figures_default_continuos_mono_figure <-  data.frame("FIGURE"=c(io_signal_figure()))
  keys_markers_default_continuos_mono_figure <-  data.frame("MARKER"=c("SIGNAL"))
  keys_markers_default_continuos_mono_figure$SUFFIX <-  c("")
  # keys_markers_default$SUFFIX <-  c("","","","", ssEnv$epiquantile ,"",ssEnv$epiquantile)
  keys_markers_default_continuos_mono_figure$EXT <-  c("bedgraph")
  keys_markers_figures_default_continuos_mono_figure <- merge(keys_markers_default_continuos_mono_figure, keys_figures_default_continuos_mono_figure, by = NULL)
  rm(keys_figures_default_continuos_mono_figure, keys_markers_default_continuos_mono_figure)

  keys_markers_figures_default <- plyr::rbind.fill(keys_markers_figures_default_discrete,
    keys_markers_figures_default_continuos,
    keys_markers_figures_default_continuos_mono_figure)
  keys_markers_figures_default$DISCRETE <- keys_markers_figures_default$EXT=="bed"

  # remove BOTH and BOTHSUM
  keys_markers_figures_default <- keys_markers_figures_default[!(keys_markers_figures_default$FIGURE=="BOTH" | keys_markers_figures_default$FIGURE=="BOTHSUM"),]

  ssEnv$keys_markers_figures_default <- keys_markers_figures_default
  ssEnv$default$keys_markers_figures_default <- keys_markers_figures_default

  figures <- arguments[["figures"]]
  unknown_figures <- figures[!figures %in% keys_markers_figures_default$FIGURE]
  arguments[["figures"]] <- unknown_figures
  if (!is.null(figures))
    keys_markers_figures_default <- keys_markers_figures_default[keys_markers_figures_default$FIGURE %in% figures,]

  markers <- arguments[["markers"]]
  unknown_markers <- markers[!markers %in% keys_markers_figures_default$MARKER]
  arguments[["markers"]] <- unknown_markers
  if (!is.null(markers))
    keys_markers_figures_default <- keys_markers_figures_default[keys_markers_figures_default$MARKER %in% markers,]

  ########################################################################################################################
  # manage AREAS and SUBAREAS #####

  keys_areas_default <- data.frame("AREA"=c("GENE","ISLAND","DMR","CHR","PROBE"))

  # PROMOTER is a window of its own here, not a shorthand read back later: the
  # positions of TSS200, TSS1500 and 1STEXON are masked once, so one gene gets one
  # measure. Grouped at read time it would reach three families of results, and a
  # gene with three p-values has none.
  keys_gene_subareas_default <- data.frame("AREA"="GENE", "SUBAREA"=c("BODY","TSS1500","TSS200","1STEXON","3UTR","5UTR","EXONBND","PROMOTER","WHOLE"))
  # ISLAND context (Illumina Relation_to_Island, 6 categories): WHOLE = whole
  # neighbourhood (core+shores+shelves, like GENE_WHOLE), ISLAND = core only,
  # OPENSEA = CpGs outside every neighbourhood. See island_opensea.R.
  keys_island_subareas_default <- data.frame("AREA"="ISLAND","SUBAREA"=c("ISLAND","N_SHORE","S_SHORE","N_SHELF","S_SHELF","OPENSEA","WHOLE"))
  keys_dmr_subareas_default <- data.frame("AREA"="DMR","SUBAREA"=c("WHOLE"))
  keys_chr_subareas_default <- data.frame("AREA"="CHR","SUBAREA"=c("WHOLE","CYTOBAND"))
  keys_probe_subareas_default <- data.frame("AREA"="PROBE","SUBAREA"=c("WHOLE"))
  keys_position_subareas_default <- data.frame("AREA"="POSITION","SUBAREA"=c("WHOLE"))

  keys_areas_subareas_default <- rbind(keys_gene_subareas_default, keys_island_subareas_default, keys_dmr_subareas_default, keys_chr_subareas_default, keys_probe_subareas_default, keys_position_subareas_default)
  ssEnv$default$keys_areas_subareas_default <- keys_areas_subareas_default

  areas <- unique(c(arguments[["areas"]],"POSITION"))
  if ("ALL" %in% areas)
  {
    keys_areas_subareas_default <- keys_areas_subareas_default
    arguments[["areas"]] <- NULL
  }
  else
  {
    unknown_areas <- areas[!areas %in% keys_areas_subareas_default$AREA]
    arguments[["areas"]] <- unknown_areas
    if (!is.null(areas))
      keys_areas_subareas_default <- keys_areas_subareas_default[keys_areas_subareas_default$AREA %in% areas,]
  }


  subareas <- unique(c(arguments[["subareas"]]))
  if("ALL" %in% subareas)
  {
    keys_areas_subareas_default <- keys_areas_subareas_default
    arguments[["subareas"]] <- NULL
  }
  else
  {
    unknown_subareas <- subareas[!subareas %in% keys_areas_subareas_default$SUBAREA]
    # could exist WHOLE associated to AREAS where exists empty SUBAREAS
    arguments[["subareas"]] <- unknown_subareas[!(unknown_subareas %in% c("WHOLE"))]
    if (!is.null(subareas))
      keys_areas_subareas_default <- keys_areas_subareas_default[keys_areas_subareas_default$SUBAREA %in% subareas,]
  }

  ########################################################################################################################
  # combine AREAS, SUBAREAS, MARKERS and FIGURES #####

  keys_areas_subareas_markers_figures <- unique(merge(keys_areas_subareas_default, keys_markers_figures_default, by = NULL))

  ########################################################################################################################

  # assign to ssEnv
  ssEnv$keys_areas_subareas <- keys_areas_subareas_default
  ssEnv$keys_markers_figures <- keys_markers_figures_default
  ssEnv$keys_areas_subareas_markers_figures <- keys_areas_subareas_markers_figures

  combine_not_empty <- function(x)
  {
    # paste0(x[x!=""], collapse = "_")
    # remove spaces
    core_name_cleaning(gsub(" ","",paste0(x[x!=""], collapse = "_")))
  }

  ssEnv$keys_areas_subareas_markers_figures$COMBINED <- apply(ssEnv$keys_areas_subareas_markers_figures[,c("MARKER","FIGURE","AREA","SUBAREA")], 1, combine_not_empty )
  ssEnv$keys_areas_subareas$COMBINED <- apply(ssEnv$keys_areas_subareas[,c("AREA","SUBAREA")], 1, combine_not_empty )
  ssEnv$keys_markers_figures$COMBINED <- apply(ssEnv$keys_markers_figures[,c("MARKER","FIGURE")], 1, combine_not_empty )

  ########################################################################################################################
  # prepare keys for pathway analysis
  keys <- ssEnv$keys_areas_subareas_markers_figures
  # select only gene
  keys <- keys[keys$AREA=="GENE",]
  keys$FIGURE <- as.character(keys$FIGURE)
  keys$SUBAREA <- as.character(keys$SUBAREA)
  # remove "COMBINED"
  keys <- keys[,c("AREA","SUBAREA","MARKER","FIGURE")]
  # where FIGURE is HYPO or HYPER for FIGURE set to HYPER_HYPO
  selector <- (keys$FIGURE=="HYPER" | keys$FIGURE=="HYPO")
  if (any(selector))
    keys[selector,"FIGURE"] <- "HYPER_HYPO"
  keys <- unique(keys)
  parts <- ssEnv$keys_areas_subareas[ssEnv$keys_areas_subareas[,"SUBAREA"]!="WHOLE","SUBAREA"]
  keys_parts <- keys[keys$SUBAREA %in% parts,]
  if (nrow(keys_parts)>0)
    keys_parts$SUBAREA <- "ALL_SUBAREAS"
  keys_whole <- keys[keys$SUBAREA=="WHOLE",]
  keys <- rbind(keys_parts,keys_whole)
  keys <- keys[!duplicated(keys),]
  # remove only HYPER or only HYPO
  keys <- keys[!(keys$FIGURE=="HYPER" | keys$FIGURE=="HYPO"),]
  if(nrow(keys)!=0)
    keys$COMBINED <- paste0(keys$AREA,keys$SUBAREA,keys$MARKER,keys$FIGURE, sep="_")
  # remove BOTH and BOTHSUM
  keys <- keys[!(keys$FIGURE=="BOTH" | keys$FIGURE=="BOTHSUM"),]
  ssEnv$keys_for_pathway <- keys

  ########################################################################################################################
  # prepare keys for enrichment analysis report format
  key_enrichment_format <- data.frame("type"="Pathway", "label"="STRINGdb","column_of_id"="term","column_of_description"="description", "column_of_adj_pvalue"="fdr","column_of_enrichment"="fold_enrichment","column_of_pvalue"="p_value", "column_of_genes"="")
  key_enrichment_format <- rbind(key_enrichment_format,data.frame("type"="Pathway", "label"="Phenolyzer_STRINGdb","column_of_id"="term","column_of_description"="description", "column_of_adj_pvalue"="fdr","column_of_enrichment"="fold_enrichment","column_of_pvalue"="p_value", "column_of_genes"=""))
  key_enrichment_format <- rbind(key_enrichment_format, data.frame("type"="Pathway","label"="WebGestalt","column_of_id"="geneSet","column_of_description"="description","column_of_adj_pvalue"="FDR","column_of_enrichment"="enrichmentRatio","column_of_pvalue"="pValue", "column_of_genes"=""))
  key_enrichment_format <- rbind(key_enrichment_format, data.frame("type"="Pathway","label"="Phenolyzer_WebGestalt","column_of_id"="geneSet","column_of_description"="description","column_of_adj_pvalue"="FDR","column_of_enrichment"="enrichmentRatio","column_of_pvalue"="pValue", "column_of_genes"=""))
  key_enrichment_format <- rbind(key_enrichment_format, data.frame("type"="Phenotype","label"="phenolyzer","column_of_id"="Description","column_of_description"="Description","column_of_adj_pvalue"="Score","column_of_enrichment"="Score","column_of_pvalue"="", "column_of_genes"=""))
  key_enrichment_format <- rbind(key_enrichment_format,data.frame("type"="Pathway", "label"="pathfindR","column_of_id"="ID","column_of_description"="Term_Description", "column_of_adj_pvalue"="highest_p","column_of_enrichment"="Fold_Enrichment", "column_of_pvalue"="highest_p", "column_of_genes"="Up_regulated+Down_regulated"))
  key_enrichment_format <- rbind(key_enrichment_format,data.frame("type"="Pathway", "label"="ctdR","column_of_id"="ChemicalID","column_of_description"="ChemicalName", "column_of_adj_pvalue"="padj","column_of_enrichment"="foldEnrichment","column_of_pvalue"="pvalue", "column_of_genes"=""))

  ssEnv$key_enrichment_format <- key_enrichment_format

  ssEnv$key_missed_areas_subareas <- ssEnv$keys_areas_subareas[0,]

  core_update_session_info(ssEnv)

  return(arguments)
}
