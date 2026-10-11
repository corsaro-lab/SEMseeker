test_that("anno_create_position_pivots", {

  tempFolder <- sem_test_folder()
  unlink(tempFolder, recursive = TRUE, force = TRUE)
  # message(tempFolder)
  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median", start_fresh =TRUE)

  ####################################################################################
  tt <- SEMseeker:::core_get_meth_tech(signal_data)
  ####################################################################################

  ss_reference <- subset(mySampleSheet,Sample_Group=="Reference")
  # ss_control <- subset(mySampleSheet,Sample_Group=="Control")[1:3,]
  # ss_case <- subset(mySampleSheet,Sample_Group=="Case")[1:3,]
  # ss_reference <- ss_reference[1:3,]
  # mySampleSheet <- rbind(ss_case,ss_control,ss_reference)
  # signal_data <- signal_data[,unique(c(ss_case$Sample_ID,ss_control$Sample_ID,ss_reference$Sample_ID))]

  if (!exists("signal_thresholds"))
  {
    signal_data <- SEMseeker:::sem_inpute_missing_values(signal_data)
    signal_thresholds <<- SEMseeker:::sem_signal_range_values(signal_data[,ss_reference$Sample_ID],batch_id)
  }
  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]

  sp <- SEMseeker:::sem_analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Case",],
    probe_features = probe_features
  )

  sp <- SEMseeker:::sem_analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Control",],
    probe_features = probe_features
  )

  ssEnv <- SEMseeker:::core_get_session_info()

  keys <- expand.grid(
    MARKER = c("DELTAS","DELTAR","LESIONS","MUTATIONS"),
    FIGURE = c("HYPER","HYPO")
  )
  area <- "POSITION"
  # every coordinate is explicit. The registry gives POSITION its
  # subarea (WHOLE) like every other area; "" was this test writing a hole by
  # hand, and the key composer now refuses it instead of filling it in.
  subarea <- "WHOLE"

  # prova con un subset di colonne che non ha le lesioni

  # test parial pivot creation
  # SEMseeker:::anno_create_position_pivots(mySampleSheet[1:10,],keys)

  # test complete
  SEMseeker:::anno_create_position_pivots(mySampleSheet,keys)

  # for (k in 1:nrow(keys))
  {
    k <- 1
    key <- keys[k,]
    marker <- as.character(key$MARKER)
    figure <- as.character(key$FIGURE)

    mutations_pivot_file_name <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS",figure,area,subarea)
    if(file.exists(mutations_pivot_file_name))
    {
      mutations_pivot <- as.data.frame(polars::pl$read_parquet(mutations_pivot_file_name))

      io_pivot_file_name <- SEMseeker:::io_pivot_file_name_parquet(marker,figure,area,subarea)
      testthat::expect_true(file.exists(io_pivot_file_name))
      if(file.exists(io_pivot_file_name))
        pivot <- as.data.frame(polars::pl$read_parquet(io_pivot_file_name))
      else
        pivot <- NULL

      pivot <- pivot[,-c(1:3)]
      mutations_pivot <- mutations_pivot[,-c(1:3)]
      testthat::expect_true(nrow(pivot)<nprobes)
      testthat::expect_true(nrow(pivot)>0)

      if (marker!="LESIONS")
      {
        testthat::expect_true(all(colnames(pivot) %in% (mySampleSheet$Sample_ID)))
        testthat::expect_true(all(colnames(pivot) == colnames(mutations_pivot)))

        testthat::expect_true(ncol(pivot)==ncol(mutations_pivot))
        testthat::expect_true(nrow(pivot)==nrow(mutations_pivot))

        pivot[is.na(pivot)] <- 0
        mutations_pivot[is.na(mutations_pivot)] <- 0

        pivot[pivot > 0] <- 1
        mutations_pivot[mutations_pivot > 0] <- 1

        testthat::expect_true(all(as.data.frame(pivot) == as.data.frame(mutations_pivot)))
      }
    }
  }

  ####################################################################################

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})
