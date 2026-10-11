test_that("io_is_coord_format detects CHR/START columns", {
  df_coord <- data.frame(
    CHR = c("chr1", "chr2"), START = c(1000L, 2000L), END = c(1001L, 2001L),
    s1 = c(0.8, 0.3), s2 = c(0.5, 0.9)
  )
  expect_true(SEMseeker:::io_is_coord_format(df_coord))

  df_illumina <- data.frame(s1 = c(0.8, 0.3), s2 = c(0.5, 0.9),
                              row.names = c("cg00000029", "cg00000165"))
  expect_false(SEMseeker:::io_is_coord_format(df_illumina))
})

test_that("io_is_coord_format accepts lowercase and alias column names", {
  expect_true(SEMseeker:::io_is_coord_format(data.frame(chr = "1", start = 100, s1 = 0.5)))
  expect_true(SEMseeker:::io_is_coord_format(data.frame(chrom = "X", chromStart = 500, s1 = 0.5)))
  expect_true(SEMseeker:::io_is_coord_format(data.frame(seqnames = "1", pos = 100, s1 = 0.5)))
})

test_that("io_coord_to_semseeker converts wide coord df to probe-indexed matrix", {
  df <- data.frame(
    CHR   = c("chr1", "chrX", "chr22"),
    START = c(10000L, 543200L, 900100L),
    END   = c(10001L, 543201L, 900101L),
    s1    = c(0.8, 0.5, 0.2),
    s2    = c(0.3, 0.7, 0.9),
    stringsAsFactors = FALSE
  )
  result <- SEMseeker:::io_coord_to_semseeker(df)

  expect_equal(nrow(result), 3L)
  expect_equal(ncol(result), 2L)  # only s1, s2
  expect_equal(rownames(result), c("1_10000", "X_543200", "22_900100"))
  expect_false("CHR" %in% colnames(result))
  expect_false("START" %in% colnames(result))
  expect_false("END" %in% colnames(result))
})

test_that("io_coord_to_semseeker strips chr prefix consistently", {
  df <- data.frame(CHR = c("CHR1", "Chr2", "chr3"),
                   START = c(100L, 200L, 300L),
                   s1 = c(0.5, 0.5, 0.5))
  res <- SEMseeker:::io_coord_to_semseeker(df)
  expect_equal(rownames(res), c("1_100", "2_200", "3_300"))
})

test_that("io_coord_to_semseeker works without END column", {
  df <- data.frame(CHR = "chr1", START = 5000L, s1 = 0.6)
  res <- SEMseeker:::io_coord_to_semseeker(df)
  expect_equal(rownames(res), "1_5000")
  expect_equal(ncol(res), 1L)
})

test_that("io_probe_id_to_coord round-trips with io_coord_to_semseeker", {
  df <- data.frame(
    CHR = c("chr1","chrX","chr22"), START = c(100L, 200L, 300L),
    s1 = c(0.1, 0.2, 0.3)
  )
  probe_ids <- rownames(SEMseeker:::io_coord_to_semseeker(df))
  back      <- SEMseeker:::io_probe_id_to_coord(probe_ids)

  expect_equal(back$CHR,   c("1", "X", "22"))
  expect_equal(back$START, c(100L, 200L, 300L))
  expect_equal(back$END,   c(101L, 201L, 301L))
})

test_that("io_coord_probe_features builds minimal probe_features df", {
  ids <- c("1_10000", "X_543200", "22_900100")
  pf  <- SEMseeker:::io_coord_probe_features(ids)

  expect_equal(nrow(pf), 3L)
  expect_true(all(c("PROBE", "CHR", "START", "END") %in% colnames(pf)))
  expect_equal(pf$PROBE, ids)
  expect_equal(pf$CHR,   c("1", "X", "22"))
  expect_equal(pf$START, c(10000L, 543200L, 900100L))
  expect_equal(pf$END,   c(10001L, 543201L, 900101L))
})

test_that("io_normalize_signal_input passes through Illumina matrix unchanged", {
  df <- data.frame(s1 = c(0.8, 0.3), row.names = c("cg00000029", "cg00000165"))
  result <- SEMseeker:::io_normalize_signal_input(df)
  expect_identical(result, df)
})

test_that("io_normalize_signal_input converts coordinate data frame", {
  df <- data.frame(CHR = "chr1", START = 10000L, s1 = 0.8)
  result <- SEMseeker:::io_normalize_signal_input(df)
  expect_equal(rownames(result), "1_10000")
  expect_false(SEMseeker:::io_is_coord_format(result))
})
