# 2026-06-09: smart split of slash-separated AREA names with
# prefix recovery. See R/smart_split_area_name.R for the heuristics.

test_that(".anno_smart_split_area_name returns name unchanged when no '/'", {
  expect_equal(SEMseeker:::.anno_smart_split_area_name("TP53"),               "TP53")
  expect_equal(SEMseeker:::.anno_smart_split_area_name("HLA-A"),              "HLA-A")
  expect_equal(SEMseeker:::.anno_smart_split_area_name("ANKHD1-EIF4EBP3"),    "ANKHD1-EIF4EBP3")
  expect_equal(SEMseeker:::.anno_smart_split_area_name(""),                   "")
})

test_that("Strategy 1: dash prefix recovery (HLA family pattern)", {
  expect_equal(
    SEMseeker:::.anno_smart_split_area_name("HLA-A/B/C"),
    c("HLA-A", "HLA-B", "HLA-C")
  )
  expect_equal(
    SEMseeker:::.anno_smart_split_area_name("GENE-A-B/C/D"),
    c("GENE-A-B", "GENE-A-C", "GENE-A-D")
  )
})

test_that("Strategy 2: alphabetic-leading prefix recovery (KRT pattern)", {
  expect_equal(
    SEMseeker:::.anno_smart_split_area_name("KRT8/18"),
    c("KRT8", "KRT18")
  )
  expect_equal(
    SEMseeker:::.anno_smart_split_area_name("ATP5/F/D"),
    c("ATP5", "ATPF", "ATPD")
  )
})

test_that("No double-prepending when suffix already starts with prefix", {
  # HBA1/HBA2: prefix = "HBA", but "HBA2" already starts with "HBA"
  expect_equal(
    SEMseeker:::.anno_smart_split_area_name("HBA1/HBA2"),
    c("HBA1", "HBA2")
  )
  # HLA-A/HLA-B: prefix = "HLA-", "HLA-B" already starts with "HLA-"
  expect_equal(
    SEMseeker:::.anno_smart_split_area_name("HLA-A/HLA-B"),
    c("HLA-A", "HLA-B")
  )
})

test_that("NA and zero-length inputs are passed through", {
  expect_true(is.na(SEMseeker:::.anno_smart_split_area_name(NA_character_)))
  expect_equal(SEMseeker:::.anno_smart_split_area_name(""), "")
})

test_that("Falls back to raw split when no prefix is recoverable", {
  # Pure-digits leading: no alpha prefix to recover. Fallback = raw parts.
  expect_equal(
    SEMseeker:::.anno_smart_split_area_name("123/456"),
    c("123", "456")
  )
})
