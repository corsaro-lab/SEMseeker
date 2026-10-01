# No session folder here on purpose.
# and these tests need no session at all — they exercise pure vocabulary
# functions. Re-assigning it truncated the vector for every file that runs
# after this one alphabetically, which is how a change confined to the
# enrichment vocabulary took down semseeker(), the association suite and the
# burden integration in CI.

test_that("neither question has a default, and the refusal is the question", {
  # AI-261. Until this release the package answered both on the caller's behalf
  # in util_keys_create() and said so nowhere. A default here is the way a
  # question never gets asked, so there is none, and the error explains what the
  # answer changes rather than listing admissible strings.
  expect_error(SEMseeker:::util_gene_region_expand(),
               "gene_region has no default")
  expect_error(SEMseeker:::util_gene_region_expand(NULL), "no default")
  expect_error(SEMseeker:::util_epimutation_direction_expand(),
               "epimutation_direction has no default")

  # The message has to carry the biology, not just the vocabulary: this is the
  # only place a caller is told that the sign is read together with the region.
  msg <- tryCatch(SEMseeker:::util_epimutation_direction_expand(),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "silences")
  expect_match(msg, "de-represses")
  expect_match(msg, "body of a gene")
  expect_match(msg, "p-value per gene is NOT defined")

  region_msg <- tryCatch(SEMseeker:::util_gene_region_expand(),
                         error = function(e) conditionMessage(e))
  # The message names the alias and the biology behind it. It no longer names the
  # windows, because the promoter is one region class now and not three of them.
  expect_match(region_msg, "PROMOTER")
  expect_match(region_msg, "gates")
})

test_that("an alias expands to the windows it is documented to mean", {
  # One window, not three: the promoter is masked as a class while the analysis
  # runs, so there is nothing left to group on the way out. Grouping here put one
  # gene into three families, which is how a gene ended up with no p-value at all.
  expect_equal(as.character(SEMseeker:::util_gene_region_expand("PROMOTER")),
               "PROMOTER")
  expect_equal(as.character(SEMseeker:::util_gene_region_expand("GENE_BODY")),
               "BODY")
  expect_equal(as.character(SEMseeker:::util_gene_region_expand("WHOLE_GENE")),
               "WHOLE")

  # The alias travels with the expansion, so the result can be stamped with the
  # answer that produced it rather than only with its consequence.
  expect_equal(attr(SEMseeker:::util_gene_region_expand("PROMOTER"), "alias"),
               "PROMOTER")
})

test_that("explicit windows override the package's definition of a region", {
  # The point of allowing the alias at all: a caller who reads PROMOTER
  # differently is not stuck with ours.
  got <- SEMseeker:::util_gene_region_expand(c("TSS200", "5UTR"))
  expect_equal(as.character(got), c("TSS200", "5UTR"))
  expect_equal(attr(got, "alias"), "TSS200+5UTR")
})

test_that("an alias mixed with explicit windows is refused, not merged", {
  # Merging them would produce a result whose stamp cannot say which definition
  # of the region was used.
  expect_error(SEMseeker:::util_gene_region_expand(c("PROMOTER", "BODY")),
               "mixes an alias with explicit windows")
})

test_that("ANY is the union of the two directions, and declares itself as one", {
  hyper <- SEMseeker:::util_epimutation_direction_expand("HYPER")
  expect_equal(as.character(hyper), "HYPER")
  expect_false(attr(hyper, "union"))

  any_dir <- SEMseeker:::util_epimutation_direction_expand("ANY")
  expect_setequal(as.character(any_dir), c("HYPER", "HYPO"))
  expect_true(attr(any_dir, "union"))

  # The union flag is what a backend needing a p-value per gene reads to refuse,
  # so it has to be on the answer itself and not re-derived by counting.
  expect_false(attr(SEMseeker:::util_epimutation_direction_expand("HYPO"), "union"))
})

test_that("a direction outside the vocabulary is refused by name", {
  expect_error(SEMseeker:::util_epimutation_direction_expand("BOTH"),
               "is not one of HYPER, HYPO, ANY")
})
