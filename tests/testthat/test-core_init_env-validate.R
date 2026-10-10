# Tests for .core_init_env_validate_args() diagnostic error message.
# These tests pin the user-facing wording of the error so a regression
# would be caught immediately.

test_that(".core_init_env_validate_args is a no-op when there are no leftover args", {
  expect_no_error(SEMseeker:::.core_init_env_validate_args(list()))
})

test_that(".core_init_env_validate_args ignores NULL and empty-character args", {
  args <- list(foo = NULL, bar = character(0))
  expect_no_error(SEMseeker:::.core_init_env_validate_args(args))
})

test_that(".core_init_env_validate_args errors on an unrecognised named arg", {
  args <- list(phenolyser = FALSE)
  expect_error(
    SEMseeker:::.core_init_env_validate_args(args),
    "Unrecognised argument"
  )
})

test_that(".core_init_env_validate_args shows the argument NAME (not just value)", {
  # Pre-the error said "ERROR: This options are not recognized: FALSE"
  # (only the value), which made the typo invisible. Now it must contain the name.
  args <- list(phenolyser = FALSE)
  err <- tryCatch(
    SEMseeker:::.core_init_env_validate_args(args),
    error = identity
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "phenolyser")
})

test_that(".core_init_env_validate_args mentions all unrecognised names", {
  args <- list(foo = "x", bar = 1L)
  err <- tryCatch(
    SEMseeker:::.core_init_env_validate_args(args),
    error = identity
  )
  expect_match(conditionMessage(err), "foo")
  expect_match(conditionMessage(err), "bar")
})

test_that(".core_init_env_validate_args truncates very long argument values", {
  args <- list(myparam = paste(rep("X", 200), collapse = ""))
  err <- tryCatch(
    SEMseeker:::.core_init_env_validate_args(args),
    error = identity
  )
  expect_match(conditionMessage(err), "myparam")
  # full 200-X string MUST NOT appear verbatim
  expect_false(grepl(paste(rep("X", 200), collapse = ""), conditionMessage(err), fixed = TRUE))
  # truncation marker should be present
  expect_match(conditionMessage(err), "\\.\\.\\.")
})
