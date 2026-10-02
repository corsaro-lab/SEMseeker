# Cross-study meta-analysis of association results

Pools the association results of several studies: for each region class
(marker, figure, area, subarea) and each instance within it, it combines
the per-study effect sizes and their standard errors and reports
fixed-effect and random-effect estimates with heterogeneity statistics.

## Usage

``` r
meta_association_across_studies(
  studies,
  result_folder,
  inference_detail,
  markers = NULL,
  figures = NULL,
  areas = NULL,
  subareas = NULL,
  statistic_parameter = "BETA",
  pvalue_column = "PVALUE_ADJ_ALL_BH",
  alpha = 0.05,
  adjustment_method = "BH",
  ...
)
```

## Arguments

- studies:

  Character vector of study result paths, optionally named. The name, or
  else the final path component, labels the study.

- result_folder:

  Directory the meta-analysis writes to.

- inference_detail:

  One row of the inference specification, identifying the request whose
  per-study results are pooled.

- markers, figures, areas, subareas:

  Optional selection within the space of the meta-analysis. `NULL` takes
  the whole space.

- statistic_parameter:

  Name of the effect-size column to pool.

- pvalue_column:

  Name of the adjusted p-value column carried alongside.

- alpha:

  Significance threshold passed to the per-study readers.

- adjustment_method:

  Multiple-testing correction passed to the per-study readers.

- ...:

  Passed to the per-study session setup.

## Value

Nothing: the call stops. When implemented, one row per region class and
instance with the pooled estimate, its confidence intervals, the
p-values and the heterogeneity statistics, plus the coordinates that
identify the class the row came from.

## Not implemented

This function refuses every call. It is exported and named because its
contract is settled, but its body was never executed once: the step that
reads the per-study results was missing, so every call died on an
unbound variable before reaching the model. Rather than leave a function
that fails on a missing object, or guess a body that cannot be run, it
refuses at the door and says so.

The contract it will honour, so that a caller can be written against it:

- `studies` is a vector of paths, one per study. The label of a study is
  the element's name when it has one, otherwise the final component of
  its path. Studies are therefore not required to sit under a common
  parent, and there is no second argument that can disagree with the
  first. Two studies whose paths end in the same component are refused:
  they would become one label and merge silently.

- `result_folder` is where the meta-analysis writes. It is the opposite
  direction from the study paths, which are read.

- the space of the meta-analysis is the **strict intersection** of the
  region classes the studies have in common, read from each study's own
  session and with no filter applied. Classes left out of the
  intersection are reported, with how many and which studies lacked
  them.

- `markers`, `figures`, `areas` and `subareas` are the caller's
  selection, checked against that space one coordinate at a time: a
  value no study measured is refused and named. The classes analysed are
  the selection intersected with the space, and the number of
  combinations that fell outside is reported. The selection is
  deliberately not passed on to the per-study sessions: filtering them
  first would compute the intersection over the narrowed space, and a
  class missing from one study would vanish from the space instead of
  falling outside it.

- estimates are pooled within one region class at a time. Combining
  across markers or areas would average effects that are not the same
  effect and return something with the shape of a meta-analysis.

## Examples

``` r
# This function refuses every call; see the "Not implemented" section for the
# contract it will honour.
invisible(NULL)
```
