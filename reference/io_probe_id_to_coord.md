# Parse synthetic probe IDs back to a CHR / START / END data frame.

Synthetic probe ID format: "\<CHR\>\_\<START\>" where CHR has no "chr"
prefix. E.g. "1_10000" → CHR = "1", START = 10000L, END = 10001L.

## Usage

``` r
io_probe_id_to_coord(probe_ids)
```

## Arguments

- probe_ids:

  Character vector of synthetic probe IDs.

## Value

data.frame with columns CHR (character), START (integer), END (integer).
