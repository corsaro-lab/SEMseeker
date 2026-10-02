# Materialise the area-level artefacts of the run Rewritten as a thin batch pass over \[io_pivot_build()\], which is now the single derivation path of the taxonomy. What used to live here - the join with the annotation, the multi-gene explode, the group-by - moved into \`io_pivot_build()\` and \[.anno_area_explode()\], so the batch pre-pass and the on-demand construction cannot diverge. What disappeared with it is the rule that chose the operator:

     if (localKeys[i, "DISCRETE"]) pivot$group_by("AREA")$sum() else pivot$group_by("AREA")$mean() 

One file per key, the operator picked by a flag and \*\*absent from the
name\*\*. Downstream, \`inference_details\$aggregation\` was validated
and then ignored: asking for \`MEDIAN\` on \`GENE_TSS1500\` returned the
mean, silently. The keys now carry the aggregation, so what is written
says which operator produced it. This pass is a convenience, not a
requirement: anything it does not materialise is built on first read by
\[io_read_pivot()\].

Materialise the area-level artefacts of the run

Rewritten as a thin batch pass over \[io_pivot_build()\], which is now
the single derivation path of the taxonomy. What used to live here - the
join with the annotation, the multi-gene explode, the group-by - moved
into \`io_pivot_build()\` and \[.anno_area_explode()\], so the batch
pre-pass and the on-demand construction cannot diverge.

What disappeared with it is the rule that chose the operator:


    if (localKeys[i, "DISCRETE"]) pivot$group_by("AREA")$sum()
    else                          pivot$group_by("AREA")$mean()

One file per key, the operator picked by a flag and \*\*absent from the
name\*\*. Downstream, \`inference_details\$aggregation\` was validated
and then ignored: asking for \`MEDIAN\` on \`GENE_TSS1500\` returned the
mean, silently. The keys now carry the aggregation, so what is written
says which operator produced it.

This pass is a convenience, not a requirement: anything it does not
materialise is built on first read by \[io_read_pivot()\].

## Usage

``` r
anno_annotate_position_pivots()
```

## Value

nothing
