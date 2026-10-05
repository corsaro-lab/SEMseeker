#' The circos layout, once: an ideogram, outside labels, and links between pairs
#'
#' Draws a circular genome plot: the cytoband ideogram of a genome build, a ring
#' of labels outside it, and one link per row joining a pair of genomic
#' intervals. It is the layout shared by every chart in the package that says
#' "these two places belong together".
#'
#' @section Why this one cannot return an object:
#' `circlize` is base graphics. There is no object to hand back: the functions
#' draw onto the current device as a side effect, so this drawing has to own the
#' device, which is the opposite of how [.plot_lollipop_build()] works. The
#' consequence for testing is stated rather than worked around: what can be
#' asserted here is that a file of plausible size appears, and everything that
#' decides *what* is drawn lives in the data-preparation function the caller uses
#' first, which is pure and is where the assertions go.
#'
#' @section The genome build is an argument, because it was a silent assumption:
#' The two functions this replaces both called `circlize::read.cytoband()` with
#' no species, which defaults to hg19. The session already records
#' `genome_build`, so a study annotated on hg38 was drawn against an hg19
#' ideogram: every gene placed at the wrong cytoband, with nothing anywhere
#' saying so. The build is now required and passed through.
#'
#' @section The gaps are computed, not counted by hand:
#' The previous code wrote `gap.after = c(rep(1, 23), 5, 5)`, which is correct
#' only for exactly 25 sectors. A different build, a filtered chromosome set or
#' a second pseudo-sector silently mismatches the sector count, and `circos.par`
#' recycles rather than complaining. The gaps are derived from the sectors
#' actually present.
#'
#' @param labels data.frame with `CHR`, `START`, `END` and `LABEL`. Drawn
#'   outside the ideogram.
#' @param links_from,links_to data.frames with `CHR`, `START`, `END` and the
#'   same number of rows: row `i` of one is joined to row `i` of the other.
#'   `NULL` for both draws the ideogram and the labels alone, which is a usable
#'   chart and not a degenerate one.
#' @param genome_build character. Passed to `circlize::read.cytoband()`, e.g.
#'   `"hg19"`, `"hg38"`.
#' @param plot_path character. File to write.
#' @param plot_format character. `"png"` or `"eps"`.
#' @param resolution_ppi numeric. Device resolution.
#' @param extra_sector character. Name of a pseudo-chromosome to append to the
#'   cytoband, for links whose far end is not a genomic position - a pathway, a
#'   chemical, a phenotype. `NULL` for none.
#' @param extra_sector_width numeric. Length of that pseudo-chromosome, in the
#'   same units as the genomic coordinates.
#' @param label_font_size,chr_font_size,label_connection_height,link_track_height
#'   Appearance, with the values the previous charts used.
#'
#' @return Invisibly the path written.
#'
#' @keywords internal
#' @noRd
.plot_circos_links_build <- function(labels, links_from, links_to,
                                     genome_build, plot_path, plot_format,
                                     resolution_ppi,
                                     extra_sector = NULL,
                                     extra_sector_width = 2e8,
                                     label_font_size = 0.3,
                                     chr_font_size = 0.3,
                                     label_connection_height = 0.03,
                                     link_track_height = 2) {

  if (!requireNamespace("circlize", quietly = TRUE))
    stop("circlize is needed to draw a circos plot and is not installed. ",
         "It is a suggested dependency: install.packages(\"circlize\").")

  required <- c("CHR", "START", "END", "LABEL")
  missing_columns <- setdiff(required, colnames(labels))
  if (length(missing_columns) > 0L)
    stop("labels is missing: ", paste(missing_columns, collapse = ", "))

  if (xor(is.null(links_from), is.null(links_to)))
    stop("links_from and links_to go together: pass both or neither.")
  if (!is.null(links_from) && nrow(links_from) != nrow(links_to))
    stop("links_from and links_to must have the same number of rows: row i of ",
         "one is joined to row i of the other. Got ", nrow(links_from),
         " and ", nrow(links_to), ".")

  cytoband <- circlize::read.cytoband(species = genome_build)$df
  if (!is.null(extra_sector))
    cytoband <- rbind(cytoband,
                      data.frame(V1 = extra_sector, V2 = 1,
                                 V3 = extra_sector_width, V4 = "", V5 = "",
                                 stringsAsFactors = FALSE))

  if (identical(plot_format, "eps"))
    grDevices::postscript(file = plot_path, width = 2480, height = 2480,
                          pointsize = 15)
  else
    # filename, not file: `file` partially matches `filename` and the package
    # turns partial-match warnings on, so the previous spelling warned on every
    # chart it drew.
    grDevices::png(filename = plot_path, width = 2480, height = 2480,
                   pointsize = 15, res = as.numeric(resolution_ppi))
  # The device is closed whatever happens below: a circos drawing that raises
  # half way through otherwise leaves the device open and every later plot in
  # the session lands in this file.
  on.exit({
    circlize::circos.clear()
    grDevices::dev.off()
  }, add = TRUE)

  sector_count <- length(unique(as.character(cytoband[[1]])))
  gaps <- c(rep(1, max(sector_count - 2L, 0L)), 5, 5)
  circlize::circos.par(gap.after = gaps, track.height = link_track_height)
  circlize::circos.genomicInitialize(cytoband, plotType = NULL)

  circlize::circos.genomicLabels(
    labels[, c("CHR", "START", "END", "LABEL")], labels.column = 4,
    side = "outside", cex = label_font_size,
    connection_height = label_connection_height)

  circlize::circos.track(
    track.index = circlize::get.current.track.index(),
    panel.fun = function(x, y) {
      circlize::circos.text(circlize::CELL_META$xcenter,
                            circlize::CELL_META$ylim[1],
                            circlize::CELL_META$sector.index,
                            niceFacing = TRUE, adj = c(0.5, 0),
                            cex = chr_font_size)
    },
    track.height = graphics::strheight("fj", cex = chr_font_size) * 2,
    bg.border = "grey", bg.col = "white", cell.padding = c(0, 0, 0, 0))

  circlize::circos.genomicIdeogram(cytoband)

  if (!is.null(links_from) && nrow(links_from) > 0L)
    circlize::circos.genomicLink(
      links_from[, c("CHR", "START", "END")],
      links_to[, c("CHR", "START", "END")],
      col = grDevices::rainbow(nrow(links_from)))

  invisible(plot_path)
}
