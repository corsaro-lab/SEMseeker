plot_area_plot_circlize <- function(areas, circle_plot_path,result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
{
  ssEnv <- core_init_env(result_folder = result_folder, maxResources = maxResources, parallel_strategy = parallel_strategy, start_fresh = FALSE, ...)

  cytoband <- circlize::read.cytoband()$df
  cytoband <- rbind(cytoband,
    data.frame(V1 = "", V2 = 1,  V3 = 2e8, V4 = "", V5 = "")
  )

  # Variables to control the font size, link length, chr font size, and link start height
  label_font_size <- 0.3
  chr_font_size <- 0.3

  label_connection_height <- 0.03

  link_start_height <- 2
  link_length <- 2

  # Save the plot as a PNG file
  if(ssEnv$plot_format == "png")
    grDevices::png(file =  circle_plot_path, width = 2480,height = 2480, pointsize  =  15, res = 400)
  if(ssEnv$plot_format == "eps")
    grDevices::postscript(file =  circle_plot_path, width = 2480,height = 2480, pointsize  =  15, res = 400)
  # grDevices::png(circle_plot_path, width = 2048, height = 2048, res = 400)

  # Set circos parameters to adjust the gap and overall appearance
  circlize::circos.par(gap.after = c(rep(1, 23), 5, 5), track.height = link_start_height)

  # Initialize the circos plot without plotting anything yet
  circlize::circos.genomicInitialize(cytoband, plotType = NULL)

  # Prepare label data frame
  areas <- areas[, c("CHR", "START", "END", "AREA")]

  # Add genomic labels with adjusted size and position
  circlize::circos.genomicLabels(areas, labels.column = 4, side = "outside", cex = label_font_size, connection_height = label_connection_height)

  # Add chromosome names with adjusted size and positioning
  circlize::circos.track(track.index = circlize::get.current.track.index(), panel.fun = function(x, y) {
    circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1], circlize::CELL_META$sector.index,
      niceFacing = TRUE, adj = c(0.5, 0), cex = chr_font_size)
  }, track.height = strheight("fj", cex = chr_font_size) * 2, bg.border = "grey", bg.col= "white", cell.padding = c(0, 0, 0, 0))

  # Add the ideogram track
  circlize::circos.genomicIdeogram(cytoband)

  # Add genomic links with unique colors
  my_col <- grDevicesrainbow(nrow(results))
  circlize::circos.genomicLink(results[, 1:3], results[, 6:8], col = my_col)

  # Clear the circos plot after drawing
  circlize::circos.clear()

  # Close the PNG device
  grDevices::dev.off()
}
