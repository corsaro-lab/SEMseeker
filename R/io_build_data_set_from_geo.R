#' io_build_data_set_from_geo
#'
#' @param GEOgse geo accession dataset identification
#' @param result_folder where sample sheet and files will be saved
#' @param downloadFiles 0 means download all files from Gene Expression Omnibus (GEO),
#'  different than zero means how many to download
#' @param ... additional arguments passed to internal helpers
#'
#' @return samplesheet, and sample's file saved and samplesheet csv
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' io_build_data_set_from_geo(
#'   GEOgse        = "GSE55763",
#'   downloadFiles = 1,
#'   result_folder = tempdir()
#' )
#' }
#' @export
io_build_data_set_from_geo <-  function(GEOgse, downloadFiles = 0,result_folder, ...) {


  ssEnv <- core_init_env( result_folder= result_folder, ...)
  # check if GEOQuery package is installed
  if(!requireNamespace("GEOquery", quietly = TRUE))
  {
    message("GEOquery package is not installed. Please install GEOquery package to use this function")
    return()
  }

    if(file.exists(GEOgse))
    {
      gse <- GEOquery::getGEO(filename = GEOgse, GSEMatrix = TRUE)
      phenoDataName <- basename(GEOgse)
      samplesheet <- gse@phenoData@data
    }
    else
    {
      gse <- GEOquery::getGEO(GEO = GEOgse, GSEMatrix = TRUE)
      phenoDataName <-  paste(GEOgse, "_series_matrix.txt.gz", sep = "")
      samplesheet <- gse[[phenoDataName]]@phenoData@data
    }

    samplesheet$Sample_ID <- samplesheet$geo_accession
    samplesheet$Sample_Group <- ""
    samplesheet$Sentrix_ID <- ""
    samplesheet$Sentrix_Position <- ""

    utils::write.table(
      samplesheet,
      paste(result_folder, "/", "final_samplesheet.csv", sep = ""),
      row.names = FALSE,
      sep="\t",
      quote = FALSE
    )


    io_dir_check_and_create(result_folder,"/")
    if(downloadFiles==0)
      downloadFiles <- nrow(samplesheet)

    downloadedFiles <- 0
    for(sample in seq_len(nrow(samplesheet)))
    {
      sample <- 1
      core_log_event(sample)
      fileName <- samplesheet [sample,"supplementary_file"]
      if(fileName=="" || fileName=="NONE")
        next
      tempData <-  gsub("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/", "", fileName)
      tempData <-  gsub("suppl/", "", tempData)
      tempData <-  gsub("_Grn.idat.gz", " ", tempData)
      tempData <-  gsub("_Red.idat.gz", " ", tempData)
      tempData <-  noquote(strsplit(tempData, "/"))
      tempData <-  strsplit(tempData[[1]], " ")
      tempData <-  tempData[length(tempData)] [[1]]
      tempData <-  unlist(strsplit(tempData, "_"))
      Sentrix_ID <- tempData[length(tempData)-1]
      Sentrix_Position <- tempData[length(tempData)]
      samplesheet [sample,"Sentrix_ID"] <- Sentrix_ID
      samplesheet [sample,"Sentrix_Position"] <- Sentrix_Position

      if (- downloadedFiles + downloadFiles > 0 ) {

        fileName <-  gsub("_Grn.idat.gz", "_COLOR.idat.gz", fileName)
        fileName <-  gsub("_Red.idat.gz", "_COLOR.idat.gz", fileName)

        for (i in 1:2)
        {

          color <- c("_Grn.idat.gz","_Red.idat.gz")[i]
          # color <- "_Grn.idat.gz"
          localFileName <- paste(Sentrix_ID,"_", Sentrix_Position, color, sep="")
          localFileName <- paste(result_folder, "/", localFileName, sep = "")
          localFileNameUnzipped <- gsub(".gz","",localFileName)
          if(!file.exists(localFileName) &&  !file.exists(localFileNameUnzipped))
          {
            fileName <-  gsub("_COLOR.idat.gz",color, fileName)
            utils::download.file(
              url = fileName,
              destfile =  localFileName,
              # method =  "internal",
              quiet = FALSE,
              mode = "wb",
              cacheOK = TRUE,
              extra = getOption("download.file.extra"),
              headers = NULL
            )
            # GEOquery dropped the gunzip re-export in recent versions
            # (R.utils::gunzip is the original delegate). Calling the source
            # directly avoids the runtime "not an exported object" error.
            R.utils::gunzip(localFileName, overwrite = TRUE)
          }
        }
        downloadedFiles <- downloadedFiles + 1
      }
    }

    #remove commas from  values
    samplesheet <- apply(samplesheet, 2, function(x){gsub(",","_",x)})

    utils::write.table(
      samplesheet,
      paste(result_folder, "/", "final_samplesheet.csv", sep = ""),
      row.names = FALSE,
      sep=",",
      quote = FALSE
    )

    return(as.data.frame(samplesheet))
}
