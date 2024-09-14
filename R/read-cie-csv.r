#' Read '.CSV' files from CIE
#'
#' Reads a CSV spectral data file and its companion JSON file with metadata as
#' published by \emph{International Commission on Illumination} (CIE) and then 
#' imports wavelengths and spectral values into one the classes for spectral
#' data defined in package 'photobiology.
#' 
#' @details
#' The CSV file contains only numbers encoded as character strings, and the
#' JSON file contains extensive metadata. The type of spectral data is encoded
#' as part of the file name. If the original file name of the CSV file is passed
#' as argument to parameter \code{file}, the function can retrieve all data and
#' metadata, enough to return an R object of the correct class. The JSON file
#' must be located in the same folder.
#'  
#' @param file.name character string
#' @param label character string, but if \code{NULL} metadata read from the
#'   JSON file is used, and if \code{NA} the \code{"what.measured"} attribute 
#'   is not set, and if a character string is passed, it is used to set 
#'   the \code{"what.measured"} attribute.
#' @param simplify logical If \code{TRUE} and the read file contained a single
#'   spectrum, extract the spectral object from the collection.
#'
#' @return Depending on the contents of the file, a \code{source_spct} object, a
#'   \code{response_spct} object, or a \code{chroma_spct} object, containing
#'   both data and metadata.
#'  
#' @export
#' 
#' @references \url{https://cie.co.at/data-tables}
#' 
#' @keywords misc
#' 
#' @examples
#' 
#' file.name <- 
#'   system.file("extdata", "CIE_illum_C.csv", 
#'               package = "photobiologyInOut", mustWork = TRUE)
#' 
#' CIE_illum_C.spct <- read_CIE_csv(file.name)
#' CIE_illum_C.spct
#'    
#' file.name <- 
#'   system.file("extdata", "CIE_sle_photopic.csv", 
#'               package = "photobiologyInOut", mustWork = TRUE)
#' 
#' CIE_sle_photopic.spct <- read_CIE_csv(file.name)
#' CIE_sle_photopic.spct
#'    
read_CIE_csv <- function(file.name,
                         label = NULL,
                         simplify = FALSE) {
  
  metadata.file.name <- paste(file.name, "metadata.json", sep = "_")
  base.name <- basename(file.name)
  
  if (grepl("^CIE_illum|^CIE_std_illum", base.name)) {
    member.class <- "source_spct"
    spct.data.var <- "s.e.irrad"
    f <- "identity"
    target <- NA_real_
  } else if (grepl("^CIE_sle|^CIE_cfb_sle", base.name)) {
    member.class <- "response_spct"
    spct.data.var <- "s.e.response"
    f <- "identity"
    target <- NA_real_
  }
  
  cie_metadata.ls <- jsonlite::fromJSON(metadata.file.name)
  
  col.names <- cie_metadata.ls$datatableInfo$columnHeaders$title
  
  cie.tb <- 
    utils::read.csv(file.name, 
                    header = FALSE,
                    col.names = col.names)
  
  z <- 
    photobiology::split2mspct(cie.tb, 
                              member.class = member.class,
                              spct.data.var = spct.data.var,
                              w.length.var = "lambda")
  
  scaling <- list(multiplier = 1, 
                  f = f, 
                  range = wl_range(z[[1]]), 
                  target = target, 
                  cols = spct.data.var)
  
  photobiology::setScaled(z, scaled = scaling)
  
  label.text <- paste("file:", base.name, label,
                 "\nCIE:", cie_metadata.ls$titles$title)
  
  comment.text <- 
    with(cie_metadata.ls,
         paste(identifier$identifierType, ": ", identifier$identifier,
               "; ", alternateIdentifiers$alternateIdentifierType, ": ", 
               alternateIdentifiers$alternateIdentifier,
               "; units: ", 
               paste(datatableInfo$columnHeaders$unit, collapse = ", "),
               "\n", publisher,
               ". ", rightsList$rightsIdentifier,
               sep = "")       
    )
  
  z <- photobiology::msmsply(z, 
                             photobiology::setWhatMeasured, 
                             what.measured = label.text)
  z <- photobiology::msmsply(z, 
                             photobiology::setHowMeasured, 
                             how.measured = 
                               cie_metadata.ls$relatedItems$titles[[1]])
  z <- photobiology::msmsply(z, 
                             `comment<-`,
                             value = comment.text)
  
  if (simplify && length(z) == 1) {
    z[[1]]
  } else {
    z
  }
  
 }

# CIE_D55_2018.mspct <- read_CIE_csv("./inst-not/CIE-JSON/CIE_illum_D55.csv")
# CIE_D75_2018.mspct <- read_CIE_csv("./inst-not/CIE-JSON/CIE_illum_D75.csv")
# CIE_LEDs_2018.mspct <- read_CIE_csv("./inst-not/CIE-JSON/CIE_illum_LEDs.csv")
# 
# autoplot(CIE_D55_2018.mspct) + geom_point()
# autoplot(CIE_D75_2018.mspct) + geom_point()
# autoplot(CIE_LEDs_2018.mspct)
