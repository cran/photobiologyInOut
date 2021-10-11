#' Read File downloaded from ASTER data base.
#' 
#' Reads and parses the header of a test file as available through the
#' ASTER reflectance database. The Name field is retrieved and copied to
#' attribute "what.measured". The header of the file is preserved as a comment.
#' 
#' @param file character string
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date is extracted from the
#'   file header.
#' @param geocode A data frame with columns \code{lon} and \code{lat} used to
#'   set attribute \code{"where.measured"}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Ignored.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#' @param npixels integer Number of pixels in spectral data.
#'   
#' @return A raw_spct object.
#' @export
#' @references 
#' 
#' \url{https://speclib.jpl.nasa.gov}
#' 
#' Baldridge, A.; Hook, S.; Grove, C. & Rivera, G. (2009) The ASTER spectral 
#' library version 2.0. Remote Sensing of Environment. 113, 711-715
#' 
#' @note The header in these files has metadata information, but mostly on the
#'   origin of the data. For a date and/or geocode are to be added to the return
#'   object it must be supplied by the user.  as well as the date-time. Some
#'   metadata is extracted and added as attributes, while the whole header is
#'   copied to the \code{comment} attribute.
#' 
#' @examples
#' 
#'  file.name <- 
#'    system.file("extdata", "drygrass-spectrum.txt", 
#'                package = "photobiologyInOut", mustWork = TRUE)
#'                 
#'  fred.spct <- read_ASTER_txt(file = file.name, npixels = Inf)
#'  
#'  fred.spct
#'  getWhatMeasured(fred.spct)
#'  cat(comment(fred.spct))
#'  
read_ASTER_txt <- function(file,
                           date = NULL,
                           geocode = NULL,
                           label = NULL,
                           tz = NULL,
                           locale = readr::default_locale(),
                           npixels = Inf) {
  
  file_header <- scan(file = file, nlines = 26, 
                      skip = 0, what = "character", sep = "\n", quiet = TRUE)
  
  label.file <- paste("File: ", basename(file), sep = "")
  label.file <- paste(file_header[[1]], label.file, sep = "\n") 
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
 
  z <- utils::read.table(
    file = file,
    col.names = c("w.length", "Rpc"),
    skip = 26
  )
  
  z[["w.length"]] <- z[["w.length"]] * 1e3 # um -> nm
  
  z <- photobiology::as.reflector_spct(z, Rfr.type = "total")

  comment(z) <-
    paste(paste("ASTER database file '", basename(file), "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")

  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  attr(z, "file.header") <- file_header
  z
}

