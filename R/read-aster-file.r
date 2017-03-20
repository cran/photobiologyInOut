#' Read File downloaded from ASTER data base.
#' 
#' Reads and parses the header of a test file as available through the
#' ASTER reflectance database. The Name field is retrieved and copied to
#' attribute "what.measured". The header of the file is preserved as a comment.
#' 
#' @param file character string
#' @param date a \code{POSIXct} object, but if \code{NULL} the date stored in
#'   file is used, and if \code{NA} no date variable is added
#' @param geocode A data frame with columns \code{lon} and \code{lat}.
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
#' @references \url{http://www.r4photobiology.info} \url{https://speclib.jpl.nasa.gov/}
#' Baldridge, A.; Hook, S.; Grove, C. & Rivera, G. (2009) The ASTER spectral 
#' library version 2.0. Remote Sensing of Environment. 113, 711-715
#' 
#' @note The header in these files has very little information, so the user
#' needs to supply the number of pixels in the array as well as the date-time.
#' The file contains a date in milliseconds but as the Raspberry Pi board
#' contains no real-time clock, it seems to default to number of milliseconds
#' since the Pi was switched on.
#' 
read_ASTER_txt <- function(file,
                           date = NULL,
                           geocode = NULL,
                           label = NULL,
                           tz = NULL,
                           locale = readr::default_locale(),
                           npixels = 2048) {
  
  file_header <- scan(file = file, nlines = 26, 
                      skip = 0, what = "character", sep = "\n")
  
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
  label <- paste(sub("Name: ", "", file_header[[1]]), label, sep = "\n") 
  photobiology::setWhatMeasured(z, label)
  z
}

