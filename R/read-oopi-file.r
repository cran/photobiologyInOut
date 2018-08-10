#' Read File Saved by Ocean Optics' Raspberry Pi software.
#' 
#' Reads and parses the header of a raw data file as output by the server
#' running on a Raspberry Pi board to extract the whole header remark field. The
#' time field is retrieved and decoded.
#' 
#' @param file character string
#' @param date a \code{POSIXct} object, but if \code{NULL} the date stored in
#'   file is used, and if \code{NA} no date variable is added
#' @param geocode A data frame with columns \code{lon} and \code{lat}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone is not saved to the file.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#' @param npixels integer Number of pixels in spectral data.
#'   
#' @return A raw_spct object.
#' @export
#' @references \url{http://www.r4photobiology.info} \url{http://oceanoptics.com/} \url{https://www.raspberrypi.org/}
#' 
#' @note The header in these files has very little information, so the user
#' needs to supply the number of pixels in the array as well as the date-time.
#' The file contains a date in milliseconds but as the Raspberry Pi board
#' contains no real-time clock, it seems to default to number of milliseconds
#' since the Pi was switched on.
#' 
read_oo_pidata <- function(file,
                           date = NULL,
                           geocode = NULL,
                           label = NULL,
                           tz = NULL,
                           locale = readr::default_locale(),
                           npixels = 2048) {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
  label <- paste("File:", basename(file), label)
  
  file_header <- scan(file = file, nlines = 4, 
                      skip = 0, what = "character", sep = "\n")
  
  if (is.null(date)) {
    date <- sub("Saved at time: ", "", file_header[1], fixed = TRUE)
    # needs further decoding if possible
  }
  
  z <- readr::read_tsv(
    file = file,
    col_names = c("w.length", "counts"),
    skip = 5,
    n_max = npixels,
    col_types = readr::cols(),
    locale = locale
  )
  
  z <- photobiology::as.raw_spct(z)

  comment(z) <-
    paste(paste("Ocean Optics Raspeberry Pi raw counts file '", basename(file), "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")

  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  attr(z, "file.header", file_header)
  z
}

