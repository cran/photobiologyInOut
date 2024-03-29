#' Read '.DTA' File Saved by Macam's Software.
#'
#' Reads and parses the header of a processed data file as output by the PC
#' program to extract the time and date fields and a user label if present,
#' and then imports wavelengths and spectral energy irradiance values. 
#'
#' @param file character string
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date is extracted from the
#'   file header.
#' @param geocode A data frame with columns \code{lon} and \code{lat} used to
#'   set attribute \code{"where.measured"}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone used for interpreting times saved in the
#'   file header.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#'
#' @return A \code{source_spct} object.
#' @export
#' @references \url{https://www.irradian.co.uk/}
#' 
#' @examples
#' 
#'  file.name <- 
#'    system.file("extdata", "spectrum.DTA", 
#'                package = "photobiologyInOut", mustWork = TRUE)
#'                 
#'  macam.spct <- read_macam_dta(file = file.name)
#'  
#'  macam.spct
#'  getWhenMeasured(macam.spct)
#'  getWhatMeasured(macam.spct)
#'  cat(comment(macam.spct))
#' 
read_macam_dta <- function(file,
                           date = NULL,
                           geocode = NULL,
                           label = NULL,
                           tz = NULL,
                           locale = readr::default_locale()) {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  file_header <- scan(file = file, nlines = 3, skip = 0,
                      what = "character", quiet = TRUE)
  if (is.null(date)) {
    date <- lubridate::dmy(sub(pattern = "@", replacement = "",
                               x = file_header[1], fixed = TRUE),
                           tz = tz)
    time <- lubridate::hms(sub(pattern = '@', replacement = "",
                               x = file_header[2], fixed = TRUE))
    date <- date + time
  }
  z <- scan(file = file,
                   what = list(w.length = double(), s.e.irrad = double()),
                   skip = 3, quiet = TRUE)

  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.source_spct(z, time.unit = "second")
  options(old.opts)

  comment(z) <-
    paste(paste("MACAM file '", basename(file), "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  how <- "Measured with a double monochromator scanning spectroradiometer."
  photobiology::setHowMeasured(z, how)
  attr(z, "file.header") <- file_header
  z
}
