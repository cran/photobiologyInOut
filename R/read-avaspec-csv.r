#' Read \code{.csv} File Saved by Aavnates' Software for AvaSpec.
#' 
#' Reads and parses the header of a processed data file as output by the
#' program Avaspec and then imports wavelength and spectral irradiance
#' values. The file header has little useful metadata information.
#' 
#' @param file character string
#' @param date a \code{POSIXct} object, but if \code{NULL} the date stored in
#'   file is used, and if \code{NA} no date variable is added
#' @param geocode A data frame with columns \code{lon} and \code{lat}.
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
#' @return A source_spct object.
#' @export
#' @author Pedro J. Aphalo
#' @references \url{http://www.r4photobiology.info}
#' @keywords misc
#' 
read_avaspec_csv <- function(file,
                             date = NULL,
                             geocode = NULL,
                             label = NULL,
                             tz = NULL,
                             locale = readr::default_locale()) {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  if (is.null(label)) {
    label <- paste("File:", file)
  }
  file_header <- scan(file = file, nlines = 6, skip = 0, what = "character")
  # watt / cm ?
  if (length(grep("Watt/cm", file_header[2], fixed = TRUE))) {
    mult <- 10e-4
  } else {
    mult <- NA
  }
  
  z <- readr::read_csv(file = file,
                       col_names = c("w.length", "s.e.irrad"),
                       skip = 6,
                       locale = locale)
  dots <- list(~s.e.irrad * mult)
  z <- dplyr::mutate_(z, .dots = stats::setNames(dots, "s.e.irrad"))
  
  old.opts <- options("photobiology.strict.range" = NA)
  z <- photobiology::as.source_spct(z, time.unit = "second")
  options(old.opts)
  
  comment(z) <-
    paste(paste("Avantes AvaSpec irradiance file '", file, "' imported on ", 
                lubridate::now(tz = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  z
}
