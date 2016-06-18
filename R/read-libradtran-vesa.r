#' Read libRadtran output file.
#' 
#' Reads and parses the header of a text file output by libRadtran for a solar 
#' spectrum simulation to extract the header and spectral data. The time and 
#' date fields are converted into a datetime object.
#' 
#' @param file character string
#' @param date a \code{POSIXct} object, but if \code{NULL} the date stored in
#'   file is used, and if \code{NA} no date variable is added
#' @param geocode A data frame with columns \code{lon} and \code{lat}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone is by default read from the file.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#' @param simplify logical Remove redundant columns from returned object.
#'   
#'   
#' @return a source_spct object, possibly containing several spectra in long
#'  form and a datetime column.
#'   
#' @references \url{http://www.r4photobiology.info}
#'
#' @note Tested only libRadtran version 2.0
#' 
#' @export
#' 
read_libradtran_vesa <- function(file, 
                                 date = NULL,
                                 geocode = NULL,
                                 label = NULL,
                                 tz = NULL,
                                 locale = readr::default_locale(),
                                 simplify = TRUE) {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  if (is.null(label)) {
    label <- paste("File:", file)
  }
  z <- readr::read_table(file = file,
                         col_names = c("w.length", "day", "time", 
                                       "s.e.irrad.dir", "s.e.irrad.diff"),
                         col_types = "dccdd",
                         locale = locale)
  dots <- list(~lubridate::ymd_hms(paste(day, time), tz = tz),
               ~s.e.irrad.dir * 1e-3,
               ~s.e.irrad.diff * 1e-3,
               ~s.e.irrad.dir + s.e.irrad.diff)
  dots.names <- c("datetime", "s.e.irrad.dir", "s.e.irrad.diff", "s.e.irrad")
  z <- dplyr::mutate_(z, .dots = stats::setNames(dots, dots.names))

  datetimes <- unique(z[["datetime"]])
  num.spectra <- length(datetimes)
  if (simplify && num.spectra == 1) {
    z <- dplyr::select_(z, "w.length", 
                        lazyeval::interp(~starts_with(x), x = "s.e.irrad"))
  } else if (simplify && num.spectra > 1) {
    z <- dplyr::select_(z, "w.length", "datetime", 
                        lazyeval::interp(~starts_with(x), x = "s.e.irrad"))
  }
  photobiology::setSourceSpct(z, time.unit = "second", multiple.wl = num.spectra)
  comment(z) <- paste("libRadtran file '", file,
                      "' imported on ", lubridate::now(tz = "UTC"), " UTC",
                      sep = "")
  photobiology::setWhenMeasured(z, datetimes)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, paste("libRadtran spectral simulation", label))
  z
}

# Allows function to work with either dplyr 0.4 (which ignores value of
# starts_with), and 0.5 which exports it as a proper function
if (packageVersion("dplyr") > '0.4.3') {
  starts_with <- function(...) dplyr::starts_with(...)
}

