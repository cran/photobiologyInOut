#' Read libRadtran's uvspec output file from batch job.
#' 
#' Reads and parses the header and body of a text file output by a script used 
#' to run libRadtran's uvspec in a batch joib for a set of solar spectrum 
#' simulations. The header and time and date fields are converted into a 
#' datetime object.
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
#' @param multiplier numeric A multiplier for conversion into W m-2 nm-1, as
#'   the units of expression of the output from "uvspec" depend on the units
#'   in which the extraterrestrial solar spectrum data is expressed.
#' @param simplify logical Remove redundant columns from returned object.
#'   
#' @return a source_spct object, possibly containing several spectra in long
#'  form and a datetime column.
#'   
#' @references \url{http://www.r4photobiology.info} \url{http://www.libradtran.org}
#'
#' @export
#' 
read_uvspec_disort_vesa <- function(file, 
                                    date = NULL,
                                    geocode = NULL,
                                    label = NULL,
                                    tz = NULL,
                                    locale = readr::default_locale(),
                                    multiplier = 1e-6,
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
               ~s.e.irrad.dir * multiplier,
               ~s.e.irrad.diff * multiplier,
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
                      "' imported on ", lubridate::now(tzone = "UTC"), " UTC",
                      sep = "")
  photobiology::setWhenMeasured(z, datetimes)
  photobiology::setWhereMeasured(z, geocode)
  if (!is.na(label)) {
    photobiology::setWhatMeasured(z, paste("libRadtran spectral simulation", label))
  }
  z
}

#' Read libRadtran's uvspec output file.
#' 
#' Read and parse a text file output by libRadtran's uvspec routine for a solar 
#' spectrum simulation. The output of uvspec depends among other things on the
#' solver used. We define a family of functions, each function for a different
#' solver.
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
#' @param multiplier numeric A multiplier for conversion into W m-2 nm-1, as
#'   the units of expression of the output from "uvspec" depend on the units
#'   in which the extraterrestrial solar spectrum data is expressed.
#' @param qty character "uvspec" returns both irradiance and intensity with
#'   solver "disort".
#'   
#' @note Currently only "irradiance" is suported as qty argument as intensity is
#'   not supported by classes and methods in package 'photobiology'.
#'   
#' @return A source_spct object.
#'   
#' @references \url{http://www.r4photobiology.info} \url{http://www.libradtran.org}
#'
#' @note Tested only with libRadtran version 2.0
#' 
#' @export
#' 
read_uvspec_disort <- function(file, 
                               date = NULL,
                               geocode = NULL,
                               label = NULL,
                               tz = NULL,
                               locale = readr::default_locale(),
                               multiplier = 1e-3,
                               qty = "irradiance") {
  stopifnot(qty == "irradiance")
  # intensity may be supported in the future
  if (is.null(tz)) {
    tz <- locale$tz
  }
  if (is.null(label)) {
    label <- paste("File:", file)
  }
  z <- readr::read_table(file = file,
                         col_names = c("lambda", 
                                       "edir", "edn", "eup", 
                                       "uavgdir", "uavgdn", "uavgup"),
                         col_types = "ddddddd",
                         locale = locale)
  dots <- list(~lambda,
               ~edir * multiplier,
               ~edn * multiplier,
               ~(edir + edn) * multiplier)
  dots.names <- c("w.length", "s.e.irrad.dir", "s.e.irrad.diff", "s.e.irrad")
  z <- dplyr::transmute_(z, .dots = stats::setNames(dots, dots.names))
  
  photobiology::setSourceSpct(z, time.unit = "second", multiple.wl = 1)
  comment(z) <- paste("libRadtran file '", file,
                      "' imported on ", lubridate::now(tzone = "UTC"), " UTC",
                      sep = "")
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  if (!is.na(label)) {
    photobiology::setWhatMeasured(z, paste("libRadtran spectral simulation", label))
  }
  z
}

# Allows function to work with either dplyr 0.4 (which ignores value of
# starts_with), and 0.5 which exports it as a proper function
starts_with <- function(...) dplyr::starts_with(...)

