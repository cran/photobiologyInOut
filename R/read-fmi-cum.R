#' Read daily cummulated solar spectrum data file(s).
#'
#' Read one or more cumulated daily spectral irradiance file as output by Anders
#' Lindors' model based on libRadTrans. The file naming conventions needed
#' are fairly strict, and file name should contain the date in a format
#' suitable for decoding by the function suplied as \code{date.f}.
#'
#' @param file Either a path to a file, a connection, or literal data (either a
#'   single string or a raw vector).
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
#' @param .skip Number of lines to skip before reading data.
#' @param .n_max Maximum number of records to read.
#' @param .date.f A function for extracting a date-time from the file name
#'   passed as charecter sring to its first argument and which returns a
#'   \code{POSIXct} object.
#'
#' @return \code{read_fmi_cum()} returns a \code{source_spct} object with
#'   \code{time.unit} attribute set to \code{"day"} and \code{when.measured}
#'   attribute set to the date-time extracted from the file name.
#'
#' @note See \code{\link[readr]{read_table}} for details of acceptable values
#'  for \code{file}.
#'
#' @export
#'  
read_fmi_cum <- function(file,
                         date = NULL,
                         geocode = NULL,
                         label = NULL,
                         tz = NULL,
                         locale = readr::default_locale(),
                         .skip = 3,
                         .n_max = -1,
                         .date.f = lubridate::ymd) {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
  label <- paste("File:", basename(file), label)

  z <- readr::read_table(
    file = file,
    col_names = c("w.length", "s.e.irrad"),
    col_types = "dd",
    skip = .skip,
    n_max = .n_max,
    locale = locale
  )
  # convert wavelength in Ångstrom to nm
  if (min(z$w.length) > 1000) {
    dots <- list(~w.length / 10)
    z <- dplyr::mutate_(z, .dots = stats::setNames(dots, "w.length"))
  }
  z <- photobiology::as.source_spct(z, time.unit = "day")
  if (is.null(date)) {
    date <- file.mtime(file, tz = tz)[1]
  }
    photobiology::setWhenMeasured(z, date)
    photobiology::setWhereMeasured(z, geocode)
    photobiology::setWhatMeasured(z, label)
  z
}

#' @rdname read_fmi_cum
#'
#' @param files list or vector of paths each one with the same requirements as
#'    described for argument \code{file}.
#'
#' @return \code{read_m_fmi_cum} returns a collection of \code{source_mspct}.
#'
#' @export
#'
read_m_fmi_cum <- function(files,
                           date = NULL,
                           geocode = NULL,
                           label = NULL,
                           tz = NULL,
                           .skip = 3,
                           .n_max = -1,
                           .date.f = lubridate::ymd) {
  # extract ISO formatted dates from file names

  list.of.spectra <- list()

  for (f in files) {
    data.name <- sub("^.*/", "", f) # remove path!
    data.name <- gsub("-", "_", data.name, fixed = TRUE) # sanitize
#    data.name <- paste(data.name, ".spct", sep = "")
    list.of.spectra[[data.name]] <-
      read_fmi_cum(
        file = f,
        date = date,
        geocode = geocode,
        label = label,
        tz = tz,
        .skip = .skip,
        .n_max = .n_max,
        .date.f = .date.f
      )
  }
  photobiology::source_mspct(list.of.spectra)
}
