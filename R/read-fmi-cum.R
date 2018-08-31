#' Read daily cummulated solar spectrum data file(s).
#'
#' Read one or more cumulated daily spectral irradiance file as output by Anders
#' Lindors' model based on libRadTrans. Dates are read from the file header and
#' parsed with the function suplied as \code{date.f}.
#'
#' @param file Either a path to a file, a connection, or literal data (either a
#'   single string or a raw vector).
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date is extracted from the
#'   file header.
#' @param geocode A data frame with columns \code{lon} and \code{lat} used to
#'   set attribute \code{"where.measured"}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone used for interpreting times saved in the file
#'   header.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#' @param .skip Number of lines to skip before reading data---i.e. the number of
#'   rows in the header.
#' @param .n_max Maximum number of records to read.
#' @param .date.f A function for extracting a date-time from the file header
#'   passed as charecter sring to its first argument and which returns a
#'   \code{POSIXct} object.
#'
#' @return \code{read_fmi_cum()} returns a \code{source_spct} object with
#'   \code{time.unit} attribute set to \code{"day"} and \code{when.measured}
#'   attribute set to the date-time extracted from the header at the top of
#'   the read file.
#'
#' @note See \code{\link[readr]{read_table}} for details of acceptable values
#'   for \code{file}.
#'
#' @examples
#' 
#' file.name <- system.file("extdata", "2014-08-21_cum.hel", 
#'                          package = "photobiologyInOut", mustWork = TRUE)
#' fmi.spct <- read_fmi_cum(file = file.name)
#'   
#' @export
#' 
read_fmi_cum <- function(file,
                         date = NULL,
                         geocode = NULL,
                         label = NULL,
                         tz = "UTC",
                         locale = readr::default_locale(),
                         .skip = 3,
                         .n_max = -1,
                         .date.f = lubridate::ymd) {
  file_header <- readr::read_lines(file = file, 
                                   n_max = .skip, 
                                   progress = FALSE)
  
  if (!grepl("(J/m2/nm)", file_header[.skip], fixed = TRUE)) {
    warning("Skipping file with unrecognized format!")
    return(source_spct())
  }
  
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
#  file_date <- file.mtime(file, tz = tz)[1]
  
  if (is.null(date)) {
    date.char <- scan(text = file_header[2], nlines = 1, what = "z")[2]
    date.time <- .date.f(date.char, tz = tz) # tz needed as otherwise a "Date" is returned
  }
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  z <- readr::read_table(
    file = file,
    col_names = c("w.length", "s.e.irrad"),
    col_types = readr::cols(w.length = readr::col_double(),
                            s.e.irrad = readr::col_double()),
    skip = .skip,
    n_max = .n_max,
    locale = locale,
    progress = FALSE
  )
  # convert wavelength in Ångstrom to nm
  if (min(z$w.length) > 1000) {
    dots <- list(~w.length / 10)
    z <- dplyr::mutate_(z, .dots = stats::setNames(dots, "w.length"))
  }
  z <- photobiology::as.source_spct(z, time.unit = "day")
  photobiology::setWhenMeasured(z, date.time)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  attr(z, "file.header") <- file_header
  z
}

#' @rdname read_fmi_cum
#'
#' @param files list or vector of paths each one with the same requirements as
#'    described for argument \code{file}.
#'
#' @return \code{read_m_fmi_cum} returns a \code{source_mspct} containing one
#' \code{source_spct} object for each one of the multiple files read.
#'
#' @export
#'
read_m_fmi_cum <- function(files,
                           date = NULL,
                           geocode = NULL,
                           label = NULL,
                           tz = "UTC",
                           .skip = 3,
                           .n_max = -1,
                           .date.f = lubridate::ymd) {
  # extract ISO formatted dates from file names

  list.of.spectra <- list()

  for (f in files) {
    data.name <- sub("^.*/", "", f) # remove path!
    data.name <- gsub("-", "_", data.name, fixed = TRUE) # sanitize
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

#' Read multiple solar spectra from a data file.
#'
#' Read spectral irradiance file as output by Anders Lindors' model based on
#' libRadTrans for hourly simulation.
#'
#' @param file Either a path to a file, a connection, or literal data (either a
#'   single string or a raw vector).
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
#' @param .skip Number of lines to skip before reading data.
#' @param .n_max Maximum number of records to read.
#'
#' @return \code{read_fmi2mspct()} returns a \code{source_mspct} object 
#'   containing \code{source_spct} objects as members, \code{time.unit} 
#'   attribute set to \code{"second"} and \code{when.measured}
#'   attribute set to the date-time values extracted from the file body.
#'
#' @note See \code{\link[readr]{read_table}} for details of acceptable values
#'  for \code{file}.
#'
#' @export
#'  
read_fmi2mspct <- function(file,
                           geocode = NULL,
                           label = NULL,
                           tz = NULL,
                           locale = readr::default_locale(),
                           .skip = 3,
                           .n_max = -1) {
  
  get_one_spct <- function(x, i, j, k) {
    header <- x[i]
    date.char <- stringr::word(header, start = 2L, end = 3L)
    # lubridate::ymd_hms() fails to parse single digit hours
    # tz = tz is needed to get POSIXct instead Date which ends as POSIXlt
    when.measured <- 
      lubridate::ymd(stringr::word(header, start = 2L, end = 2L), tz = tz) +
      lubridate::hms(stringr::word(header, start = 3L, end = 3L))
    zenith.angle <- as.numeric(stringr::word(header, start = 5L, end = 5L))
    
    z <-
      readr::read_table2(
        file = file,
        col_names = c("minutes", "w.length", "s.e.irrad"),
        col_types = readr::cols(minutes = readr::col_number(),
                                w.length = readr::col_double(),
                                s.e.irrad = readr::col_double()),
        skip = i,
        n_max = k - j,
        locale = locale,
        progress = FALSE
      )[ , -1]
    # convert wavelength in Ångstrom to nm
    # convert irradiance to W m-2 nm-1
    if (min(z$w.length) > 1000) {
      dots <- list(~w.length / 10, ~s.e.irrad / 1000)
      z <- dplyr::mutate_(z, .dots = stats::setNames(dots, c("w.length", "s.e.irrad")))
    }
    
    z <- photobiology::as.source_spct(z, time.unit = "second")
    photobiology::setWhenMeasured(z, when.measured)
    photobiology::setWhereMeasured(z, geocode)
    photobiology::setWhatMeasured(z, label)
    attr(z, "file.header") <- header
    z<- list(z)
    names(z) <- date.char
    z
  }
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  whole_file <- readr::read_lines(file, progress = FALSE)
  header.start.idxs <- grep("# ", whole_file, fixed = TRUE)
  data.start.idxs <- header.start.idxs + 1L
  data.end.idxs <- grep("end", whole_file, fixed = TRUE)
  if (length(header.start.idxs) == 0 ||
      length(header.start.idxs) != length(data.start.idxs) ||
      length(header.start.idxs) != length(data.end.idxs)) {
    warning("Unmatched delimiters or no delimiters found!")
    return(source_mspct())
  }
  
  zz <- list()
  for (idx in seq(along.with = header.start.idxs)) {
    zz <- c(zz,
            get_one_spct(x = whole_file,
                         i = header.start.idxs[idx],
                         j = data.start.idxs[idx],
                         k = data.end.idxs[idx]
            ))
  }
  zz <- photobiology::as.source_mspct(zz)
  comment(zz) <- label
  zz
}
