#' Read '.PRN' File(s) Saved by LI-COR's PC1800 Program.
#' 
#' Reads and parses the header of a processed data file as output by the PC1800
#' program to extract the whole header remark field and also check whether data
#' is in photon or energy based units. The time field is ignored as it does not
#' contain year information. This instrument is no longer being manufactured.
#' 
#' @param file Path to file as a character string.
#' @param date a \code{POSIXct} object, but if \code{NULL} the date stored in
#'   file is used, and if \code{NA} no date variable is added.
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
#' @return \code{read_licor_prn()} returns a \code{source_spct} object with
#'   \code{time.unit} attribute set to \code{"second"} and \code{when.measured}
#'   attribute set to the date-time extracted from the file name, or supplied.
#' @export
#' @references \url{http://www.r4photobiology.info} \url{http://www.licor.com}
#' 
#' @keywords misc
#'   
#' @note The LI-1800 spectroradiometer does not store the year as part of the
#' data, only month, day, and time of day. Because of this, in the current
#' version, if \code{NULL} is the argument to date, year is set to 0000.
#' 
read_licor_prn <- function(file,
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
  file_header <- scan(
    file = file,
    nlines = 7,
    skip = 0,
    what = "character",
    sep = "\n"
  )
  
  if (is.null(date)) {
    line05 <- sub("Date:", "", file_header[5])
    date <- lubridate::parse_date_time(line05, "mdHM", tz = tz)
  }
  
  if (!is.na(match("(QNTM)", file_header[2], nomatch = FALSE))) {
    unit.in <- "photon"
  } else {
    unit.in <- "energy"
  }
  
  if (unit.in == "photon") {
    col_names <- c("w.length", "s.q.irrad")
    mult <- 1e-6 # umol -> mol
  } else if (unit.in == "energy") {
    col_names <- c("w.length", "s.e.irrad")
    mult <- 1.0 # joule -> joule
  }
  
  # does not decode first column correctly if it includes values >= 1000
  # z <-
  #   readr::read_table(file,
  #                     col_names = col_names,
  #                     col_types = "dd",
  #                     skip = 7,
  #                     locale = locale)
  
  z <- 
    utils::read.table(file,
                      header = FALSE,
                      dec = ".",
                      row.names = NULL,
                      col.names = col_names,
                      colClasses = "double",
                      skip = 7
                      )
  if (mult != 1) {
    dots <- list(~s.q.irrad * mult)
    z <- dplyr::mutate_(z, .dots = stats::setNames(dots, "s.q.irrad"))
  }
  
  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.source_spct(z, time.unit = "second")
  options(old.opts)
  comment(z) <-
    paste(paste("LICOR LI-1800 file '", file, "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")
  
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  z
}

#' @rdname read_licor_prn
#' @param files A list or vector of character strings.
#' @export
#' @return Function \code{read_m_licor_prn()} returns a source_mspct object
#'   containing one spectrum per file read.
#'   
#' @details Function \code{read_m_licor_prn()} calls \code{red_licor_file()} 
#'   for each file in \code{files}. See \code{\link[readr]{read_table}} for
#'   a description of valid arguments for \code{files}.
#' 
read_m_licor_prn <- function(files,
                             date = NULL,
                             geocode = NULL,
                             label = NULL,
                             tz = Sys.timezone(location = FALSE),
                             locale = readr::default_locale()) {
  list.of.spectra <- list()
  for (f in files) {
    spct.name <- tolower(sub(".PRN", "", f))
    
    list.of.spectra[[spct.name]] <-
      read_licor_prn(
        file = f,
        date = date,
        geocode = geocode,
        label = label,
        tz = tz,
        locale = locale
      )
  }
  photobiology::source_mspct(list.of.spectra)
}

