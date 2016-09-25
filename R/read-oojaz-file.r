#' Read Absolute Irradiance File Saved by Ocean Optics' Jaz spectrometer.
#' 
#' Reads and parses the header of a processed data file as output by
#' Jaz instruments to extract the whole header remark field The time field is
#' retrieved.
#' 
#' @param file character string.
#' @param date a \code{POSIXct} object, but if \code{NULL} the date stored in 
#'   file header is used, and if \code{NA} the "when.measured" attribute is not
#'   set.
#' @param geocode A data frame with columns \code{lon} and \code{lat}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone used for interpreting times saved in the file
#'   header.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#'   
#' @note Although the parameter is called \code{date} a date time is accepted 
#'   and expected. Time resolution is 1 s.
#'   
#' @return A source_spct object.
#' @export
#' @references \url{http://www.r4photobiology.info} \url{http://oceanoptics.com/}
#' @keywords misc
#' 
read_oo_jazirrad <- function(file,
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
  line01 <-
    scan(
      file = file,
      nlines =  1,
      skip = 0,
      what = "character"
    )
  if (line01[1] != "Jaz") {
    warning("Input file was not created by a Jaz spectrometer as expected: skipping")
    return(source_spct())
  }
  if (line01[3] != "Irradiance") {
    warning("Input file does not contain data labeled as 'Irradiance' as expected: skipping")
    return(source_spct())
  }
  file_header <-
    scan(
      file = file,
      nlines = 18,
      skip = 0,
      what = "character",
      sep = "\n"
    )
  
  npixels <- as.integer(sub("Number of Pixels in Processed Spectrum: ", "", 
                            file_header[15], fixed = TRUE))
  
  if (is.null(date)) {
    line03 <- sub("Date: [[:alpha:]]{3} ", "", file_header[3])
    date <-
      lubridate::parse_date_time(line03, "mdHMSy", tz = tz)
  }
  
  #  data_header <- scan(file = file, nlines = 1, skip = 20, what = "character")
  
  z <- readr::read_tsv(
    file = file,
    col_names = TRUE,
    skip = 19,
    n_max = npixels,
    locale = locale
  )
  dots <- list(~W, ~P)
  z <- dplyr::select_(z, .dots = stats::setNames(dots, c("w.length", "s.e.irrad")))
  dots <- list(~s.e.irrad * 1e-2) # uW cm-2 nm-1 -> W m-2 nm-1
  z <- dplyr::mutate_(z, .dots = stats::setNames(dots, "s.e.irrad"))
  
  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.source_spct(z, time.unit = "second")
  options(old.opts)

  comment(z) <-
    paste(paste("Ocean Optics Jaz irradiance file '", file, "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")
    photobiology::setWhenMeasured(z, date)
    photobiology::setWhereMeasured(z, geocode)
    photobiology::setWhatMeasured(z, label)
  z
}

#' @rdname read_oo_jazirrad
#' @return A raw_spct object.
#' @export
#' 
read_oo_jazdata <- function(file,
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
  line01 <-
    scan(
      file = file,
      nlines =  1,
      skip = 0,
      what = "character"
    )
  if (line01[1] != "Jaz") {
    warning("Input file was not created by a Jaz spectrometer as expected: skipping")
    return(source_spct())
  }
  if (line01[2] != "Data") {
    warning("Input file does not contain data labeled as 'Data' as expected: skipping")
    return(source_spct())
  }
  file_header <-
    scan(
      file = file,
      nlines = 16,
      skip = 0,
      what = "character",
      sep = "\n"
    )
  
  npixels <- as.integer(sub("Number of Pixels in Processed Spectrum: ", "", 
                            file_header[16], fixed = TRUE))
  
  if (is.null(date)) {
    line03 <- sub("Date: [[:alpha:]]{3} ", "", file_header[3])
    date <-
      lubridate::parse_date_time(line03, "mdHMSy", tz = tz)
  }
  
  spectrometer.sn <- sub("Spectrometers: ", "", file_header[8])
  
  inst.descriptor <-
    list(
      time = date,
      w = NULL,
      sr.index = NA_integer_,
      ch.index = NA_integer_,
      spectrometer.name = line01[1],
      spectrometer.sn = spectrometer.sn,
      bench.grating = NA_character_,
      bench.filter = NA_character_,
      bench.slit = NA_character_,
      min.integ.time = NA_integer_,
      max.integ.time = NA_integer_,
      max.counts = NA_integer_,
      wavelengths = NA_real_,
      bad.pixs = numeric(),
      inst.calib = list()
    )
  
  yes.no <- c(Yes = 1L, No = 0L)
  sn.tag <- paste(" \\(", spectrometer.sn, "\\)", sep = "")
  inst.settings <- 
    list(
      time = date,
      w = NULL,
      sr.index = NA_integer_,
      ch.index = NA_integer_,
      correct.elec.dark =  yes.no[
        sub(sn.tag, "", sub("Correct for Electrical Dark: ", "", file_header[12]))],
      correct.non.lin = yes.no[
        sub(sn.tag, "", sub("Correct for Detector Non-linearity: ",
            "", file_header[14]))],
      correct.stray.light = yes.no[
        sub(sn.tag, "", sub("Correct for Stray Light: ", "", 
            file_header[15]))],
      boxcar.width = 
        sub(sn.tag, "", sub("Boxcar Smoothing: ", "", 
            file_header[11])),
      integ.time = 
        sub(sn.tag, "", sub("Integration Time \\(usec\\): ", "", 
            file_header[9])),
      num.scans = 
        sub(sn.tag, "", sub("Spectra Averaged: ", "", 
            file_header[10]))
    )
  
  #  data_header <- scan(file = file, nlines = 1, skip = 20, what = "character")
  
  z <- readr::read_tsv(
    file = file,
    col_names = TRUE,
    skip = 17,
    n_max = npixels,
    locale = locale
  )
  dots <- list(~W, ~S)
  z <- dplyr::select_(z, .dots = stats::setNames(dots, c("w.length", "counts")))
  
  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.raw_spct(z)
  options(old.opts)
  
  comment(z) <-
    paste(paste("Ocean Optics Jaz raw counts file '", file, "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")

  attr(z, "linearized") <- inst.settings$correct.non.lin
  photobiology::setInstrDesc(z, inst.descriptor)
  photobiology::setInstrSettings(z, inst.settings)
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  z
}
