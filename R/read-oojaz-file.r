#' Read Files Saved by Ocean Optics' Jaz spectrometer.
#' 
#' Reads and parses the header of processed data text files output by
#' Jaz instruments extracting the spectral data from the body of the file 
#' and the metadata, including time and date of measurement from the header.
#' Jaz modular spectrometers were manufactured by Ocean Optics.
#' The company formerly named Ocean Optics is now called Ocean Insight.
#' 
#' @details Function \code{read_oo_jazirrad} can read processed irradiance
#' output files. Function \code{read_oo_jazpc} can read processed transmittance
#' and reflectance output files (expressed as \%s). Function 
#' \code{read_oo_jazdata} can read raw-counts data.
#' 
#' @param file character string.
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date is extracted from the
#'   file header.
#' @param geocode A data frame with columns \code{lon} and \code{lat} used to
#'   set attribute \code{"where.measured"}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone used for interpreting times saved in the file
#'   header.
#' @param locale The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#'   
#' @note Although the parameter is called \code{date} a date time is accepted 
#'   and expected. Time resolution is < 1 s if seconds are entered with a
#'   decimal fraction, such as "2021-10-05 10:10:10.1234".
#'   
#' @return A source_spct object, a filter_spct object, a reflector_spct object
#'   or a raw_spct object.
#' 
#' @export
#' 
#' @references \url{https://www.oceaninsight.com/}
#' 
#' @examples
#' 
#'  file.name <- 
#'    system.file("extdata", "spectrum.jaz", 
#'                package = "photobiologyInOut", mustWork = TRUE)
#'                 
#'  jaz.spct <- read_oo_jazpc(file = file.name)
#'  
#'  jaz.spct
#'  getWhenMeasured(jaz.spct)
#'  getWhatMeasured(jaz.spct)
#'  cat(comment(jaz.spct))
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
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  line01 <-
    scan(
      file = file,
      nlines =  1,
      skip = 0,
      what = "character", 
      quiet = TRUE
    )
  if (line01[1] != "Jaz") {
    warning("Input file was not created by a Jaz spectrometer as expected: skipping")
    return(photobiology::source_spct())
  }
  if (line01[3] != "Irradiance") {
    warning("Input file does not contain data labeled as 'Irradiance' as expected: skipping")
    return(photobiology::source_spct())
  }
  file_header <-
    scan(
      file = file,
      nlines = 18,
      skip = 0,
      what = "character",
      sep = "\n", 
      quiet = TRUE
    )
  NonASCII <- tools::showNonASCII(file_header)
  if (length(NonASCII) > 0L) {
    warning("Found non-ASCII characters in file header: ", 
            NonASCII,
            "replacing with ' '.")
    file_header <- iconv(file_header, to = "ASCII", sub = " ")
  }
  
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
    col_types = readr::cols(),
    progress = FALSE,
    locale = locale
  )
  dots <- list(~W, ~P)
  z <- dplyr::select(z, 
                     w.length = "W",
                     s.e.irrad = "S")
  z[["s.e.irrad"]] <- z[["s.e.irrad"]] * 1e-2 # uW cm-2 nm-1 -> W m-2 nm-1
  
  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.source_spct(z, time.unit = "second")
  options(old.opts)

  comment(z) <-
    paste(paste("Ocean Optics Jaz irradiance file '", basename(file), "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  how <- "Measured with an array spectrometer."
  photobiology::setHowMeasured(z, how)
  attr(z, "file.header") <- file_header
  z
}

#' @rdname read_oo_jazirrad
#'
#' @param qty.in character string, one of "Tpc" (spectral transmittance, \%), "A"
#'   (spectral absorbance), or "Rpc" (spectral reflectance, \%).
#' @param Tfr.type character string, either "total" or "internal".
#' @param Rfr.type character string, either "total" or "specular".
#'
#' @export
#' 
read_oo_jazpc <- function(file,
                          qty.in = "Tpc",
                          Tfr.type = c("total", "internal"),
                          Rfr.type = c("total", "specular"),
                          date = NULL,
                          geocode = NULL,
                          label = NULL,
                          tz = NULL,
                          locale = readr::default_locale()) {
  stopifnot(qty.in %in% c("Tpc", "Rpc", "A"))
  
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  line01 <-
    scan(
      file = file,
      nlines =  1,
      skip = 0,
      what = "character", 
      quiet = TRUE
    )
  if (line01[1] != "Jaz") {
    warning("Input file was not created by a Jaz spectrometer as expected: skipping")
    return(photobiology::source_spct())
  }
  file_header <-
    scan(
      file = file,
      nlines = 16,
      skip = 0,
      what = "character",
      sep = "\n", 
      quiet = TRUE
    )
  
  npixels <- as.integer(sub("Number of Pixels in Processed Spectrum: ", "", 
                            file_header[16], fixed = TRUE))
  
  if (is.null(date)) {
    line03 <- sub("Date: [[:alpha:]]{3} ", "", file_header[3])
    date <-
      lubridate::parse_date_time(line03, "mdHMSy", tz = tz)
  }
  
  #  data_header <- scan(file = file, nlines = 1, skip = 20, what = "character")
  
  z <- readr::read_tsv(
    file = file,
    col_names = TRUE,
    skip = 17,
    n_max = npixels,
    col_types = readr::cols(),
    progress = FALSE,
    locale = locale
  )

  z <- z[ , c("W", "P")]
  colnames(z) <- c("w.length", qty.in)

  old.opts <- options("photobiology.strict.range" = NA_integer_)
  if (qty.in == "Rpc") {
    z <- photobiology::as.reflector_spct(z, Rfr.type = Rfr.type)
  } else {
    z <- photobiology::as.filter_spct(z, Tfr.type = Tfr.type)
  }
  options(old.opts)
  
  comment(z) <-
    paste(paste("Ocean Optics Jaz file '", basename(file), "' imported as ", qty.in, " on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  attr(z, "file.header") <- file_header
  z
}

#' @rdname read_oo_jazirrad
#' 
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
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else {
    label <- paste(label.file, label, sep = "\n")
  }
  
  line01 <-
    scan(
      file = file,
      nlines =  1,
      skip = 0,
      what = "character", 
      quiet = TRUE
    )
  if (line01[1] != "Jaz") {
    warning("Input file was not created by a Jaz spectrometer as expected: skipping")
    return(photobiology::source_spct())
  }
  if (line01[2] != "Data") {
    warning("Input file does not contain data labeled as 'Data' as expected: skipping")
    return(photobiology::source_spct())
  }
  file_header <-
    scan(
      file = file,
      nlines = 16,
      skip = 0,
      what = "character",
      sep = "\n", 
      quiet = TRUE
    )
  
  npixels <- as.integer(sub("Number of Pixels in Processed Spectrum: ", "", 
                            file_header[16], fixed = TRUE))
  
  if (is.null(date)) {
    line03 <- sub("Date: [[:alpha:]]{3} ", "", file_header[3])
    date <-
      lubridate::parse_date_time(line03, "mdHMSy", tz = tz)
  }
  
  spectrometer.sn <- sub("Spectrometers: ", "", file_header[8])
  sn.tag <- paste(" \\(", spectrometer.sn, "\\)", sep = "")
  yes.no <- c(Yes = 1L, No = 0L)

  integ.time <- as.numeric(
    sub(sn.tag, "", sub("Integration Time \\(usec\\): ", "", 
                        file_header[9])))
  
  num.scans <- as.integer(
    sub(sn.tag, "", sub("Spectra Averaged: ", "", 
                        file_header[10])))
  num.scans <- ifelse(num.scans < 1L, 1L, num.scans)
  
  boxcar.width <- as.integer( 
    sub(sn.tag, "", sub("Boxcar Smoothing: ", "", 
                        file_header[11])))
  boxcar.width <- ifelse(boxcar.width < 1L, 1L, boxcar.width)

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
      boxcar.width = boxcar.width,
      integ.time = integ.time,
      num.scans = num.scans,
      tot.time = integ.time * num.scans
    )
  
  #  data_header <- scan(file = file, nlines = 1, skip = 20, what = "character")
  
  z <- readr::read_tsv(
    file = file,
    col_names = TRUE,
    skip = 17,
    n_max = npixels,
    col_types = readr::cols(),
    progress = FALSE,
    locale = locale
  )

  z <- dplyr::select(z, 
                     w.length = "W",
                     counts = "S")
  
  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.raw_spct(z)
  options(old.opts)
  
  comment(z) <-
    paste(paste("Ocean Optics Jaz raw counts file '", basename(file), "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")

  attr(z, "linearized") <- inst.settings$correct.non.lin
  photobiology::setInstrDesc(z, inst.descriptor)
  photobiology::setInstrSettings(z, inst.settings)
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  attr(z, "file.header") <- file_header
  z
}
