#' Read File Saved by Ocean Optics' SpectraSuite.
#' 
#' Reads and parses the header of a processed data file as output by
#' SpectraSuite to extract the whole header remark field. The time field is
#' retrieved and decoded.
#' 
#' @param file character string
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date is extracted from the
#'   file header.
#' @param geocode A data frame with columns \code{lon} and \code{lat} used to
#'   set attribute \code{"where.measured"}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone is by default read from the file.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#'   
#' @return A source_spct object.
#' @export
#' @references \url{http://www.r4photobiology.info} \url{http://oceanoptics.com/}
#' @keywords misc
#' 
read_oo_ssirrad <- function(file,
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
  
  line01 <- scan(file = file, nlines =  1, skip = 0, what="character")
  if (line01[1] != "SpectraSuite") {
    warning("Input file was not created by SpectrSuite as expected: skipping")
    return(source_spct())
  }
  file_header <- scan(file = file, nlines = 16, 
                      skip = 0, what="character", sep = "\n")
  
  npixels <- as.integer(sub("Number of Pixels in Processed Spectrum: ", "", 
                            file_header[16], fixed = TRUE))
  
  if (is.null(date)) {
    line03 <- sub("Date: [[:alpha:]]{3} ", "", file_header[3])
    if (is.null(tz)) {
      tz <- sub("^(.{16})([[:upper:]]{3,4})(.{5})$", "\\2", line03)
      if (nchar(tz) == 4) {
        tz <- sub("S", "", tz)
      }
    }
    date <- lubridate::parse_date_time(line03, "mdHMSy", tz = tz, locale = "C")
  }
  
  z <- readr::read_tsv(
    file = file,
    col_names = c("w.length", "s.e.irrad"),
    skip = 17,
    n_max = npixels,
    col_types = readr::cols(),
    locale = locale
  )
  
  dots <- list(~s.e.irrad * 1e-2) # uW cm-2 nm-1 -> W m-2 nm-1
  z <- dplyr::mutate_(z, .dots = stats::setNames(dots, "s.e.irrad"))

  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.source_spct(z, time.unit = "second")
  options(old.opts)

  comment(z) <-
    paste(paste("Ocean Optics Spectra Suite irradiance file '", file, "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")
  
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  attr(z, "file.header") <- file_header
  z
}

#' @rdname read_oo_ssirrad
#' @return A raw_spct object.
#' @export
#' 
read_oo_ssdata<- function(file,
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
  
  line01 <- scan(file = file, nlines =  1, skip = 0, what="character")
  if (line01[1] != "SpectraSuite") {
    warning("Input file was not created by SpectrSuite as expected: skipping")
    return(raw_spct())
  }
  file_header <- scan(file = file, nlines = 16, 
                      skip = 0, what="character", sep = "\n")
  
  npixels <- as.integer(sub("Number of Pixels in Processed Spectrum: ", "", 
                            file_header[16], fixed = TRUE))
  
  if (is.null(date)) {
    line03 <- sub("Date: [[:alpha:]]{3} ", "", file_header[3])
    if (is.null(tz)) {
      tz <- sub("^(.{16})([[:upper:]]{3,4})(.{5})$", "\\2", line03)
      if (nchar(tz) == 4) {
        tz <- sub("S", "", tz)
      }
    }
    date <- lubridate::parse_date_time(line03, "mdHMSy", tz = tz)
  }
  
  z <- readr::read_tsv(
    file = file,
    col_names = c("w.length", "counts"),
    skip = 17,
    n_max = npixels,
    col_types = readr::cols(),
    locale = locale
  )
  
  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.raw_spct(z)
  options(old.opts)

  comment(z) <-
    paste(paste("Ocean Optics Spectra Suite raw counts file '", basename(file), "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")

  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  attr(z, "file.header") <- file_header
  z
}

