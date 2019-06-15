#' Read File Saved by Ocean Optics' Raspberry Pi software.
#' 
#' Reads and parses the header of a raw data file as output by the server
#' running on a Raspberry Pi board to extract the whole header remark field. The
#' time field is retrieved and decoded.
#' 
#' @param file character string
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date is set to the file
#'   modification date.
#' @param geocode A data frame with columns \code{lon} and \code{lat} used to
#'   set attribute \code{"where.measured"}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone is not saved to the file.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#' @param npixels integer Number of pixels in spectral data.
#' @param spectrometer.sn character The serial number of the spectrometer needs
#'   to be supplied by the user as it is not included in the file header.
#'   
#' @return A raw_spct object.
#' 
#' @export
#' @references \url{https://www.r4photobiology.info} \url{https://oceanoptics.com/} \url{https://www.raspberrypi.org/}
#' 
#' @note The header in these files has very little information, so the user
#' needs to supply the number of pixels in the array as well as the date-time.
#' The file contains a time in milliseconds but as the Raspberry Pi board
#' contains no real-time clock, it seems to default to number of milliseconds
#' since the Pi was switched on. If no argument is passed to date this
#' attribute is set to the file modification date obtained with file.mtime().
#' This date-time gives an upper limit to the real time of measurement as in
#' some operating systems it is reset when the file is copied or even without
#' any good apparent reason.
#' 
read_oo_pidata <- function(file,
                           date = NULL,
                           geocode = NULL,
                           label = NULL,
                           tz = NULL,
                           locale = readr::default_locale(),
                           npixels = 2048,
                           spectrometer.sn = "FLMS00673") {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  first.line <- scan(file = file, nlines = 1, 
                     skip = 0, what = "character", sep = "\n")
  
  skip.n <- ifelse(grepl("sequence", first.line), 1L, 0L)
  
  file_header <- scan(file = file, nlines = 4, 
                      skip = skip.n, what = "character", sep = "\n")
  NonASCII <- tools::showNonASCII(file_header)
  if (length(NonASCII) > 0L) {
    warning("Found non-ASCII characters in file header: ", 
            NonASCII,
            "replacing with ' '.")
    file_header <- iconv(file_header, to = "ASCII", sub = " ")
  }
  
  if (is.null(date)) {
     # we use file modification time lacking anything better
     date <- file.mtime(file)
  }
  
  integ.time <- as.numeric(sub("Integration time: ", "", file_header[2L]))
  
  num.scans <- as.integer(sub("Scans to average: ", "", file_header[3L])) 
  num.scans <- ifelse(num.scans < 1L, 1L, num.scans)
  
  boxcar.width <- as.integer(sub("Boxcar smoothing: ", "", file_header[4L]))
  boxcar.width <- ifelse(boxcar.width < 1L, 1L, boxcar.width)
  
  inst.descriptor <-
    list(
      time = date,
      w = NULL,
      sr.index = NA_integer_,
      ch.index = NA_integer_,
      spectrometer.name = "OO with Raspberry Pi",
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
      correct.elec.dark =  NA_integer_,
      correct.non.lin = NA_integer_,
      correct.stray.light = NA_integer_,
      boxcar.width = boxcar.width,
      integ.time = integ.time,
      num.scans = num.scans,
      tot.time = integ.time * num.scans
    )
  
  z <- readr::read_tsv(
    file = file,
    col_names = c("w.length", "counts"),
    skip = 5 + skip.n,
    n_max = npixels,
    col_types = readr::cols(),
    locale = locale
  )
  
  z <- photobiology::as.raw_spct(z)

  comment(z) <-
    paste(paste("Ocean Optics Raspeberry Pi raw counts file '", basename(file), "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")

  photobiology::setInstrDesc(z, inst.descriptor)
  photobiology::setInstrSettings(z, inst.settings)
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  how <- "Measured with an array spectrometer."
  photobiology::setHowMeasured(z, how)
  attr(z, "file.header") <- file_header
  z
}

