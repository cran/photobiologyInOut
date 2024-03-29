#' Read '.csv' File Saved by Avantes' Software for AvaSpec.
#' 
#' Reads and parses the header of a processed data file as output by the
#' program Avaspec and then imports wavelength and spectral irradiance
#' values. The file header has little useful metadata information.
#' 
#' @param file character string
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date is extracted from the
#'   file header.
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
#'   
#' @return A source_spct object.
#' @export
#' @references \url{https://www.avantes.com/}
#' @keywords misc
#' 
#' @examples
#' 
#'  file.name <- 
#'     system.file("extdata", "spectrum-avaspec.csv", 
#'                 package = "photobiologyInOut", mustWork = TRUE)
#'                 
#'  avaspec.spct <- read_avaspec_csv(file = file.name)
#'  
#'  avaspec.spct
#'  getWhatMeasured(avaspec.spct)
#'  cat(comment(avaspec.spct))
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
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  file_header <- scan(file = file, nlines = 6, skip = 0, 
                      what = "character", quiet = TRUE)
  # watt / cm ?
  if (length(grep("Watt/cm", file_header[2], fixed = TRUE))) {
    mult <- 10e-4
  } else {
    mult <- NA
  }
  
  z <- readr::read_csv(file = file,
                       col_names = c("w.length", "s.e.irrad"),
                       skip = 6,
                       col_types = readr::cols(),
                       progress = FALSE,
                       locale = locale)
  z[ , "s.e.irrad"] <- z[ , "s.e.irrad"] * mult

  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.source_spct(z, time.unit = "second")
  options(old.opts)
  
  comment(z) <-
    paste(paste("Avantes AvaSpec irradiance file '", basename(file), "' imported on ", 
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

#' @rdname read_avaspec_csv
#' 
#' @param path Path to the xls/xlsx file
#' 
#' @export
#' 
read_avaspec_xls <- function(path,
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
  
  file_header <- readxl::read_excel(path, skip = 0, col_names = FALSE)[1:4, 1]

  column_names <- readxl::read_excel(path, skip = 4, col_names = FALSE)[1, ]
  column_names <- as.vector(as.matrix(column_names))
  column_names <- gsub("Wave.*", "w.length", column_names)
  column_names <- gsub("Absolute Irradiance.*", "s.e.irrad", column_names)
  column_names <- gsub("Photon Count.*", "s.q.irrad", column_names)

  z <- readxl::read_excel(path,
                          col_names = column_names,
                          progress = FALSE,
                          skip = 6)
  
  z <- z[!is.na(z["w.length"]),
         names(z) %in% c("w.length", "s.e.irrad", "s.q.irrad")]
  
  z["s.e.irrad"] <- z["s.e.irrad"] * 1e-2 # 1e4 * 1e-6
  z["s.q.irrad"] <- z["s.q.irrad"] * 1e-6
  
  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.source_spct(z, time.unit = "second")
  options(old.opts)
  
  comment(z) <-
    paste(paste("Avantes AvaSpec irradiance file (Excel)'", basename(path), "' imported on ", 
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

