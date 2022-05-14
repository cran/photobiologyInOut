#' Read '.TXT' File(s) Saved by LI-COR's LI-180 spectroradiometer.
#' 
#' Reads and parses the header of a data file as output by the LI-180
#' spectrometer (not to be confused with the LI-1800 spectrometer released in
#' the 1980's by LI-COR) to extract the whole header remark field and also
#' decode whether data is in photon or energy based units. This is a new
#' instrument released in year 2020.
#' 
#' @param file Path to file as a character string.
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
#' @param s.qty character The name of the spectral quantity to be read. One of
#'   "s.e.irrad" or "s.q.irrad".
#'   
#' @return \code{read_licor_espd()} returns a \code{source_spct} object with
#'   \code{time.unit} attribute set to \code{"second"} and \code{when.measured}
#'   attribute set to the date-time extracted from the file header, or supplied
#'   by the user. Spectrometer model, serial number and integration time are
#'   stored in attributes. The whole file header is saved as a \code{comment}
#'   while the footer is discarded.
#'   
#' @export
#' 
#' @references LI-COR Biosciences, Environmental.
#'   \url{https://www.licor.com/env/}
#' 
#' @note The LI-180 spectroradiometer stores little information of the
#'   instrument and settings, possibly because they cannot be altered by the
#'   user or configured. The length of the file header does not seem to be
#'   fixed, so the start of the spectral data is detected by searching for
#'   "380nm".
#'   
#' @examples
#' 
#'   file.name <- 
#'     system.file("extdata", "LI-180-irradiance.txt", 
#'                 package = "photobiologyInOut", mustWork = TRUE)
#'                 
#'   licor180.spct <- read_li180_txt(file = file.name)
#'   
#'   licor180.spct
#'   getWhenMeasured(licor180.spct)
#'   getWhatMeasured(licor180.spct)
#'   cat(comment(licor180.spct))
#'   
read_li180_txt <- function(file,
                           date = NULL,
                           geocode = NULL,
                           label = NULL,
                           tz = NULL,
                           locale = readr::default_locale(),
                           s.qty = "s.e.irrad") {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  file_header <- scan(
    file = file,
    nlines = 50,
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
  
  if (!grepl("LI-180$", file_header[1])) {
    warning("File '", file, 
            "' lacks a header with 'Model Name	LI-180' as first line.")
  }
  
  # find first spectral data line
  first.data.line <- which(grepl("^380nm", file_header))
  if (length(first.data.line) != 1L) {
    stop("No spectral data found in file: ", file)
    return(photobiology::source_spct())
  }
  file_header <- file_header[1:(first.data.line - 1L)]
  
  if (is.null(date)) {
    date.char <- sub("Time  ", "", 
                     file_header[grepl("^Time", file_header)])
    date <- lubridate::parse_date_time(date.char, "ymdHMS", tz = tz)
  }
  
  instr.desc <- 
    list(
      spectrometer.name =
        sub("Model Name	", "", 
            file_header[grepl("^Model Name	", file_header)]),
      spectrometer.sn =
        sub("Serial Number	", "", 
            file_header[grepl("^Serial Number	", file_header)]),
      bench.grating = NA_character_,
      bench.slit = NA_character_
    )
  
  instr.settings <- 
    list(
      integ.time =
        as.numeric(
          sub("I-Time	", "", 
              file_header[grepl("^I-Time	", file_header)])
        ) * 1e3,
      tot.time = NA_real_,
      num.scans = NA_integer_,
      rel.signal = NA_real_
    )

  col_names <- c("w.length", s.qty)
 
  z <- 
    utils::read.table(file,
                      header = FALSE,
                      dec = ".",
                      row.names = NULL,
                      col.names = col_names,
                      colClasses = c("character", "double"),
                      skip = first.data.line - 1L,
                      nrows = 401
                      )

  z[["w.length"]] <- as.numeric(gsub("nm", "", z[["w.length"]]))
  z[[s.qty]] <- z[[s.qty]] * 1e-3
  
  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.source_spct(z)
  options(old.opts)
  
  comment(z) <-
    paste(paste("LICOR LI-180 file '", basename(file), "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          paste(file_header, collapse = "\n"), 
          sep = "\n")
  
  photobiology::setInstrDesc(z, instr.desc)
  photobiology::setInstrSettings(z, instr.settings)
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  how <- "Measured with an array spectroradiomete."
  photobiology::setHowMeasured(z, how)
  attr(z, "file.header") <- file_header
  z
}

#' @rdname read_li180_txt
#' @param files A list or vector of character strings.
#' @export
#' @return Function \code{read_m_licor_espd()} returns a source_mspct object
#'   containing one spectrum per file read.
#'   
#' @details Function \code{read_m_licor_espd()} calls \code{red_licor_espd()} 
#'   for each file in \code{files}. See \code{\link{read.table}} for
#'   a description of valid arguments for \code{files}.
#' 
read_m_li180_txt <- function(files,
                             date = NULL,
                             geocode = NULL,
                             label = NULL,
                             tz = Sys.timezone(),
                             locale = readr::default_locale(),
                             s.qty = NULL) {
  list.of.spectra <- list()
  for (f in files) {
    spct.name <- tolower(sub(".PRN", "", f))
    
    list.of.spectra[[spct.name]] <-
      read_li180_txt(
        file = f,
        date = date,
        geocode = geocode,
        label = label,
        tz = tz,
        locale = locale,
        s.qty
      )
  }
  photobiology::generic_mspct(list.of.spectra, 
                              class = class(list.of.spectra[[1]]))
}

