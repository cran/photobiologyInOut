#' Read File Saved by Ocean Optics' SpectraSuite.
#' 
#' Reads the spectral data and in addition parses the header of a energy
#' irradiance data file as output by SpectraSuite. SpectraSuite is a program 
#' from Ocean Optics used to measure UV, visible and NIR radiation with array
#' spectrometers from the same company. OceanView replaces the no longer 
#' supported SpectraSuite program.
#' 
#' @inherit read_oo_ovirrad details return references
#' @inheritParams read_oo_ovirrad
#'   
#' @export
#' 
#' @examples
#' 
#'  file.name <- 
#'    system.file("extdata", "spectrum.SSIrrad", 
#'                package = "photobiologyInOut", mustWork = TRUE)
#'                 
#'  ooss.spct <- read_oo_ssirrad(file = file.name)
#'  
#'  ooss.spct
#'  getWhenMeasured(ooss.spct)
#'  getWhatMeasured(ooss.spct)
#'  getHowMeasured(ooss.spct)
#'  cat(comment(ooss.spct))
#' 
read_oo_ssirrad <- function(file,
                            date = NULL,
                            geocode = NULL,
                            label = NULL,
                            tz = NULL,
                            locale = readr::default_locale(),
                            range = NULL) {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  line01 <- scan(file = file, nlines =  1, skip = 0, 
                 what="character", quiet = TRUE)
  if (line01[1] != "SpectraSuite") {
    warning("Input file was not created by SpectrSuite as expected: skipping")
    return(photobiology::source_spct())
  }
  file_header <- scan(file = file, nlines = 20, 
                      skip = 0, what="character", 
                      blank.lines.skip = FALSE, # to get start of data 
                      sep = "\n", quiet = TRUE)
  NonASCII <- tools::showNonASCII(file_header)
  if (length(NonASCII) > 0L) {
    warning("Found non-ASCII characters in file header: ", 
            NonASCII,
            "replacing with ' '.")
    file_header <- iconv(file_header, to = "ASCII", sub = " ")
  }
  
  ln.idx <- which(grepl("^Number of Pixels in Processed Spectrum: ",
                        file_header))
  npixels <- as.integer(sub("Number of Pixels in Processed Spectrum: ", "", 
                            file_header[ln.idx], fixed = TRUE))
  stopifnot("Header parsing failure" = 
              !is.na(npixels) && is.integer(npixels) && length(npixels == 1))
  
  ln.idx <- which(grepl("^Spectrometers: ", file_header))
  sr.sn <- trimws(sub("Spectrometers: ", "", 
                      file_header[ln.idx], fixed = TRUE))
  
  ln.idx <- which(grepl("^User: ", file_header))
  sr.user <- trimws(sub("User: ", "",  file_header[ln.idx], fixed = TRUE))
  
  if (is.null(date)) {
    ln.idx <- which(grepl("^Date: ",
                          file_header))
    line03 <- sub("Date: [[:alpha:]]{3} ", "", file_header[ln.idx])
    if (is.null(tz)) {
      tz <- sub("^(.{16})([[:upper:]]{3,4})(.{5})$", "\\2", line03)
      if (nchar(tz) == 4) {
        tz <- sub("S", "", tz)
      }
    }
    date <- lubridate::parse_date_time(line03, "mdHMSy", tz = tz, locale = "C")
  }
  
  to.skip <- which(grepl("^>>>>>Begin", file_header))
  stopifnot("Header parsing failure" = 
              !is.na(to.skip) && is.integer(to.skip) && length(to.skip == 1))
  
  z <- readr::read_tsv(
    file = file,
    col_names = c("w.length", "s.e.irrad"),
    skip = to.skip,
    n_max = npixels,
    col_types = readr::cols(),
    progress = FALSE,
    locale = locale
  )
  
  z[["s.e.irrad"]] <- z[["s.e.irrad"]] * 1e-2 # uW cm-2 nm-1 -> W m-2 nm-1

  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.source_spct(z, time.unit = "second")
  if (!is.null(range)) {
    z <- photobiology::clip_wl(z, range)
  }
  options(old.opts)

  comment(z) <-
    paste("Ocean Optics SpectraSuite irradiance file '", basename(file), 
          "' imported on ", 
          lubridate::round_date(lubridate::now(tzone = "UTC")), " UTC ",
          "with function 'read_oo_ssirrad()'.\n",
          "R packages 'photobiologyInOut' ", 
          utils::packageVersion(pkg = "photobiologyInOut"), 
          " and 'photobiology' ",
          utils::packageVersion(pkg = "photobiology"), 
          " were used.", sep = "")
  
  how.measured <- paste("Measured by user ", sr.user, 
                        " with Ocean Optics spectrometer with s.n. ",
                        sr.sn, " and SpectraSuite software.", sep = "")
  
  photobiology::setHowMeasured(z, how.measured)
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
  
  line01 <- scan(file = file, nlines =  1, skip = 0, 
                 what="character", quiet = TRUE)
  if (line01[1] != "SpectraSuite") {
    warning("Input file was not created by SpectrSuite as expected: skipping")
    return(photobiology::raw_spct())
  }
  file_header <- scan(file = file, nlines = 16, 
                      skip = 0, what="character", 
                      sep = "\n", quiet = TRUE)
  
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
    progress = FALSE,
    locale = locale
  )
  
  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.raw_spct(z)
  options(old.opts)

  comment(z) <-
    paste(paste("Ocean Optics Spectra Suite raw counts file '", 
                basename(file), "' imported on ", 
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

