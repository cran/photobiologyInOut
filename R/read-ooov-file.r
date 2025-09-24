#' Read File Saved by Ocean Optics' OceanView.
#' 
#' Reads the spectral data and in addition parses the header of a energy
#' irradiance data file as output by OceanView. OceanView is a program from
#' Ocean Optics used to measure UV, visible and NIR radiation with array
#' spectrometers from the same company. OceanView replaces the no longer 
#' supported SpectraSuite program.
#' 
#' @param file character string Path to the file to be read, following R's
#'   use of forward slashes as separator for folder names.
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date is extracted from the
#'   file header.
#' @param geocode A data frame with columns \code{lon} and \code{lat}, and
#'   optionally \code{address} used to set attribute \code{"where.measured"}.
#' @param label character string to which to set the \code{"what.measured"} 
#'   attribute. If \code{NULL} the value of \code{basename(file)} is used, 
#'   and if \code{NA} the \code{"what.measured"} attribute is not set.
#' @param tz character A time zone recognized by R. If \code{NULL}, the default,
#'    it is extracted from `locale`.
#' @param locale	The locale controls defaults that vary from place to place. 
#'    The default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names. Its value must match that used to write the imported file, which
#'   is not necessarily the default one or the local one.
#' @param range a numeric vector of length two, or any other object for which
#'   function \code{range()} will return a range of wavelengths expressed in
#'   nanometres.
#'   
#' @details The header of the file is first decoded and parsed to extract the
#'   time of data acquisition and serial number of the spectrometer, and to
#'   locate the start of the spectral data. The fields are located by name.
#'   The spectral irradiance is re-expressed in \eqn{W m^{-2} nm^{-1}} and 
#'   returned as an object of class \code{\link[photobiology]{source_spct}} 
#'   with metadata stored in attributes \code{when.measured},
#'   \code{what.measured}, and \code{how.measured} set to values extracted from
#'   the header. The value stored in the \code{how.measured} attribute includes
#'   the User: and Serial Number: values extracted from the file header. 
#'   If an argument is passed to parameter \code{geocode}, its value is
#'   saved in attribute \code{where.measured}. The file header in whole
#'   is copied to attribute \code{file.header}. The object's \code{comment}
#'   always gives a text that includes the file name, time of import, function
#'   name and the version of packages 'photobiology' and 'photobiologyInOut' 
#'   used.
#'
#' @return A \code{source_spct} object with columns \code{w.length} with
#'   wavelengths in nanometres and \code{s.e.irrad} with spectral energy
#'   irradiance in \eqn{W m^{-2} nm^{-1}}, attributes \code{comment},
#'   \code{what.measured}, \code{when.measured}, \code{how.measured},
#'   \code{where.measured} and \code{file.header} containing metadata for
#'   the spectrum.
#'   
#' @export
#' @references \url{https://www.oceanoptics.com/}
#' @keywords misc
#' 
#' @examples
#'
#' # energy spectral irradiance file
#' 
#'  file.name <-
#'    system.file("extdata", "spectrum.OVIrrad", 
#'                package = "photobiologyInOut", mustWork = TRUE)
#'                 
#'  ooov.spct <- 
#'    read_oo_ovirrad(file = file.name,
#'                    locale = readr::locale("en", 
#'                                           decimal_mark = ",",
#'                                           grouping_mark = "",
#'                                           tz = "Europe/Warsaw"))
#'  
#'  ooov.spct
#'  getWhenMeasured(ooov.spct)
#'  getWhatMeasured(ooov.spct)
#'  getHowMeasured(ooov.spct)
#'  cat(comment(ooov.spct))
#' 
#'  ooov_clipped.spct <- 
#'    read_oo_ovirrad(file = file.name,
#'                    locale = readr::locale("en", 
#'                                           decimal_mark = ",",
#'                                           grouping_mark = "",
#'                                           tz = "Europe/Warsaw"),
#'                    range = c(280, NA))
#'  
#'  ooov_clipped.spct
#'  getWhenMeasured(ooov_clipped.spct)
#'  getWhatMeasured(ooov_clipped.spct)
#'  getHowMeasured(ooov_clipped.spct)
#'  cat(comment(ooov_clipped.spct))
#' 
read_oo_ovirrad <- function(file,
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
  if (!(grepl("^Data", line01[1]) && grepl("^from", line01[2]))) {
    warning("Input file was not created by OceanView as expected: skipping")
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
  ln.idx <- which(grepl("^Number of Pixels in Spectrum: ",
                     file_header))
  npixels <- as.integer(sub("Number of Pixels in Spectrum: ", "", 
                            file_header[ln.idx], fixed = TRUE))
  stopifnot("Header parsing failure" = 
              !is.na(npixels) && is.integer(npixels) && length(npixels == 1))
  
  ln.idx <- which(grepl("^Spectrometer: ", file_header))
  sr.sn <- trimws(sub("Spectrometer: ", "", file_header[ln.idx], fixed = TRUE))
  
  ln.idx <- which(grepl("^User: ", file_header))
  sr.user <- trimws(sub("User: ", "",  file_header[ln.idx], fixed = TRUE))
  
  if (is.null(date)) {
    ln.idx <- which(grepl("^Date: ",
                          file_header))
    line02 <- sub("Date: [[:alpha:]]{3} ", "", file_header[ln.idx])
    if (is.null(tz)) {
      tz <- sub("^(.{16})([[:upper:]]{3,4})(.{5})$", "\\2", line02)
      if (nchar(tz) == 4) {
        tz <- sub("S", "", tz)
      }
    }
    date <- lubridate::parse_date_time(line02, "mdHMSy", tz = tz, locale = "C")
  }
  
  to.skip <- which(grepl("^>>>>>Begin", file_header))
  stopifnot("Header parsing failure" = 
              !is.na(to.skip) && is.integer(to.skip) && length(to.skip == 1))
  
  z <- readr::read_table(
    file = file,
    col_names = c("w.length", "s.e.irrad"),
    skip = to.skip,
    n_max = npixels,
    guess_max = npixels,
    col_types = readr::cols(),
    progress = FALSE,
    locale = locale
  )
  
  # dots <- list(~s.e.irrad * 1e-2) # uW cm-2 nm-1 -> W m-2 nm-1
  # z <- dplyr::mutate_(z, .dots = stats::setNames(dots, "s.e.irrad"))

  z[["s.e.irrad"]] <- z[["s.e.irrad"]] * 1e-2 # uW cm-2 nm-1 -> W m-2 nm-1

  old.opts <- options("photobiology.strict.range" = NA_integer_)
  z <- photobiology::as.source_spct(z, time.unit = "second")
  if (!is.null(range)) {
    z <- photobiology::clip_wl(z, range)
  }
  options(old.opts)

  comment(z) <-
    paste("Ocean Optics OceanView irradiance file '", basename(file), 
          "' imported on ", 
          lubridate::round_date(lubridate::now(tzone = "UTC")), " UTC ",
          "with function 'read_oo_ovirrad()'.\n",
          "R packages 'photobiologyInOut' ", 
          utils::packageVersion(pkg = "photobiologyInOut"), 
          " and 'photobiology' ",
          utils::packageVersion(pkg = "photobiology"), 
          " were used.", sep = "")
  
  how.measured <- paste("Measured by user ", sr.user, 
                        " with Ocean Optics spectrometer with s.n. ",
                        sr.sn, " and OceanView software.", sep = "")
  
  photobiology::setHowMeasured(z, how.measured)
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  attr(z, "file.header") <- file_header
  z
}
