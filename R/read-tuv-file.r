#' Read TUV output file.
#' 
#' Reads and parses the header of a text file output by the TUV program to
#' extract the header and spectral data. The time field is converted to a date.
#' 
#' @param file character string
#' @param date a \code{POSIXct} object, but if \code{NULL} the date stored in
#'   file is used, and if \code{NA} no date variable is added
#' @param geocode A data frame with columns \code{lon} and \code{lat}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone is by default read from the file.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#'   
#' @return a source_spct object obtained by 'melting' the TUV file, and adding
#'   a factor \code{spct.idx}, and variables \code{zenith.angle} and
#'   \code{date}.
#'   
#' @references \url{http://www.r4photobiology.info}
#' @keywords misc
#'
#' @note Tested only with TUV versison 5.0.
#' 
#' @export
#' 
read_tuv_usrout <- function(file, 
                            date = lubridate::today(),
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
  file_header <- scan(file = file, nlines = 5, what = "character", sep = "\n" )
  hours <- scan(text = sub(pattern = "wc, nm", replacement = "",
                           x = file_header[4], fixed = TRUE))
  num.spectra <- length(hours)
  
  minutes <- trunc((hours - trunc(hours)) * 60)
  seconds <- trunc((minutes - trunc(minutes)) * 60)

  lubridate::hour(date) <- trunc(hours) 
  lubridate::minute(date) <- trunc(minutes)
  lubridate::second(date) <- trunc(seconds)
  
  angles <- scan(text = sub(pattern = "sza = ", replacement = "", 
                            x = file_header[5], fixed = TRUE))
  
  wide.df <- readr::read_table(file = file, skip = 5, 
                               col_names = c("w.length", LETTERS[1:num.spectra]),
                               locale = locale)
  
  wl.length <- length(wide.df[["w.length"]])

  z <- reshape2::melt(wide.df, id.vars = "w.length", 
                      value.name = "s.e.irrad", variable.name = "spct.idx")
  
  z[["angle"]] <- with(z, rep(angles, rep(wl.length, num.spectra)))
  z[["date"]] <- with(z, rep(as.POSIXct(date), rep(wl.length, num.spectra)))
  
  photobiology::setSourceSpct(z, time.unit = "second", multiple.wl = num.spectra)

  comment(z) <- paste(paste("TUV file '", file, "' imported on ", 
                            lubridate::now(tz = "UTC"), " UTC", sep = ""), 
                      paste(file_header, collapse = "\n"), sep = "\n")
  photobiology::setWhatMeasured(z, paste("TUV spectral simulation", label))
  photobiology::setWhenMeasured(z, geocode)
  setWhenMeasured(z, unique(z[["date"]]))
  z
}

