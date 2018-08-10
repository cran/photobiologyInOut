#' Read TUV output file.
#'
#' Reads and parses the header of a text file output by the TUV program to
#' extract the header and spectral data. The time field is converted to a date.
#'
#' @param file character string
#' @param ozone.du numeric Ozone column in Dobson units.
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
#' @return a source_spct object obtained by 'melting' the TUV file, and adding a
#'   factor \code{spct.idx}, and variables \code{zenith.angle} and \code{date}.
#'
#' @references \url{http://www.r4photobiology.info}
#'   \url{https://www2.acom.ucar.edu/modeling/tuv-download}
#' @keywords misc
#'
#' @note The ozone column value used in the simulation cannot be retrieved from
#' the file. Tested only with TUV version 5.0.
#'
#' @export
#' 
read_tuv_usrout <- function(file, 
                            ozone.du = NULL,
                            date = lubridate::today(),
                            geocode = NULL,
                            label = NULL,
                            tz = NULL,
                            locale = readr::default_locale()) {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  if (is.null(geocode)) {
    geocode <- tibble::tibble(lon = NA_real_, lat = NA_real_)
  }
  label <- paste("File:", basename(file), label)
  
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
                               col_types = readr::cols(),
                               locale = locale)
  
  wl.length <- length(wide.df[["w.length"]])

  z <- tidyr::gather(wide.df, .dots = -c("w.length"), 
                      value = "s.e.irrad", key = "spct.idx")
  
  z[["angle"]] <- with(z, rep(angles, rep(wl.length, num.spectra)))
  z[["date"]] <- with(z, rep(as.POSIXct(date), rep(wl.length, num.spectra)))
  
  photobiology::setSourceSpct(z, time.unit = "second", multiple.wl = num.spectra)

  comment(z) <- paste(paste("TUV file '", basename(file), "' imported on ", 
                            lubridate::now(tzone = "UTC"), " UTC", sep = ""), 
                      paste(file_header, collapse = "\n"), sep = "\n")
  photobiology::setWhatMeasured(z, paste("TUV spectral simulation", label))
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhenMeasured(z, unique(z[["date"]]))
  attr(z, "file.header", file_header)
  z
}

#' @rdname read_tuv_usrout
#' 
#' @export
#' 
read_tuv_usrout2mspct <- function(file, 
                                  ozone.du = NULL,
                                  date = lubridate::today(),
                                  geocode = NULL,
                                  label = NULL,
                                  tz = NULL,
                                  locale = readr::default_locale()) {
  z <- read_tuv_usrout(file = file,
                       ozone.du = ozone.du,
                       date = date,
                       geocode = geocode,
                       label = label,
                       tz = tz,
                       locale = locale)
  photobiology::subset2mspct(z)
}
  
#' Read Quick TUV output file.
#' 
#' Reads and parses the header of a text file output by the Quick TUV on-line
#' web front-end at \url{http://cprm.acom.ucar.edu/Models/TUV/Interactive_TUV/}
#' to extract the header and spectral data. The time field is converted to a
#' date.
#' 
#' @param file character string with the name of a text file.
#' @param ozone.du numeric Ozone column in Dobson units.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone is by default read from the file.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#'   
#' @return a source_spct object obtained by finding the center of wavelength
#'   intervals in the Quick TUV output file, and adding variables \code{zenith.angle} and
#'   \code{date}.
#'   
#' @references \url{http://www.r4photobiology.info} 
#' \url{http://cprm.acom.ucar.edu/Models/TUV/Interactive_TUV/}
#' 
#' @note The ozone column value used in the simulation cannot be retrieved from
#' the file. Tested only with Quick TUV versison 5.2 on 2018-07-30. This 
#' function can be expected to be robust to variations in the position of lines
#' in the imported file and resistent to the presence of extraneous text or
#' even summaries.
#' 
#' @export
#' 
read_qtuv_txt <- function(file, 
                          ozone.du = NULL,
                          label = NULL,
                          tz = NULL,
                          locale = readr::default_locale()) {
  if (is.null(tz)) {
    tz <- locale[["tz"]]
  }
  
  label <- paste("File:", basename(file), label)
  
  # make sure we read whole header even is garbage present
  file_header <- readr::read_lines(file = file, n_max = 100L)
  # find top of header
  header.start.idx <- grep("INPUT PARAMETERS:", file_header, fixed = TRUE)
  if (!length(header.start.idx)) {
    warning("File '", file, "' seems not to be output from Quick TUV Calculator.")
    return(source_spct())
  }
  
  # check that file contains spectral irradiance data
  # "SPECTRAL IRRADIANCE (W m-2 nm-1):" -> top of data
  header.end.idx <- grep("SPECTRAL IRRADIANCE", file_header, fixed = TRUE)
  if (!length(header.end.idx)) {
    warning("File '", file, "' does not contain spectral data.")
    return(source_spct())
  } else {
    data.header.line <- file_header[header.end.idx]
  }
  spct.header.idx <- header.end.idx + 2L

  # trim garbage above and below header
  file_header <- file_header[header.start.idx:header.end.idx]

  # find length of spectral data
  grid.line.idx <- grep("w-grid:", file_header, fixed = TRUE)
  temp <-
    scan(text = file_header[grid.line.idx], 
         what = list(NULL, length = 1L, wl.min = 1, wl.maX = 1))
  length.spct <- temp[["length"]] - 1L # number of wl intervals
  wl.min = temp[["wl.min"]]
  wl.max = temp[["wl.max"]]
  # decode metadata
  # read date
  date.line <- grep("idate =", file_header)
  temp <- scan(text = file_header[date.line],
                    what = list(NULL, NULL, idate = "", NULL, NULL, esfact = 1))
  date <- lubridate::ymd(temp[["idate"]])
  esfact <- temp[["esfact"]]
  
  # "solar zenith angle = " -> angle. Always present either user supplied or calculated
  zenith.angle.line <- file_header[grepl("solar zenith angle", file_header)]
  zenith.angle <- scan(text = sub("solar zenith angle =", "", zenith.angle.line, fixed = TRUE))     
  
  # "  measurement point: index            1  altitude=    0.000000". Always present.
  altitude.line <- file_header[grepl("altitude=", file_header)]
  temp <- unlist(
    scan(text = sub("measurement point: ", "", altitude.line, fixed = TRUE),
       what = list(NULL, index = 1, NULL, alt = 1))
    )
  ground.elevation <- temp["alt"]
  observer.above.ground <- temp["index"]
  # " lat=    0.000000      long=    0.000000      ut=    12.00000    ". 
  # Present for option 1 (user supplied location and time)
  geocode.line <- file_header[grepl("lat=", file_header)]
  if (length(geocode.line)) {
    temp <-
    unlist(
      scan(text = geocode.line, what = list(NULL, lat = 1, NULL, lon = 1, NULL, utc = 1))
    )
    lat <- temp["lat"]
    lon <- temp["lon"]
    hours <- temp["utc"]
    
  } else {
    lat <- lon <- hours <- NA_real_
  }
  
  geocode <- tibble::tibble(lon = lon, lat = lat)
  
  if (!is.na(hours)) {
    minutes <- trunc((hours - trunc(hours)) * 60)
    seconds <- trunc((minutes - trunc(minutes)) * 60)
    
    lubridate::hour(date) <- trunc(hours) 
    lubridate::minute(date) <- trunc(minutes)
    lubridate::second(date) <- trunc(seconds)
  }
  # assemple comment
  comment.txt <- paste(gsub(":", "", data.header.line), "\n",
                       "from file: ", file, 
                       " generated by Quick TUV on", "\n",
                       file.info(file)[["ctime"]],
                       "ozone column (DU) = ", ozone.du, "\n",
                       "zenith angle (degrees) = ", zenith.angle, "\n",
                       "altitude (km)  = ", ground.elevation, "\n",
                       "observer elev. = ", observer.above.ground)
  
  # read spectrum
  spct.tb <-
    readr::read_table2(file, skip = spct.header.idx,
                       col_types = readr::cols(.default = readr::col_double()),
                       col_names = c("w.length.s", "w.length.l",
                                     "s.e.irrad.dir",
                                     "s.e.irrad.diff.down", "s.e.irrad.diff.up",
                                     "s.e.irrad"),
                       n_max = length.spct)
  spct.tb <- stats::na.omit(spct.tb)
  z <-
    photobiology::source_spct(w.length = (spct.tb[["w.length.s"]] + spct.tb[["w.length.l"]]) / 2,
                              s.e.irrad = spct.tb[["s.e.irrad"]],
                              s.e.irrad.dir = spct.tb[["s.e.irrad.dir"]],
                              s.e.irrad.diff.down = spct.tb[["s.e.irrad.diff.down"]],
                              s.e.irrad.diff.up = spct.tb[["s.e.irrad.diff.up"]],
                              angle = zenith.angle,
                              date = rep(as.POSIXct(date), nrow(spct.tb)),
                              comment = comment.txt)
  photobiology::setWhatMeasured(z, paste("Quick TUV spectral simulation", label))
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  attr(z, "file.header", file_header)
  z
}
  
  