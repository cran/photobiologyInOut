#' Read TUV output file.
#'
#' Reads and parses the header of a text file output by the TUV program to
#' extract the header and spectral data. The time field is converted to a date.
#'
#' @param file character string
#' @param ozone.du numeric Ozone column in Dobson units.
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
#' @return a source_spct object obtained by 'melting' the TUV file, and adding a
#'   factor \code{spct.idx}, and variables \code{zenith.angle} and \code{date}.
#'
#' @references 
#' \url{https://www2.acom.ucar.edu/modeling/tropospheric-ultraviolet-and-visible-tuv-radiation-model}
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

  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  file_header <- scan(file = file, nlines = 5, what = "character", 
                      sep = "\n", quiet = TRUE)
  NonASCII <- tools::showNonASCII(file_header)
  if (length(NonASCII) > 0L) {
    warning("Found non-ASCII characters in file header: ", 
            NonASCII,
            "replacing with ' '.")
    file_header <- iconv(file_header, to = "ASCII", sub = " ")
  }
  
  hours <- scan(text = sub(pattern = "wc, nm", replacement = "",
                           x = file_header[4], fixed = TRUE), 
                quiet = TRUE)
  num.spectra <- length(hours)
  
  minutes <- trunc((hours - trunc(hours)) * 60)
  seconds <- trunc((minutes - trunc(minutes)) * 60)

  lubridate::hour(date) <- trunc(hours) 
  lubridate::minute(date) <- trunc(minutes)
  lubridate::second(date) <- trunc(seconds)
  date <- as.POSIXct(date)
  
  angles <- scan(
    text = sub(pattern = "sza = ", 
               replacement = "", 
               x = file_header[5], fixed = TRUE), 
    quiet = TRUE)
  
  wide.df <- readr::read_table(
    file = file, skip = 5, 
    col_names = c("w.length", LETTERS[1:num.spectra]),
    col_types = readr::cols(),
    progress = FALSE,
    locale = locale)
  
  wl.length <- length(wide.df[["w.length"]])

  z <- tidyr::pivot_longer(
    data = wide.df,
    cols = tidyselect::all_of(setdiff(colnames(wide.df), "w.length")),
    names_to = "spct.idx",
    values_to = "s.e.irrad") 
  z <- z[order(z[["spct.idx"]]), ]
  
  z[["angle"]] <- rep(angles, rep(wl.length, num.spectra))
  z[["date"]] <- rep(as.POSIXct(date), rep(wl.length, num.spectra))
  
  photobiology::setSourceSpct(z, 
                              time.unit = "second", 
                              multiple.wl = num.spectra)

  comment(z) <- paste(paste("TUV file '", 
                            basename(file), 
                            "' imported on ", 
                            lubridate::now(tzone = "UTC"), 
                            " UTC", 
                            sep = ""), 
                      paste(file_header, 
                            collapse = "\n"), 
                      sep = "\n")
  photobiology::setWhatMeasured(z, paste("TUV spectral simulation", label))
  photobiology::setWhereMeasured(z, geocode, simplify = TRUE)
  photobiology::setWhenMeasured(z, unique(z[["date"]]))
  attr(z, "file.header") <- file_header
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
#' web front-end at UCAR to extract the header and spectral data. The time field
#' is converted to a date.
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
#' @param added.vars character vector Accepted member values are \code{"sun.elevation"},
#'   \code{zenith.angle}, \code{"time"} and \code{"ozone.du"}.
#'   
#' @return a source_spct object obtained by finding the center of wavelength
#'   intervals in the Quick TUV output file, and adding the variables listed
#'   in \code{added.vars}. To obtain the same value as in version (<= 0.4.28)
#'   pass \code{added.vars = c("angle", "date")} in the call.
#'   
#' @references \url{https://www.acom.ucar.edu/Models/TUV/Interactive_TUV/}
#' 
#' @note The ozone column value used in the simulation cannot be retrieved from
#' the file. Tested files from Quick TUV version 5.2 on 2018-07-30 and also
#' with more recent files in early 2024. This 
#' function can be expected to be robust to variations in the position of lines
#' in the imported file and resistant to the presence of extraneous text or
#' even summaries. By default web browsers save the output returned by the
#' Quick TUV calculator as an HTML output, some of them with minimal headers
#' and other with more extensive ones. In some cases, character escapes replace
#' actual new lines. In most cases these HTML files are decoded correctly, but
#' if not, use "save as" in the browser and select "text" when saving. As a
#' last recourse, messed up files can be manually edited before import.
#' 
#' @export
#' 
read_qtuv_txt <- function(file, 
                          ozone.du = NULL,
                          label = NULL,
                          tz = NULL,
                          locale = readr::default_locale(),
                          added.vars = NULL) {
  if (is.null(tz)) {
    tz <- locale[["tz"]]
  }
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label) || label == "") {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = " ")
  }
  
  # make sure we read whole header even is garbage present
  file_header <- readr::read_lines(file = file, n_max = 60L)
  
  # find top of file header
  header.start.idx <- grep("INPUT PARAMETERS:", file_header, fixed = TRUE)
  if (!length(header.start.idx)) {
    warning("File '", file, "' seems not to be output from Quick TUV Calculator.")
    return(photobiology::source_spct())
  }
  
  # check that file contains spectral irradiance data
  # "SPECTRAL IRRADIANCE (W m-2 nm-1):" -> top of data
  header.end.idx <- grep("SPECTRAL IRRADIANCE", file_header, fixed = TRUE)
  if (!length(header.end.idx)) {
    warning("File '", file, "' does not contain spectral data.")
    return(photobiology::source_spct())
  } else {
    data.header.line <- 
      file_header[[header.end.idx]] |>
      gsub("\\\\n|\\n", "", x = _) |>
      trimws(x = _)
  }

  # find data column headers
  # used later but should have whole fil_header as argument
  spct.header.idx <- grep("LOWER WVL", file_header, fixed = TRUE)
  
  # trim garbage above and below header
  file_header <- file_header[header.start.idx:header.end.idx]

  # read wavelength grid data
  wgrid.line.idx <- grep("w-grid:", file_header, fixed = TRUE)
  temp <-
    scan(text = file_header[wgrid.line.idx], 
         what = list(NULL, length = 1L, wl.min = 1, wl.max = 1), 
         quiet = TRUE)
  length.spct <- temp[["length"]] - 1L # number of wl intervals
  wl.min = temp[["wl.min"]]
  wl.max = temp[["wl.max"]]
  
  # read vertical grid data
  zgrid.line.idx <- grep("z-grid:", file_header, fixed = TRUE)
  temp <-
    scan(text = file_header[zgrid.line.idx], 
         what = list(NULL, NULL, z.min = 1, z.max = 1), 
         quiet = TRUE)
  observer.km.asl <- temp["z.min"]

  # read date
  date.line <- grep("idate =", file_header)
  temp <- scan(text = file_header[date.line],
               what = list(NULL, NULL, idate = "", NULL, NULL, esfact = 1), 
               quiet = TRUE)
  date <- lubridate::ymd(temp[["idate"]], tz = tz)
#  esfact <- temp[["esfact"]]
  
  # "solar zenith angle = " -> angle. Always present either user supplied or calculated
  zenith.angle.line <- file_header[grepl("solar zenith angle", file_header)]
  zenith.angle <- scan(text = sub("solar zenith angle =", "", 
                                  zenith.angle.line, fixed = TRUE), 
                       quiet = TRUE)     
  
  # "  measurement point: index            1  altitude=    0.000000". Always present.
  altitude.line <- file_header[grepl("altitude=", file_header)]
  temp <- unlist(
    scan(text = sub("measurement point: ", "", altitude.line, fixed = TRUE),
         what = list(NULL, index = 1, NULL, alt = 1), 
         quiet = TRUE)
    )
  ground.km.asl <- temp["alt"]

    # " lat=    0.000000      long=    0.000000      ut=    12.00000    ". 
  # Present for option 1 (user supplied location and time)
  geocode.line <- file_header[grepl("lat=", file_header)]
  if (length(geocode.line)) {
    temp <-
    unlist(
      scan(text = geocode.line, 
           what = list(NULL, lat = 1, NULL, lon = 1, NULL, utc = 1), 
           quiet = TRUE)
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
  date <- as.POSIXct(date)
  
  # assemble comment
  comment.txt <- paste(data.header.line, "\n",
                       "computed by UCAR's Quick TUV calculator\n",
                       "read from file: ", basename(file), "\n",
                       "with modification date: ",
                       lubridate::round_date(file.info(file)[["mtime"]],
                                             unit = "second"), "\n",
                       "ozone column: ", ozone.du, " DU\n",
                       "sun elevation: ", 90 - zenith.angle, " degrees\n",
                       "ground altitude: ", ground.km.asl, " km a.s.l.\n",
                       "observer altitude: ", observer.km.asl, " km a.s.l.",
                       sep = "")

  # read spectrum
  spct.tb <-
    readr::read_table(file, skip = spct.header.idx,
                      col_types = readr::cols(.default = readr::col_double()),
                      col_names = c("w.length.s", "w.length.l",
                                    "s.e.irrad.dir",
                                    "s.e.irrad.diff.down", "s.e.irrad.diff.up",
                                    "s.e.irrad"),
                      progress = FALSE,
                      n_max = length.spct)
  spct.tb <- stats::na.omit(spct.tb)
  # convert to spectrum object
  z <-
    photobiology::source_spct(w.length = (spct.tb[["w.length.s"]] + spct.tb[["w.length.l"]]) / 2,
                              s.e.irrad = spct.tb[["s.e.irrad"]],
                              s.e.irrad.dir = spct.tb[["s.e.irrad.dir"]],
                              s.e.irrad.diff.down = spct.tb[["s.e.irrad.diff.down"]],
                              s.e.irrad.diff.up = spct.tb[["s.e.irrad.diff.up"]],
                              comment = comment.txt)
  # add optional columns with values of input params
  if ("time" %in% added.vars) {
    z[["time"]] <- rep(as.POSIXct(date), nrow(z))
  }
  if ("sun.elevation" %in% added.vars) {
    z[["sun.elevation"]] <- rep(90 - zenith.angle, nrow(z))
  }
  if ("zenith.angle" %in% added.vars) {
    z[["zenith.angle"]] <- rep(zenith.angle, nrow(z))
  }
  if ("ozone.du" %in% added.vars) {
    z[["ozone.du"]] <- rep(ozone.du, nrow(z))
  }
  if ("angle" %in% added.vars) {
    z[["angle"]] <- rep(zenith.angle, nrow(z))
  }
  if ("date" %in% added.vars) {
    z[["date"]] <- rep(as.POSIXct(date), nrow(z))
  }
  # add metadata
  photobiology::setWhatMeasured(z, paste("Solar spectrum (model simulation).", label))
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  how <- "Computer simulation with the TUV model version 5.3"
  photobiology::setHowMeasured(z, how)
  attr(z, "file.header") <- file_header
  z
}

#' Spectral irradiance from the Quick TUV calculator
#'
#' Call the Quick TUV calculator web server and return a \code{source_spct}
#' object with the simulated spectral energy irradiance data.
#'
#' @param w.length list of parameters describing the wavelengths, or a numeric
#'   vector from which the parameters will be constructed.
#' @param sun.elevation numeric Angle in degrees above the horizon. If NULL its
#'   value is computed from \code{geocode} and \code{time}, otherwise 
#'   arguments passed to these two parameters are ignored.
#' @param geocode data.frame with variables \code{lon} and \code{lat} as numeric values
#'   (degrees), and character variable \code{address}; nrow > 1, allowed for collections.
#' @param time A "vector" of POSIXct time, with any valid time zone (TZ) is
#'   allowed, default is current time.
#' @param tz character Time zone is by default read from the file.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#' @param ozone.du numeric Ozone column in Dobson units.
#' @param albedo numeric Surface albedo (= reflectance) as a fraction of one.
#' @param measurement.altitude,ground.altitude numeric Altitudes above sea level
#'   expressed in km.
#' @param clouds data.frame Parameters \code{optical.depth} (vertical), \code{top} and 
#'   \code{base} expressed in km; nrow > 1, allowed for collections.
#' @param aerosols data.frame Parameters \code{optical.depth} (total extinction), \code{ssaaer} 
#'   (cloud single scattering albedo) and \code{alpha} (wavelength dependence of 
#'   optical depth); nrow > 1, allowed for collections.
#' @param num.streams integer Number of streams used in computations, 2 or 4.
#' @param spectra named list with weights for the different components of the 
#'   spectrum.
#' @param added.vars character vector Accepted member values are \code{"sun.elevation"},
#'   \code{zenith.angle}, \code{"date"} and \code{"ozone.du"}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param server.url character The URL used to access the Quick TUV calculator
#'   server.
#' @param file character The name under which the file returned by the server is
#'   locally saved. If \code{NULL} a temporary file is used and discarded
#'   immediately. File paths are supported when valid.
#'   
#' @return In the case of \code{qtuv_s.e.irrad()}, a source_spct object obtained
#'   by finding the center of each wavelength interval in the Quick TUV output
#'   file, and adding the variables listed in \code{added.vars}. In the case of
#'   \code{qtuv_m_s.e.irrad()}, a source_mspct object containing a collection of
#'   such spectra.
#'   
#' @section Side effect: If a file name is passed as argument, the data as
#'   downloaded are saved into persistent files, one file per spectrum. The
#'   names of the saved files always end in `.txt`.
#'
#' @references \url{https://www.acom.ucar.edu/Models/TUV/Interactive_TUV/}. This
#'   URL could change in the future as well as the server URL. The formal
#'   parameter \code{server.url} was included only for use in such a case.
#'
#' @details The Quick TUV calculator, is an on-line freely accessible server
#'   running the TUV atmospheric chemistry and radiation transfer model with a
#'   simplified user interface. In this case, version 5.3 is called passing the
#'   parameter values passed as arguments in the call to
#'   \code{qtuv_s.e.irrad()}. The response is saved in a text file that is
#'   subsequently read with function \code{read_qtuv_txt()} into a
#'   \code{source_spct} object.
#'   
#'   Function \code{qtuv_m_s.e.irrad()} calls \code{qtuv_s.e.irrad()} repeatedly 
#'   accepting a numeric vector longer than one as argument, for parameters: 
#'   \code{sun.elevation}, \code{time} or \code{ozone.du}, and data frames with
#'   nrow >= 1. In a given call, only one parameter at a time can obey multiple 
#'   values, with others currently truncated to the first value.
#'
#'   The formal parameter names are informative and consistent with other 
#'   functions in the R for Photobiology Suite and differ from the short names
#'   used for the parameters in the FORTRAN code of the TUV model. In the case
#'   of \code{w.length} two ways of specifying wavelengths are supported. Some
#'   defaults also differ from those of the Quick TUV calculator.
#'   
#'   In the current implementation, \code{qtuv_m_s.e.irrad()}, accepts multiple
#'   values as arguments for only one parameter at a time. In the case of
#'   elevation, both `ground.elevation` and `measurement.elevation` can have
#'   each one or more values. When too many values are passed in the call, only
#'   the first one is used.
#' 
#' @note The Quick TUV calculator has multiple output modes that return
#'   different types of computed values. The use of output mode 5 is hard-coded
#'   in this function as other modes return summary values rather than spectral
#'   data. Package 'foqat' provides a more flexible alternative supporting other
#'   output modes in addition to mode 5.
#'   
#'   If the argument passed to \code{w.length} is a numeric vector, the range
#'   and length are used to reconstruct the accepted parameters. The returned
#'   spectrum has always a uniformly spaced wavelengths.
#'   
#'   When using this function, more detailed metadata are available than when
#'   reading an output file, as not all the simulation input parameters are
#'   listed in the output text.
#'   
#'   In interactive use of the Quick TUV Calculator, the same parameters as
#'   accepted by \code{qtuv_s.e.irradiance()} as arguments are
#'   entered through the web interface at 
#'   \url{https://www.acom.ucar.edu/Models/TUV/Interactive_TUV/}. This page 
#'   together with its documentation, can be consulted for additional
#'   information on the parameters and the model. 
#'   
#' @section Warning!: This function connects to a server managed by UCAR, the
#'   University Corporation for Atmospheric Research located in the U.S.A. to
#'   obtain simulated spectral data. UCAR manages the U.S. National Science
#'   Foundation National Center for Atmospheric Research (NSF NCAR) on behalf of
#'   NSF. As any download with the HTTPS protocol, using this function entails
#'   some risk. To minimize the risk, the returned page is saved as plain text,
#'   checked for conformity with the expected content, and if valid decoded into
#'   an R data object. When using the default argument \code{file = NULL}, the
#'   file used is a temporary one and is deleted before the function returns the
#'   call, irrespective of it being conformant or not.
#' 
#'   The administrators of the Quick TUV Calculator at UCAR suggest a maximum
#'   load of approximately 100 spectral simulations per day and user. For larger
#'   workloads they encourage the local installation of the TUV model which is
#'   open-source and freely available. A local installation, also allows access
#'   to the full set of input parameters and outputs. Currently a local instance
#'   of the TUV model can be called from R with package 'foqat'.
#'       
#' @references
#' Sasha Madronich (2017-2021) Tropospheric Ultraviolet and Visible radiation (TUV) 
#' model. \url{https://www2.acom.ucar.edu/modeling/tropospheric-ultraviolet-and-visible-tuv-radiation-model}.
#' Visited on 2024-08-29.
#' 
#' @export
#' 
qtuv_s.e.irrad <- 
  function(w.length = list(wStart = 280, 
                           wStop = 420, 
                           wIntervals = 140), 
           sun.elevation = NULL,
           geocode = data.frame(lon = 0, 
                                lat = 51.5, 
                                address = "Greenwich"), 
           time = lubridate::now(),
           tz = NULL,
           locale = readr::default_locale(),
           ozone.du = 300, 
           albedo = 0.1, 
           ground.altitude = 0, 
           measurement.altitude = NULL, 
           clouds = data.frame(optical.depth = 0.00, 
                               base = 4.00, 
                               top = 5.00),
           aerosols = data.frame(optical.depth = 0.235,
                                 ssaaer = 0.990, 
                                 alpha = 1.000),
           num.streams = 2,
           spectra = list(direct = 1.0, 
                          diffuse.down = 1.0, 
                          diffuse.up = 0),
           added.vars = NULL,
           label = "",
           server.url = "https://www.acom.ucar.edu/cgi-bin/acom/TUV/V5.3/tuv",
           file = NULL) {
    # check parameters
    num.streams <- as.integer(abs(num.streams))
    if (num.streams == 2L) {
      num.streams <- -2L
    }
    stopifnot("'num.streams' must be 2 or 4" = num.streams %in% c(-2L, 4L))
    
    if (is.null(measurement.altitude) || measurement.altitude < ground.altitude) {
      measurement.altitude <- ground.altitude
    }
    if (any(ground.altitude > 10)) {
      stop("Altitude must be expressed in km!")
    }
    # Select inputMode
    if (is.null(sun.elevation)) {
      # use geocode and time
      inputMode <- 0L
      zenith <- 0 # ignored
    } else {
      inputMode <- 1L
      stopifnot(sun.elevation >= -20 && sun.elevation <= 90)
      zenith <- sun.elevation - 90
    }
    
    # always use the same outputMode
    outputMode <- 5L # spectral irradiance 
    if (is.null(tz)) {
      tz <- locale[["tz"]]
    }
    # timeStamp uses the time zone of the time passed by the user
    tz <- lubridate::tz(time)
    timeStamp <- strftime(time, 
                          format = "%H:%M:%S", 
                          tz = tz)
    # computations are always based on "UTC" (="GMT")
    time <- lubridate::with_tz(time, tzone = "UTC")
    date <- strftime(time, 
                     format = "%Y%m%d", 
                     tz = "UTC")
    time.h <- lubridate::hour(time) + 
      lubridate::minute(time) / 60 + 
      lubridate::second(time) /3600
    
    if (is.numeric(w.length)) {
      # numeric vectors need to be converted to parameter values
      wIntervals = length(w.length) - 1
      wl.step <- (max(w.length) - min(w.length)) / wIntervals
      wStart = min(w.length) - wl.step / 2 
      wStop = max(w.length) + wl.step / 2 
      wIntervals = length(w.length) - 1
    } else if (is.list(w.length) && 
      # parameters supplied by user
               all(c("wStart", "wStop", "wIntervals") %in% 
                   names(w.length))) {
      wStart <- w.length$wStart
      wStop <- w.length$wStop
      wIntervals <- w.length$wIntervals
    } else {
      stop("Invalid defintion of wavelengths in 'w.length'")
    }
    # build URL to call Quick TUV
    url <- paste0(c(server.url, 
                    "?wStart=", wStart, 
                    "&wStop=", wStop, 
                    "&wIntervals=", wIntervals, 
                    "&inputMode=", inputMode, 
                    "&latitude=", geocode$lat, 
                    "&longitude=", geocode$lon, 
                    "&date=", date, 
                    "&timeStamp=", timeStamp,  
                    "&zenith=", zenith, 
                    "&ozone=", ozone.du, 
                    "&albedo=", albedo, 
                    "&gAltitude=", ground.altitude, 
                    "&mAltitude=", measurement.altitude, 
                    "&taucld=", clouds$optical.depth, 
                    "&zbase=", clouds$base, 
                    "&ztop=", clouds$top, 
                    "&tauaer=", aerosols$optical.depth, 
                    "&ssaaer=", aerosols$ssaaer, 
                    "&alpha=", aerosols$alpha, 
                    "&time=", time.h, 
                    "&outputMode=", 5, 
                    "&nStreams=", num.streams, 
                    "&dirsun=", spectra$direct, 
                    "&difdn=", spectra$diffuse.down, 
                    "&difup=", spectra$diffuse.up), 
                  collapse='')
    
    if (is.null(file)) {
      # NOTE: I need to check if textConnection could be used instead
      qtuv.file <- tempfile()
      on.exit(unlink(qtuv.file))
    } else {
      if (grepl("\\.txt$", file)) {
        qtuv.file <- file
      } else {
        qtuv.file <- paste(file, ".txt", sep = "")
      }
    }
    utils::download.file(url, 
                         destfile = qtuv.file, 
                         quiet = !getOption("photobiology.verbose", 
                                            default = FALSE))
    if (file.exists(qtuv.file)) {
      z <- read_qtuv_txt(file = qtuv.file,
                         ozone.du = ozone.du,
                         label = label,
                         tz = "UTC",
                         locale = locale,
                         added.vars = added.vars)
      where_measured(z) <- geocode
      when_measured(z) <- time
      how_measured(z) <- paste(how_measured(z), "Using", num.streams, "streams.")
    } else {
      warning("No file was returned by the server")
      z <- source_spct()
    }
    attr(z, "qtuv.url") <- url
    z
  }

#' @rdname qtuv_s.e.irrad
#' 
#' @export
#' 
qtuv_m_s.e.irrad <- 
  function(w.length = list(wStart=280, 
                           wStop=420, 
                           wIntervals=140), 
           sun.elevation = NULL,
           geocode = data.frame(lon = 0, 
                                lat = 51.5, 
                                address = "Greenwich"), 
           time = lubridate::now(),
           tz = NULL,
           locale = readr::default_locale(),
           ozone.du = 300, 
           albedo = 0.1, 
           ground.altitude = 0, 
           measurement.altitude = NULL, 
           clouds = data.frame(optical.depth = 0.00, 
                               base = 4.00, 
                               top = 5.00),
           aerosols = data.frame(optical.depth = 0.235,
                                 ssaaer = 0.990, 
                                 alpha = 1.000),
           num.streams = 2,
           spectra = list(direct = 1.0, 
                          diffuse.down = 1.0, 
                          diffuse.up = 0),
           added.vars = NULL,
           label = "",
           server.url = "https://www.acom.ucar.edu/cgi-bin/acom/TUV/V5.3/tuv",
           file = NULL) {
    # check parameters
    num.streams <- as.integer(abs(num.streams))
    if (num.streams == 2L) {
      num.streams <- -2L
    }
    stopifnot("'num.streams' must be 2 or 4" = num.streams %in% c(-2L, 4L))
    
    # ensure consistent length of altitudes
    if (is.null(measurement.altitude)) {
      measurement.altitude <- ground.altitude
    } else if (length(measurement.altitude) > length(ground.altitude)) {
      ground.altitude <- rep_len(ground.altitude, length.out = length(measurement.altitude))
    } else if (length(ground.altitude) > length(measurement.altitude)) {
      measurement.altitude <- rep_len(measurement.altitude, length.out = length(ground.altitude))
    }
    altitudes <- data.frame(measurement = measurement.altitude,
                            ground = ground.altitude)
    
    if (length(ozone.du) > 1L) {
      if (length(ozone.du) > 25) {
        message("Please, do not overload the Quick TUV calculator")
      }
      z <- list()
      for (ozone in ozone.du) {
        member.name <- paste("ozone", ozone, sep = ".")
        if (!is.null(file) && is.character(file)) {
          file.name <- paste(gsub("\\.txt$", "", file), 
                             "-", member.name, ".txt", sep = "")
        } else {
          file.name <- NULL
        }
        z[[paste(member.name)]] <-
          qtuv_s.e.irrad(w.length = w.length, 
                         sun.elevation = sun.elevation[[1]],
                         geocode = geocode[1, ], 
                         time = time[[1]],
                         tz = tz,
                         locale = locale,
                         ozone.du = ozone, 
                         albedo = albedo, 
                         ground.altitude = altitudes[1, "ground"], 
                         measurement.altitude = altitudes[1, "measurement"], 
                         clouds = clouds[1, ],
                         aerosols = aerosols[1, ],
                         num.streams = num.streams,
                         spectra = spectra,
                         added.vars = added.vars,
                         label = label,
                         file = file.name)
      }
    } else if (length(sun.elevation) > 1L) {
      if (length(sun.elevation) > 25) {
        message("Please, do not overload the Quick TUV calculator")
      }
      z <- list()
      for (one.elevation in sun.elevation) {
        member.name <- paste("sun.elevation", one.elevation, sep = ".")
        if (!is.null(file) && is.character(file)) {
          file.name <- paste(gsub("\\.txt$", "", file), 
                             "-", member.name, ".txt", sep = "")
        } else {
          file.name <- NULL
        }
        z[[paste(member.name)]] <-
          qtuv_s.e.irrad(w.length = w.length, 
                         sun.elevation = one.elevation,
                         geocode = geocode[1, ], 
                         time = time[[1]],
                         tz = tz,
                         locale = locale,
                         ozone.du = ozone.du[[1]], 
                         albedo = albedo, 
                         ground.altitude = altitudes[1, "ground"], 
                         measurement.altitude = altitudes[1, "measurement"], 
                         clouds = clouds[1, ],
                         aerosols = aerosols[1, ],
                         num.streams = num.streams,
                         spectra = spectra,
                         added.vars = added.vars,
                         label = label,
                         file = file.name)
      }
    } else if (length(time) > 1L) {
      if (length(time) > 25) {
        message("Please, do not overload the Quick TUV calculator")
      }
      z <- list()
      # use index as for converts POSIXct into numeric
      for (i in seq_along(time)) {
        member.name <- paste("time", format(time[[i]]), sep = ".")
        if (!is.null(file) && is.character(file)) {
          file.name <- paste(gsub("\\.txt$", "", file), 
                             "-", make.names(member.name), ".txt", sep = "")
        } else {
          file.name <- NULL
        }
        z[[paste(member.name)]] <-
          qtuv_s.e.irrad(w.length = w.length, 
                         sun.elevation = sun.elevation[[1]],
                         geocode = geocode[1, ], 
                         time = time[[i]],
                         tz = tz,
                         locale = locale,
                         ozone.du = ozone.du[[1]], 
                         albedo = albedo, 
                         ground.altitude = altitudes[1, "ground"], 
                         measurement.altitude = altitudes[1, "measurement"], 
                         clouds = clouds[1, ],
                         aerosols = aerosols[1, ],
                         num.streams = num.streams,
                         spectra = spectra,
                         added.vars = added.vars,
                         label = label,
                         file = file.name)
      }
    } else if (nrow(geocode) > 1L) {
      if (nrow(geocode) > 25) {
        message("Please, do not overload the Quick TUV calculator")
      }
      if (!"address" %in% colnames(geocode) || 
          any(is.na(geocode[["address"]]))) {
        geocode[["address"]] <- 
          paste("geo", 1L:nrow(geocode), sep = ".")
      } else {
        geocode[["address"]] <- make.unique(geocode[["address"]])
      }
      z <- list()
      # use index to walk through data frame rows
      for (i in seq_along(geocode[[1]])) {
        member.name <- paste("geocode", geocode[i, "address"], sep = ".")
        if (!is.null(file) && is.character(file)) {
          file.name <- paste(gsub("\\.txt$", "", file), 
                             "-", make.names(member.name), ".txt", sep = "")
        } else {
          file.name <- NULL
        }
        z[[paste(member.name)]] <-
          qtuv_s.e.irrad(w.length = w.length, 
                         sun.elevation = sun.elevation[[1]],
                         geocode = geocode[i, ], 
                         time = time[[1]],
                         tz = tz,
                         locale = locale,
                         ozone.du = ozone.du[[1]], 
                         albedo = albedo, 
                         ground.altitude = altitudes[1, "ground"], 
                         measurement.altitude = altitudes[1, "measurement"], 
                         clouds = clouds[1, ],
                         aerosols = aerosols[1, ],
                         num.streams = num.streams,
                         spectra = spectra,
                         added.vars = added.vars,
                         label = label,
                         file = file.name)
      }
    } else if (nrow(clouds) > 1L) {
      if (nrow(clouds) > 25) {
        message("Please, do not overload the Quick TUV calculator")
      }
      if (!"label" %in% colnames(clouds) || 
          any(is.na(clouds[["label"]]))) {
        clouds[["label"]] <- 
          paste("geo", 1L:nrow(clouds), sep = ".")
      } else {
        clouds[["label"]] <- make.unique(clouds[["label"]])
      }
      z <- list()
      # use index to walk through data frame rows
      for (i in seq_along(clouds[[1]])) {
        member.name <- paste("clouds", clouds[i, "label"], sep = ".")
        if (!is.null(file) && is.character(file)) {
          file.name <- paste(gsub("\\.txt$", "", file), 
                             "-", make.names(member.name), ".txt", sep = "")
        } else {
          file.name <- NULL
        }
        z[[paste(member.name)]] <-
          qtuv_s.e.irrad(w.length = w.length, 
                         sun.elevation = sun.elevation[[1]],
                         geocode = geocode[1, ], 
                         time = time[[1]],
                         tz = tz,
                         locale = locale,
                         ozone.du = ozone.du[[1]], 
                         albedo = albedo, 
                         ground.altitude = altitudes[1, "ground"], 
                         measurement.altitude = altitudes[1, "measurement"], 
                         clouds = clouds[i, ],
                         aerosols = aerosols[1, ],
                         num.streams = num.streams,
                         spectra = spectra,
                         added.vars = added.vars,
                         label = label,
                         file = file.name)
      }
    } else if (nrow(aerosols) > 1L) {
      if (nrow(aerosols) > 25) {
        message("Please, do not overload the Quick TUV calculator")
      }
      if (!"label" %in% colnames(aerosols) || 
          any(is.na(aerosols[["label"]]))) {
        aerosols[["label"]] <- 
          paste("label", 1L:nrow(aerosols), sep = ".")
      } else {
        aerosols[["label"]] <- make.unique(aerosols[["label"]])
      }
      z <- list()
      # use index to walk through data frame rows
      for (i in seq_along(aerosols[[1]])) {
        member.name <- paste("aerosols", aerosols[i, "label"], sep = ".")
        if (!is.null(file) && is.character(file)) {
          file.name <- paste(gsub("\\.txt$", "", file), 
                             "-", make.names(member.name), ".txt", sep = "")
        } else {
          file.name <- NULL
        }
        z[[paste(member.name)]] <-
          qtuv_s.e.irrad(w.length = w.length, 
                         sun.elevation = sun.elevation[[1]],
                         geocode = geocode[1, ], 
                         time = time[[1]],
                         tz = tz,
                         locale = locale,
                         ozone.du = ozone.du[[1]], 
                         albedo = albedo, 
                         ground.altitude = altitudes[1, "ground"], 
                         measurement.altitude = altitudes[1, "measurement"], 
                         clouds = clouds[1, ],
                         aerosols = aerosols[i, ],
                         num.streams = num.streams,
                         spectra = spectra,
                         added.vars = added.vars,
                         label = label,
                         file = file.name)
      }
    } else if (nrow(altitudes) > 1L) {
      if (nrow(altitudes) > 25) {
        message("Please, do not overload the Quick TUV calculator")
      }
      if (!"label" %in% colnames(altitudes) || 
          any(is.na(altitudes[["label"]]))) {
        altitudes[["label"]] <- 
          paste("label", 1L:nrow(altitudes), sep = ".")
      } else {
        altitudes[["label"]] <- make.unique(altitudes[["label"]])
      }
      z <- list()
      # use index to walk through data frame rows
      for (i in seq_along(altitudes[[1]])) {
        member.name <- paste("altitudes", altitudes[i, "label"], sep = ".")
        if (!is.null(file) && is.character(file)) {
          file.name <- paste(gsub("\\.txt$", "", file), 
                             "-", make.names(member.name), ".txt", sep = "")
        } else {
          file.name <- NULL
        }
        z[[paste(member.name)]] <-
          qtuv_s.e.irrad(w.length = w.length, 
                         sun.elevation = sun.elevation[[1]],
                         geocode = geocode[1, ], 
                         time = time[[1]],
                         tz = tz,
                         locale = locale,
                         ozone.du = ozone.du[[1]], 
                         albedo = albedo, 
                         ground.altitude = altitudes[i, "ground"], 
                         measurement.altitude = altitudes[i, "measurement"], 
                         clouds = clouds[1, ],
                         aerosols = aerosols[1, ],
                         num.streams = num.streams,
                         spectra = spectra,
                         added.vars = added.vars,
                         label = label,
                         file = file.name)
      }
    }else {
      z <- list(
        qtuv_s.e.irrad(w.length = w.length, 
                       sun.elevation = sun.elevation[[1]],
                       geocode = geocode[1, ], 
                       time = time[[1]],
                       tz = tz,
                       locale = locale,
                       ozone.du = ozone.du[[1]], 
                       albedo = albedo, 
                       ground.altitude = ground.altitude[[1]], 
                       measurement.altitude = measurement.altitude[[1]], 
                       clouds = clouds[1, ],
                       aerosols = aerosols[1, ],
                       num.streams = num.streams,
                       spectra = spectra,
                       added.vars = added.vars,
                       label = label,
                       file = file)
      )
    }
    source_mspct(z)
  }

#' Clouds descriptor
#' 
#' Constructor of a named list of parameter values to be used as argument to
#' parameter \code{clouds} of function \code{\link{qtuv_s.e.irrad}()}.
#' 
#' @param cloud.type character One of "clear.sky", "cirrus", "stratocumulus" or
#'   "overcast". 
#'   
#' @details This function provide a rough approximation for parameter values.
#'   In reality there is large variation in the cloud optical depths (COD) and
#'   in the elevation at which clouds are located, within each type of cloud.
#'   The TUV model assumes a continuous uniform cloud layer, thus the normally
#'   discontinuous cover of cumulus clouds cannot be described.
#'   
#' @return A one-row data frame with members "optical.depth", "base", "top" and
#'   "label".
#'   
#' @export
#' 
#' @examples
#' 
#' qtuv_clouds("clear.sky")
#' qtuv_clouds("cirrus")
#' qtuv_clouds(c("clear.sky", "cirrus"))
#' 
qtuv_clouds <- function(cloud.type = "clear.sky") {
  z <- data.frame()
  for (cld in cloud.type) {
    cld.df <- 
      switch(cld,
             clear.sky = data.frame(optical.depth = 0.00, 
                                    base = 4.00, 
                                    top = 5.00,
                                    label = "clear.sky"),
             cirrus = data.frame(optical.depth = 5.00, 
                                 base = 6.00, 
                                 top = 8.00,
                                 label = "cirrus"),
             stratocumulus = data.frame(optical.depth = 10.00, 
                                        base = 4.00, 
                                        top = 5.00,
                                        label = "stratocumulus"),
             overcast = data.frame(optical.depth = 20.00, 
                                   base = 0.50, 
                                   top = 2.00,
                                   label = "overcast"),
             data.frame(optical.density = 0.00, 
                        base = 4.00, 
                        top = 5.00,
                        label = "clear.sky")
      )
    
    z <- rbind(z, cld.df)
   }
  z
}
