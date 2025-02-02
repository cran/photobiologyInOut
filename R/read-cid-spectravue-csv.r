#' Read File Saved by CID's SpectraVue.
#'
#' Read wavelength and spectral data from Measurements.CSV files exported from
#' CID Bio-Sciences' \strong{SpectraVue CI-710s} (not the older CI-710!) leaf
#' spectrometer, importing them into R. Available metadata is also extracted
#' from the file. \code{read_cid_spectravue_csv()} only accepts "row oriented"
#' CSV files. These may contain multiple spectra, one per row.
#'
#' @param file character
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date and time are extracted
#'   from the file.
#' @param geocode A data frame with columns \code{lon} and \code{lat} used to
#'   set attribute \code{"where.measured"}. If \code{NULL}, the default, the
#'   geocode is extracted from the file, if present, and if \code{NA} the
#'   "where.measured" attribute is not set.
#' @param label character string. If \code{NULL}, the default, the value of the
#'   "tag" present in the file is used, and if \code{NA} the "what.measured"
#'   attribute is not set.
#' @param tz character Time zone is by default that of the machine's locale.
#' @param locale The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#' @param range	numeric A vector of length two, or any other object for which
#'   function \code{range()} will return range of wavelengths expressed in
#'   nanometres.
#' @param simplify logical If TRUE, single spectra are returned as individual
#'   spectra instead of collections of length one.
#' @param absorbance.to character Affects only absorbance measurements:
#'   \code{"object"}, \code{"all"} or \code{"A"}.
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error (\code{TRUE}) instead of a warning (\code{FALSE}), or the test
#'   is disabled (\code{NA}).
#' @param ... additional arguments passed to the constructor of the
#'   `filter_spct` object.
#'
#' @details SpectraVue's row-wise spectral \code{Measurements.CSV} files contain
#'   columns with metadata on the right edge, followed by columns with data for
#'   each of the 2048 pixels or wavelengths. The value in column "Mode"
#'   indicates the quantity measured, decoded into variables \code{Tpc},
#'   \code{Rpc} or \code{A}. The data the rows in the CSV file are read and
#'   stored in \code{filter_spct},  \code{reflector_spct} or \code{object_spct}
#'   objects. These objects are collected into a single \code{filter_mspct},
#'   \code{reflector_mspct}, \code{object_mspct} or \code{generic_spct} object
#'   and returned.
#'
#'   Spectral data outside the range 400 nm to 1000 nm are very noisy and thus
#'   outside the valid range for the measurements. Out-of-range spectral data
#'   can also be cause by calibration drift. Consequently, reading of data is
#'   done always with the range check disabled, while whether a check is used
#'   before returning the collection of spectra depends on the argument passed
#'   to \code{strict.range} which by default is set to disable checks. This is
#'   done, because in most cases measurements from this instrument tend to
#'   require further processing before they comply with theoretical expectations
#'   of \eqn{Tfr + Rfr + Afr = 1}.
#'
#' @section Internal \emph{vs.} total transmittance and absorbance: Spectravue
#'   returns transmittance values labelled with \code{Transmittance} as
#'   \code{Mode}. Transmittance values are not \emph{total} as most of the
#'   scattered light transmitted is not detected. Absorbance (\code{Abs}) values
#'   returned labelled with \code{Absorbance} as \code{Mode} are for absorbance
#'   computed from the Transmittance. This estimate of absorbance overestimates
#'   real absorbance in the case of scattering materials like plant leaves. It
#'   is best to save spectral data acquired as absorbance into objects of class
#'   `object_spct` containing reflectance and transmittance, the default, as
#'   this preserves all the acquired data.
#'
#' @section Specular \emph{vs.} total reflectance: This function assumes that
#'   SpectraVue returns close to \emph{total} reflectance readings. Given the optics
#'   of the instrument this is likely only an approximation.
#'
#' @note SpectraVue creates three or four \code{.CSV} files for each series of
#'   measurements saved. Of these files, this function reads the one with name
#'   ending in \code{Measurements.CSV}. The first part of the file name gives
#'   the time of the session, but as the files can contain multiple spectra
#'   measured at different times, the time metadata is extracted separately for
#'   each spectrum. We provide a default argument for \code{range} that discards
#'   data for short and long wavelengths because values outside this range are
#'   according to the instrument's manual outside the usable range and in
#'   practice extremely noisy.
#'
#' @section Warning!!:
#' CID's support has answered on 2022-05-19 that the extremely biased (plainly
#' wrong!) values of transmittance measured by this instrument are due to a
#' design flaw and that they are working on a solution for the problem. In
#' practice, reflectance seems biased but usable as an instrument-specific
#' quantity with arbitrary units. Transmittance and absorbance seem useless as
#' values are wrong by about an order of magnitude.
#' 
#' @return An object of class \code{filter_spct}, \code{relector_spct},
#'   \code{object_spct} or \code{generic_mspct}.
#'
#' @export
#' @references \emph{https://cid-inc.com/}
#'
#' @examples
#'
#' # read file containing a single reflectaance spectrum
#'
#'  file.name <-
#'    system.file("extdata", "cid-spectravue-Rpc-Measurements.csv",
#'                package = "photobiologyInOut", mustWork = TRUE)
#'
#'  cid_filter.spct <-
#'    read_cid_spectravue_csv(file = file.name)
#'  summary(cid_filter.spct)
#'
#'  cid_filter.spct <-
#'    read_cid_spectravue_csv(file = file.name, simplify = FALSE)
#'  summary(cid_filter.spct)
#'
#'  # read data measured as absorbance (A, Rpc and Tpc)
#'
#'  file.name <-
#'    system.file("extdata", "cid-spectravue-Abs-Measurements.csv",
#'                package = "photobiologyInOut", mustWork = TRUE)
#'  cid.object_spct <-
#'    read_cid_spectravue_csv(file = file.name)
#'  summary(cid.object_spct)
#'
#'  cid_A.filter_spct <-
#'    read_cid_spectravue_csv(file = file.name, absorbance.to = "A")
#'  summary(cid_A.filter_spct)
#' 
read_cid_spectravue_csv <- 
  function(file,
           date = NULL,
           geocode = NULL,
           label = NULL,
           tz = NULL,
           locale = readr::default_locale(),
           range = c(380, 1100),
           simplify = TRUE,
           absorbance.to = "object",
           strict.range = NA,
           ...) {
  if (!grepl("Measurements.csv$", file, ignore.case = TRUE)) {
    warning("Only processed measurements CSV files from CID's SpectraVue CI-710s",
            "leaf spectrometer are supported.")
    return(photobiology::filter_mspct())
  }
  
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
  if (is.character(date)) {
    date <- anytime::anytime(date, tz = tz)
    message("Dates and times from file overriden by user!")
  }
  
  # read whole file
  data <- utils::read.csv(file = file, header = FALSE, as.is = TRUE)
  
  # extract metadata and do conversions
  headers <- tibble::as_tibble(data[-1, 1:9])
  colnames(headers) <- c("Date", "Mode", "integ.time", 
                         "boxcar", "num.scans", "Tag",
                         "Version", "CalibrationID", "Coordinates") 
  
  headers[["integ.time"]] <- as.double(headers[["integ.time"]])
  headers[["boxcar"]] <- as.integer(headers[["boxcar"]])
  headers[["num.scans"]] <- as.integer(headers[["num.scans"]])
  headers[["var.name"]] <- c(Transmittance = "Tpc",
                             Reflectance = "Rpc",
                             Absorbance = "A",
                             "Abs Reflectance" = "Rpc",
                             "Abs Transmittance" = "Tpc")[headers[["Mode"]]]
  
  headers[["when.measured"]] <- 
    if(is.null(date)) {
      lubridate::mdy_hms(headers[["Date"]], tz = tz)
    } else {
      date
    }
  
  # I do not have information about the format used for Coordinates
  headers[["where.measured"]] <-
    if (is.null(geocode) || is.na(geocode)) {
      SunCalcMeeus::na_geocode()
    } else {
      geocode
    }
  
  # fill missing Tags
  selector <- is.na(headers[["Tag"]])
  headers[["Tag"]][selector] <- "NN"
  headers[["Tag"]] <- paste(headers[["Mode"]], headers[["Tag"]], sep = ".")
  headers[["Tag"]] <- make.names(headers[["Tag"]], unique = TRUE)
    
  # extract wavelengths and spectral data
  # file rows have a trailing "," and leftmost columns with metadata
  data <- t(unname(data[ , -c(1:10, 2059)])) # data becomes a matrix
  
  w.length <- data[ , 1] 
  data <- data[ , -1, drop = FALSE]

  zz <- photobiology::generic_mspct()
  # we use a while loop as the number of spectra (= rows) per iteration varies
  col <- 1L
  while (col <= ncol(data)) {
    if (absorbance.to == "object" && headers[["var.name"]][col] == "A") {
      # A is always followed by Rpc and Tpc
      z <- data.frame(w.length, data[ , col + 1], data[ , col + 2]) # data is a matrix!
      names(z) <- c("w.length", headers[["var.name"]][col + 1:2])

      old.opts <- options("photobiology.strict.range" = NA_integer_)
      z <- photobiology::as.object_spct(z, 
                                        Tfr.type = "total", 
                                        Rfr.type = "total", 
                                        ...)
      spct.name <- headers[["Tag"]][col]
      col <- col + 2
    } else {
      z <- list(w.length, data[ , col]) # data is a matrix!
      names(z) <- c("w.length", headers[["var.name"]][col])
      
      old.opts <- options("photobiology.strict.range" = NA_integer_)
      if (headers[["var.name"]][col] == "A") {
        z <- photobiology::as.filter_spct(z, 
                                          Tfr.type = "internal", 
                                          strict.range = NA_integer_,
                                          ...)
        if (absorbance.to != "all") {
          col <- col + 2
        } 
      } else if (headers[["var.name"]][col] %in% c("Afr", "Tfr", "Apc", "Tpc")) {
        z <- photobiology::as.filter_spct(z, 
                                          Tfr.type = "total", 
                                          strict.range = NA_integer_,
                                          ...)
      } else if (headers[["var.name"]][col] %in% c("Rfr", "Rpc")) {
        z <- photobiology::as.reflector_spct(z, 
                                             Rfr.type = "total",
                                             strict.range = NA_integer_,
                                             ...)
      } else {
        stop("Assertion failed for column named: ", headers[["var.name"]][col])
      }
      spct.name <- headers[["Tag"]][col]
    }
    z <- photobiology::clip_wl(x = z, range = range)
    options(old.opts)
    
    comment(z) <-
      paste(paste("CID Bio-Science SpectraVue CI-710s file: '", basename(file), 
                  "' with Mode '", headers[["Mode"]], "'\nImported on ", 
                  lubridate::now(tzone = "UTC"), " UTC\n", sep = ""),
            paste(headers[col, ], collapse = "\n"), sep = "\n")
    
    photobiology::setWhenMeasured(z, headers[col, "when.measured"])
    photobiology::setWhereMeasured(z, SunCalcMeeus::na_geocode())
    photobiology::setWhatMeasured(z, paste(label, headers[["Tag"]][col]))
    photobiology::setHowMeasured(z, "CID Bio-Science SpectraVue CI-710s leaf spectrometer")
    
    instr.descriptor <- list(time = headers[["when.measured"]][col],
                             spectrometer.name = "CI-710s",
                             spectrometer.sn = NA_character_,
                             bench.grating = "default",
                             bench.filter = "default",
                             bench.slit = "default",
                             max.counts = 65535,
                             n.pixels = 2048)
    
    photobiology::setInstrDesc(z, instr.descriptor)
    
    instr.settings <- list(calibration.id =  headers[["CalibrationID"]][col],
                           pix.selector = TRUE,
                           integ.time = headers[["integ.time"]][col] * 1e3, # ms -> us
                           num.scans = headers[["num.scans"]][col],
                           tot.time = headers[["integ.time"]][col] * 
                             headers[["num.scans"]][col] * 1e3, # ms -> us
                           rel.signal = NA_real_,
                           boxcar.width = headers[["boxcar"]][col],
                           linearized = TRUE,
                           dark.subtracted = TRUE)
    
    photobiology::setInstrSettings(z, instr.settings)
    
    zz[[spct.name]] <- z
    
    col <- col + 1
  }
  
  if (simplify && length(zz) == 1) {
      zz[[1]]
  } else {
    switch(photobiology::shared_member_class(zz)[[1]],
           object_spct = photobiology::as.object_mspct(zz, strict.range = strict.range),
           filter_spct = photobiology::as.filter_mspct(zz, strict.range = strict.range),
           reflector_spct = photobiology::as.reflector_mspct(zz, strict.range = strict.range),
           zz)
  }
}
