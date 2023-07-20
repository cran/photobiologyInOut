#' Read File Saved by Wasatch's Enlighten.
#' 
#' Read wavelength and spectral data from the data section of a file as
#' output by Enlighten importing them into R. Parse the header of a file to
#' extract the acquisition time, instrument name and serial number, as well
#' additional metadata related to the instrument and its settings. Function
#' \code{read_wasatch_csv()} only accepts "column oriented" CSV files.
#' 
#' @param file character 
#' @param s.qty character, possibly named. The name of the quantity using the
#'   conventions accepted used in package 'photobiology' that is to be imported
#'   from column "Processed" from the file.
#' @param extra.cols character What to do non-processed data columns if present 
#'   in file. One of \code{"keep"}, \code{"drop.pixel"}, \code{"drop"} or 
#'   \code{"split"}.
#' @param ... additional arguments passed to the constructor of the spectrum
#'   object.
#' @param scale.factor numeric vector of length 1, or length equal to the number
#'   of rows (= detector pixels). Numeric multiplier applied to returned
#'   spectral values.
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date and time are extracted from the
#'   file header.
#' @param geocode A data frame with columns \code{lon} and \code{lat} used to
#'   set attribute \code{"where.measured"}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone is by default that of the machine's locale.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#' @param simplify logical If TRUE, single spectra are returned as individual 
#'   spectra instead of collections of length one. 
#'   
#' @details Enlighten's column-wise CSV files contain at least two columns,
#'   Wavelength and Processed. In the header the Technique used is recorded.
#'   Additional data columns can be present. Column Pixel contains the pixel
#'   index in the array as integers. Columns Raw, Dark and Reference contain
#'   detector counts data. Technique is used to guess the type of spectrum
#'   stored in the column named Processed, which can be detector counts or
#'   derived values. By default the data are read into a single spectrum object
#'   and all columns retained, but only the data in Processed are interpreted as
#'   spectral data corresponding to the class of the object. If passed
#'   \code{extra.cols = "drop"}, only Wavelength and Processed are copied to the
#'   returned object, while if passed \code{extra.cols = "drop.pixel"} only the
#'   contents of column Pixel are discarded. If passed \code{extra.cols =
#'   "split"} all columns containing spectral data are each read into a separate
#'   spectrum, these are collected and a "generic_mspct" object containing them
#'   returned. \code{extra.cols} can be a named vector of mappings, of length at
#'   least one but possibly longer. If longer a "generic_mspct" is returned,
#'   otherwise a spectrum object as inferred from the name each column is mapped
#'   to.
#'   
#' @note Enlighten, the free software from Wasatch Photonics can save spectra in
#'   a variety of additional formats: different types of CSV files, plain text
#'   and JSON. Plain text files contain no metadata or even column headers and
#'   if the need arises can be read with R function \code{read.table()}. JSON
#'   files contain the most detailed metadata.
#' 
#' @return An object of a class derived from \code{generic_spct} such as
#'   \code{raw_spct} or \code{filter_spct}. \code{generic_spct} is derived from
#'   tibble and data frame.
#'   
#' @export
#' @references \url{https://wasatchphotonics.com/}
#'   \url{https://wasatchphotonics.com/product-category/software/}
#' 
#' @section Acknowledgements:
#' We thank Ruud Niesen from Photon Mission 
#' (\url{https://photonmission.com/}) for organizing the loan of the
#' spectrometer used to produce the various files needed for the development of 
#' this function.
#' 
#' @examples
#' 
#'  file.name <- 
#'    system.file("extdata", "enlighten-wasatch-scope.csv",
#'                package = "photobiologyInOut", mustWork = TRUE)
#'               
#'  wasatch.raw.spct <- 
#'    read_wasatch_csv(file = file.name)
#'  summary(wasatch.raw.spct)
#' 
#'  wasatch.raw.spct <- 
#'    read_wasatch_csv(file = file.name, s.qty = "counts")
#'  summary(wasatch.raw.spct)
#' 
#'  wasatch.raw.spct <- 
#'    read_wasatch_csv(file = file.name, s.qty = c(Processed = "counts"))
#'  summary(wasatch.raw.spct)
#' 
#'  wasatch.raw.spct <- 
#'    read_wasatch_csv(file = file.name, extra.cols = "drop")
#'  summary(wasatch.raw.spct)
#' 
read_wasatch_csv <- function(file,
                             date = NULL,
                             geocode = NULL,
                             label = NULL,
                             tz = NULL,
                             locale = readr::default_locale(),
                             s.qty = NULL,
                             extra.cols = "keep",
                             scale.factor = 1,
                             simplify = TRUE,
                             ...) {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
  if (is.character(date)) {
    date <- anytime::anytime(date, tz = tz)
  }
  
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  line01 <- scan(file = file, nlines =  1, skip = 0,
                 what="character", quiet = TRUE)
  if (line01[1] != "ENLIGHTEN") {
    warning("Input file was not created by ENLIGHTEN as expected: skipping")
    return(photobiology::generic_spct())
  }
  file_header <- scan(file = file, nlines = 40, skip = 0,
                      what="character", sep = "\n", quiet = TRUE)
  NonASCII <- tools::showNonASCII(file_header)
  if (length(NonASCII) > 0L) {
    warning("Found non-ASCII characters in file header: ", 
            NonASCII,
            "replacing with '.'.")
    file_header <- iconv(file_header, to = "ASCII", sub = " ")
  }
  
  npixels <-
    as.integer(sub("Pixel Count,", "", 
                   file_header[grepl("Pixel Count,", file_header)], 
                   fixed = TRUE))
  
  if (is.null(date)) {
    time.stamp <- sub("Timestamp,", "", 
                  file_header[grepl("Timestamp,", file_header)])
    if (is.null(tz)) {
      tz <- tz(lubridate::now())
      warning("Using computer's current TZ setting: ", tz)
    }
    date <- lubridate::ymd_hms(time.stamp, tz = tz)
  }
  
  wp.technique <-
    sub("Technique,", "", 
        file_header[grepl("Technique,", file_header)], 
        fixed = TRUE)

  measurement.id <-
    sub("Measurement ID,", "", 
        file_header[grepl("Measurement ID,", file_header)], 
        fixed = TRUE)
  
  integ.time <- 
    as.numeric(sub("Integration Time,", "", 
                   file_header[grepl("Integration Time,", file_header)], 
                   fixed = TRUE)) # ms
  
  num.scans <- 
    as.numeric(sub("Scan Averaging,", "", 
                   file_header[grepl("Scan Averaging,", file_header)], 
                   fixed = TRUE)) # number of scans?
  
  boxcar <-
    as.numeric(sub("Boxcar,", "", 
                   file_header[grepl("Boxcar,", file_header)], 
                   fixed = TRUE)) # ???
  
  sensor.temperature <-
    as.numeric(sub("Temperature,", "", 
                   file_header[grepl("^Temperature,", file_header)], 
                   fixed = TRUE))
  
  laser.temperature <-
    as.numeric(sub("Laser Temperature,", "", 
                   file_header[grepl("Laser Temperature,", file_header)], 
                   fixed = TRUE))

  inst.sn <-
    sub("Serial Number,", "", 
        file_header[grepl("Serial Number,", file_header)], 
        fixed = TRUE)
  
  inst.model <-
    sub("Model,", "", 
        file_header[grepl("Model,", file_header)], 
        fixed = TRUE)
  
  inst.slit <-
    sub("Slit Width,", "", 
        file_header[grepl("Slit Width,", file_header)], 
        fixed = TRUE) # um
  
  dark.subtracted <-
    grepl("Dark substracted", 
          file_header[grepl("Note,", file_header)], 
          fixed = TRUE) # T/F
  
  # search for data header
  i <- 30L
  while(!grepl("wavelength", tolower(file_header[i]))) {
    i <- i + 1L
  }
  
  # inspect column names
  column.names <- file_header[i]
  column.names <- scan(text = column.names, what = "", 
                       sep = ",", quiet = TRUE)

  x <- utils::read.csv(file = file,
                       header = TRUE,
                       skip = i,
                       nrows = npixels)
  
  # As a default we try to guess the quantity from the technique used
  if (is.null(s.qty)) {
    data.cols <- 
      switch(wp.technique,
             "Scope"               = c(Processed = "counts"),
             "Relative Irradiance" = c(Processed = "counts"),
             "Absorbance"          = c(Processed = "A"),
             "Transmission"        = c(Processed = "Tpc"),
             "Raman"               = c(Processed = "counts")
    )
  } else {
    # Support user-supplied mappings
    data.cols <- s.qty
    #  but also a plain variable name with "Processed" as implied target
    if (is.null(names(data.cols))) {
      names(data.cols) <- "Processed"
    }
  }
  
  if (extra.cols == "split" && length(data.cols) == 1L) {
    other.cols <- 
        c(Reference = "counts",
          Raw = "counts",
          Dark = "counts")[c( "Reference", "Raw","Dark") %in% names(x)]

    data.cols <- c(data.cols, other.cols)
    extra.cols <- "drop"
  }
  
  zz <- photobiology::generic_mspct()
  for (col in names(data.cols)) {
    z <- x
    processed.name <- data.cols[col]
    
    # we rename only columns we recognize
    column.names.map <- c(Wavelength = "w.length", data.cols[col]) 
    
    #  stopifnot(all(names(column.names.map) %in% names(z))) # assertion check
    
    selector <- colnames(z) %in% names(column.names.map)
    colnames(z)[selector] <- column.names.map[colnames(z)[selector]]
    
    z[[processed.name]] <- z[[processed.name]] * scale.factor
    
    old.opts <- options("photobiology.strict.range" = NA_integer_)
    if (processed.name %in% c("s.e.irrad", "s.q.irrad")) {
      z <- photobiology::as.source_spct(z, time.unit = "second", ...)
    } else if (processed.name == "counts") {
      z <- photobiology::as.raw_spct(z, ...)
    } else if (processed.name == "cps") {
      z <- photobiology::as.cps_spct(z, ...)
    } else if (processed.name %in% c("Afr", "Tfr", "Apc", "Tpc", "A")) {
      z <- photobiology::as.filter_spct(z, ...)
    } else if (processed.name %in% c("Rfr", "Rpc")) {
      z <- photobiology::as.reflector_spct(z, ...)
    } else {
      z <- photobiology::as.generic_spct(z, ...)
    }
    options(old.opts)
    
    if (extra.cols %in% c("drop", "drop.all")) {
      z <- photobiology::drop_user_cols(z)
    } else if (extra.cols == "drop.pixel") {
      z <- photobiology::drop_user_cols(z, 
                                        keep.also = setdiff(colnames(z), "Pixel"))
    } else if (extra.cols != "keep") {
      warning("Bad 'extra.cols' argument '", extra.cols, "', assuming 'keep'.")
    }
    
    comment(z) <-
      paste(paste("Wasatch Enlighten file: '", basename(file), 
                  "' of technique '", wp.technique, "'\nImported on ", 
                  lubridate::now(tzone = "UTC"), " UTC\n", sep = ""),
            paste(file_header[1:(i-1)], collapse = "\n"), sep = "\n")
    
    photobiology::setWhenMeasured(z, date)
    photobiology::setWhereMeasured(z, geocode)
    photobiology::setWhatMeasured(z, label)
    
    instr.descriptor <- list(time = date,
                             spectrometer.name = inst.model,
                             spectrometer.sn = inst.sn,
                             bench.grating = "default",
                             bench.filter = NA_character_,
                             bench.slit = inst.slit,
                             max.counts = 64000)
    
    photobiology::setInstrDesc(z, instr.descriptor)
    
    instr.settings <- list(wp.technique = wp.technique,
                           measurement.id = measurement.id,
                           pix.selector = TRUE,
                           integ.time = integ.time * 1e3, # ms -> us
                           num.scans = num.scans,
                           tot.time = integ.time * num.scans * 1e3, # ms -> us
                           rel.signal = NA_real_,
                           boxcar.width = boxcar,
                           linearized = TRUE,
                           dark.subtracted = dark.subtracted,
                           sensor.temperature = sensor.temperature,
                           laser.temperature = laser.temperature)
    
    photobiology::setInstrSettings(z, instr.settings)
    
    attr(z, "file.header") <- file_header[1:(i-1)]
    
    zz[[col]] <- z
    
  }
  
  if (simplify && length(zz) == 1) {
    zz[[1]]
  } else {
    switch(photobiology::shared_member_class(zz)[[1]],
           raw_spct = photobiology::as.raw_mspct(zz),
           cps_spct = photobiology::as.cps_mspct(zz),
           source_spct = photobiology::as.source_mspct(zz),
           object_spct = photobiology::as.object_mspct(zz),
           filter_spct = photobiology::as.filter_mspct(zz),
           reflector_spct = photobiology::as.reflector_mspct(zz),
           zz)
  }
}
