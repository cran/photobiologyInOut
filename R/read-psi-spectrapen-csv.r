#' Read '.CSV' File Saved by PSI's Software.
#'
#' Reads and parses the header of a processed .CSV file as output by the by the
#' PSI (Photon Systems Instruments, Czech Republic) SpectraPen miniature
#' spectrometer.
#'
#' @param file character string.
#' @param start.row integer The first line to read, counting from the top of the
#'   file.
#' @param date a \code{POSIXct} object to use to set the \code{"when.measured"}
#'   attribute. If \code{NULL}, the default, the date is extracted from the file
#'   header.
#' @param geocode A data frame with columns \code{lon} and \code{lat} used to
#'   set attribute \code{"where.measured"}. If \code{NULL}, the default, the
#'   location is extracted from the file header.
#' @param label character string, but if \code{NULL} the label from \code{file}
#'   header is used, if the label is missing, the index is used, and if
#'   \code{NA} the "what.measured" attribute is not set.
#' @param tz character Time zone used for interpreting times saved in the file
#'   header.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use
#'   \code{\link[readr]{locale}} to create your own locale that controls things
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#'
#' @return A \code{source_mspct} object.
#' @export
#' @references \url{https://psi.cz/}
#'
#' @examples
#'  # fetch path to example file to read
#'  file.name <-
#'    system.file("extdata", "spectrum-psi-spectrapen-SP.csv",
#'                package = "photobiologyInOut", mustWork = TRUE)
#'
#'  spectrapen.mspct <- read_spectrapen_csv(file = file.name)
#'
#'  spectrapen.mspct
#'  getWhenMeasured(spectrapen.mspct)
#'  getWhatMeasured(spectrapen.mspct)
#'  cat(comment(spectrapen.mspct))
#' 
read_spectrapen_csv <- function(file,
                                start.row = 1,
                                date = NULL,
                                geocode = NULL,
                                label = NULL,
                                tz = NULL,
                                locale = readr::default_locale()) {
  skip.idx <- start.row - 1
  if (is.null(tz)) {
    tz <- locale$tz
  }
 
  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  file_header <- scan(file = file, nlines = 1, skip = skip.idx,
                      what = "character", quiet = TRUE)
  stopifnot(file_header[1] == "Irradiance")
  if (grepl("[\u03BCW/cm2/nm]", file_header[2])) { # greek mu
    qty.in <- "energy"
    spct.data.var <- "s.e.irrad"
    scale.factor <- 1e-2
  } else if (grepl("[\u03BCE/m2/s/nm]", file_header[2])) { # greek mu
    qty.in <- "photon"
    spct.data.var <- "s.q.irrad"
    scale.factor <- 1e-6
  }
  repeat {
    spct_times <- scan(file = file, nlines = 1, skip = skip.idx + 1, sep = ",",
                       what = "character", quiet = TRUE)
    if (length(spct_times) == 0L) {
      skip.idx <- skip.idx + 1L
    } else {
      break()
    }
  }
  spct_times <- scan(file = file, nlines = 1, skip = skip.idx + 1, sep = ",",
                     what = "character", quiet = TRUE)
  stopifnot(spct_times[1] == "Time")
  spct_times <- spct_times[c(-1, -length(spct_times))] # comma at end of line adds a field
  spct_times <- lubridate::dmy_hms(spct_times, tz = tz)
  spct_idxs <- scan(file = file, nlines = 1, skip = skip.idx + 2, sep = ",",
                    what = "character", quiet = TRUE)
  stopifnot(spct_idxs[1] == "Index")
  spct_idxs <- spct_idxs[c(-1, -length(spct_idxs))] # comma at end of line adds a field
  spct_names <- scan(file = file, nlines = 1, skip = skip.idx + 3, sep = ",",
                     what = "character", quiet = TRUE)
  stopifnot(spct_names[1] == "Name")
  spct_names <- spct_names[c(-1, -length(spct_names))] # comma at end of line adds a field
  to_rename <- spct_names == ""
  spct_names[to_rename] <- paste("spct", spct_idxs[to_rename], sep = ".")
  spct_geocodes <- scan(file = file, nlines = 1, skip = skip.idx + 4, sep = ",",
                        what = "character", quiet = TRUE)
  stopifnot(spct_geocodes[1] == "GPS")
  spct_geocodes <- spct_geocodes[c(-1, -length(spct_geocodes))] # comma at end of line adds a field
  skip.idx <- skip.idx + 5
  repeat {
    marker <- scan(file = file, nlines = 1, skip = skip.idx, sep = ",",
                   what = "character", quiet = TRUE)[1]
    skip.idx <- skip.idx + 1
    if (marker == "[nm]") break()
  }
  
  spct.df <- utils::read.csv(file, header = FALSE, 
                             col.names = c("w.length", spct_names),
                             skip = skip.idx, 
                             nrows = 256)
  
  z.mspct <- split2mspct(spct.df, 
                         member.class = "source_spct",
                         w.length.var = "w.length",
                         spct.data.var = spct.data.var)
  
  comment.text <- paste("PSI SpectraPen file '", basename(file), "' imported on ", 
                        lubridate::now(tzone = "UTC"), " UTC", sep = "")
  
  for (i in seq_along(z.mspct)) {
    tmp.spct <- z.mspct[[i]]
    tmp.spct[[spct.data.var]] <- tmp.spct[[spct.data.var]] * scale.factor
    when_measured(tmp.spct) <- spct_times[i]
    what_measured(tmp.spct) <- spct_names[i]
    where_measured(tmp.spct) <- na_geocode()
    how_measured(tmp.spct) <- "PSI SpectraPen handheld array spectrometer"
    comment(tmp.spct) <- comment.text
    z.mspct[[i]] <- tmp.spct
  }
  
  z.mspct
}
