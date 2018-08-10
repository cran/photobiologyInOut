#' Read '.CSV' FReD database.
#'
#' Reads a CSV data file downloaded from the FReD (Floral Reflectance Database)
#' and then imports wavelengths and spectral reflectance values and flower ID.
#'
#' @param file character string
#' @param date a \code{POSIXct} object, but if \code{NULL} the date stored in
#'   file is used, and if \code{NA} no date variable is added
#' @param geocode A data frame with columns \code{lon} and \code{lat}.
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
#' @return A reflectance_spct object.
#' @export
#' @references \url{http://www.reflectance.co.uk}
#' Arnold SEJ, Faruq S, Savolainen V, McOwan PW, Chittka L, 2010 FReD: The 
#' Floral Reflectance Database — A Web Portal for Analyses of Flower Colour. 
#' PLoS ONE 5(12): e14287. doi:10.1371/journal.pone.0014287
#' 
#' @keywords misc
#'
read_FReD_csv <- function(file,
                           date = NA,
                           geocode = NULL,
                           label = NULL,
                           tz = NULL,
                           locale = readr::default_locale()) {
  if (is.null(tz)) {
    tz <- locale$tz
  }
  
  label <- paste("File:", basename(file), label)
  
  z <- readr::read_csv(file = file, 
                      col_names = c("flower.id", "w.length", "Rfr"),
                      col_types = readr::cols(),
                      locale = locale,
                      skip = 0)

  z <- photobiology::as.reflector_spct(z, Rfr.type = "total")

  if (length(unique(z[["flower.id"]])) > 1) {
    warning("Spectrum contains data for multiple flower IDs")
  }
  
  comment(z) <-
    paste(paste("FReD file '", basename(file), "' imported on ", 
                lubridate::now(tzone = "UTC"), " UTC", sep = ""),
          sep = "\n")
  photobiology::setWhenMeasured(z, date)
  photobiology::setWhereMeasured(z, geocode)
  photobiology::setWhatMeasured(z, label)
  z
}
