#' Read '.CSV' file(s) dowloaded from YoctoPuce modules.
#' 
#' Reads and parses the header of a processed data file as output by the PC400
#' or PC200W programmes extracting variable names, units and quantities from the
#' header. Uses the comment attribute to store the metadata.
#' 
#' @param file Path to file as a character string.
#' @param geocode A data frame with columns \code{lon} and \code{lat}.
#' @param label character string, but if \code{NULL} the value of \code{file} is
#'   used, and if \code{NA} the "what.measured" attribute is not set.
#' @param data_skip integer Number of records (rows) to skip from the actual
#'   data block.
#' @param n_max integer Maximum number of records to read.
#' @param locale	The locale controls defaults that vary from place to place. The
#'   default locale is US-centric (like R), but you can use 
#'   \code{\link[readr]{locale}} to create your own locale that controls things 
#'   like the default time zone, encoding, decimal mark, big mark, and day/month
#'   names.
#'   
#' @return \code{read_yoctopuce_csv()} returns a \code{tibble::tibble} object.
#' @export
#' @references \url{http://www.r4photobiology.info} \url{http://www.yoctopuce.com/}
#' 
#' @note This function should be able to read data log files from any YoctoPuce
#'   USB interface module with data logging capabilities.
#'   
read_yoctopuce_csv <- function(file,
                               geocode = NULL,
                               label = NULL,
                               data_skip = 0, 
                               n_max = Inf, 
                               locale = readr::default_locale()) {
  data_skip <- as.integer(data_skip)
  
  label <- paste("File:", basename(file), label)
  
  z <- readr::read_delim(file = file,file,
                  delim = ";",
                  col_types = readr::cols(), 
                  n_max = n_max + data_skip)
  
  names(z) <- make.names(names(z), unique = TRUE)
  
  if ("UNIX.time" %in% names(z)) {
    z <- dplyr::select(z, -"UNIX.time")
  }
  
  if (data_skip >= 1) {
    z <- dplyr::slice(z, -(1:data_skip))
  }

  comment(z) <- paste("Data from a YoctoPuce module.", label)
  z
}
