#' Read '.DAT' file(s) saved by modern Campbell Scientific loggers.
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
#' @return \code{read_csi_dat()} returns a \code{tibble::tibble} object.
#' @export
#' @references \url{http://www.r4photobiology.info} \url{http://www.campbellsci.eu/}
#' 
#' @note This function is not useful for .DAT and .PRN files from old CSI
#'   loggers and software. Those were simple files, lacking metadata, which was
#'   stored in separate .FLD files.
#'   
read_csi_dat <- function(file, geocode = NULL, label = NULL,
                         data_skip = 0, n_max = Inf, 
                         locale = readr::default_locale()) {
  label <- paste("File:", basename(file), label)

  head_line <-  scan(file, 
                     what = character(), nlines = 1, skip = 0,
                     sep = ",")
  head_line = paste(head_line, collapse = "\t")
  col_names <- scan(file, 
                    what = character(), nlines = 1, skip = 1,
                    sep = ",")
  units <- scan(file, 
                what = character(), nlines = 1, skip = 2,
                sep = ",", encoding = "UTF-8")
  units <- sub("^$", "NA", units)
  qty <- scan(file, 
              what = character(), nlines = 1, skip = 3,
              sep = ",")
  metadata <- paste(col_names, units, qty, sep = "\t", collapse = "\n")
  comment.txt <- paste(label, head_line, metadata, sep = "\n---\n")
  # readr::read_csv automatically recognizes dates as well as numeric
  # values allowing this very simple and flexible implementation.
  z <-
    readr::read_csv(file, skip = 4 + data_skip, n_max = n_max,
             col_types = "", col_names = col_names)
  comment(z) <- comment.txt
  z
}
