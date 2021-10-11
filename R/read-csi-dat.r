#' Read '.DAT' file(s) saved by modern Campbell Scientific loggers.
#' 
#' Reads and parses the header of a processed data file as output by the PC400
#' or PC200W programmes extracting variable names, units and quantities from the
#' header. Uses the comment attribute to store the metadata.
#' 
#' @param file Path to file as a character string.
#' @param geocode A data frame with columns \code{lon} and \code{lat} used to
#'   set attribute \code{"where.measured"}.
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
#' @param na character Vector of strings to interpret as missing values. Set 
#'   this option to character() to indicate no missing values.
#' @param ... Further named arguments currently passed to \code{read_csv()}.
#'   
#' @return \code{read_csi_dat()} returns a \code{tibble::tibble} object.
#' @export
#' @references \url{https://www.campbellsci.eu/}
#' 
#' @note This function is not useful for .DAT and .PRN files from old CSI
#'   loggers and software. Those were simple files, lacking metadata, which was
#'   stored in separate .FLD files.
#'   
read_csi_dat <- function(file, 
                         geocode = NULL, 
                         label = NULL,
                         data_skip = 0, 
                         n_max = Inf, 
                         locale = readr::default_locale(),
                         na = c("", "NA", "NAN"),
                         ...) {

  label.file <- paste("File: ", basename(file), sep = "")
  if (is.null(label)) {
    label <- label.file
  } else if (!is.na(label)) {
    label <- paste(label.file, label, sep = "\n")
  }
  
  file_header <- readr::read_lines(file, n_max = 4, skip = 0)
  head_line <-  scan(text = file_header[1L], what = character(), 
                     sep = ",", quiet = TRUE)
  head_line = paste(head_line, collapse = "\t")
  col_names <- scan(text = file_header[2L], what = character(), 
                    sep = ",", quiet = TRUE)
  units <- scan(text = file_header[3L], what = character(), 
                sep = ",", encoding = "UTF-8", quiet = TRUE)
  units <- sub("^$", "NA", units)
  qty <- scan(file, 
              what = character(), nlines = 1, skip = 3,
              sep = ",", quiet = TRUE)
  metadata <- paste(col_names, units, qty, sep = "\t", collapse = "\n")
  comment.txt <- paste(label, head_line, metadata, sep = "\n---\n")
  # readr::read_csv automatically recognizes dates as well as numeric
  # values allowing this very simple and flexible implementation.
  z <-
    readr::read_csv(file = file, 
                    skip = 4 + data_skip, 
                    n_max = n_max, 
                    col_names = col_names,
                    col_types = readr::cols(),
                    progress = FALSE,
                    ...)

  attr(z, "file.header") <- file_header
  comment(z) <- comment.txt
  z
}
