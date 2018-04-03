#' @section Warning!: Most of the file formats supported are not standardized,
#'   and are a moving target because of changes in instrument firmware and
#'   support software. In addition the output format, especially with models,
#'   can depend on settings that users can alter. So do check that import is
#'   working as expected, and if not, please please raise a bug or support
#'   ticket and upload one example of an incorrectly decoded file.
#' 
#' @note From version 0.4.4 the time zone (tz) used for decoding dates
#'   and times in files imported defaults to "UTC". In most cases you will need
#'   to pass the tz (or the locale) where the file was created as an argument to
#'   the functions!
#'   
#' @references 
#' Aphalo, Pedro J. (2015) The r4photobiology suite. UV4Plants Bulletin, 2015:1,
#' 21-29. \url{https://doi.org/10.19232/uv4pb.2015.1.14}.
#'
#' @import photobiology
#' 
"_PACKAGE"
