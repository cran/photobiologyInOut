#' Read multiple foreign files with spectral data
#' 
#' Read spectra from a homogeneous list of files based on a path and a list of
#' filenames or a path and a search pattern for files. The imported spectra are
#' returned as a single object of one of the collection of spectra classes from
#' package 'photobiology'.
#' 
#' @details This function iterates over a list of file names reading them with
#'   the function passed as argument to `.fun` and combines the spectra as a
#'   collection of spectra of a class suitable for the spectral objects returned
#'   by the argument to `.fun`. This function can either return for each file
#'   read either a single spectrum as an object of class `generic_spct` or a
#'   class derived from it, or a collection of spectra of class `generic_mspct`
#'   or a class derived from it. The class of the returned object depends on
#'   the class of the member spectra.
#' 
#' @param path character A path point to the location of the files.
#' @param list character A vector or list of character strings pointing to
#'    files relative to \code{path},
#' @param pattern character A search pattern to select files within \code{path}.
#'    See \code{\link{list.files}} which is used internally. Argument ignored
#'    is list is non-null.
#' @param .fun function One of the functions exported by this package for
#'    reading spectral data.
#' @param ... Named arguments passed ot the call to \code{.fun}.
#' 
#' @return An object of class `generic_mspct` or a class derived from it, 
#'   containing a collection of member spectra of class `generic_spct` or of 
#'   one of the classes derived from it.
#' 
#' @export
#' 
#' 
read_foreign2mspct <- function(path = ".", 
                               list = NULL, 
                               pattern = NULL, 
                               .fun, ...) {
  if (is.null(list)) {
    list <- list.files(path = path, pattern = pattern)
  }
  tags <- make.names(basename(path = list), unique = TRUE)
  names(list) <- tags 
  zz <- generic_mspct()
  for (t in tags) {
    z <- .fun(file = list[[t]], ...)
    if (is.generic_spct(z)) {
      # append one spectrum
      zz[[t]] <- z
    } else if (is.generic_mspct(z)) {
      # concatenate collections of spectra
      zz <- c(zz, z)
    } else {
      warning("Skipping object returned by '.fun' for ", list[[t]], 
              "as it is neither a spectrum nor a collection of spectra.")
    }
  }
  switch(photobiology::shared_member_class(zz)[[1]],
         raw_spct = photobiology::as.raw_mspct(zz),
         cps_spct = photobiology::as.cps_mspct(zz),
         source_spct = photobiology::as.source_mspct(zz),
         object_spct = photobiology::as.object_mspct(zz),
         filter_spct = photobiology::as.filter_mspct(zz),
         reflector_spct = photobiology::as.reflector_mspct(zz),
         zz)
}