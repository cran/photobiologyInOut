#' Color reproduction index
#'
#' Wrapper on function \code{\link[colorSpec]{computeCRI}} from package
#' 'colorSpec' that accepts \code{\link[photobiology]{source_spct}} objects.
#'
#' @param spct source_spct A single light source spectrum.
#' @param adapt logical If TRUE, then a special chromatic adaption is performed,
#'   see Details in \code{\link[colorSpec]{computeCRI}}.
#' @param attach logical If TRUE, then a large list of intermediate calculations
#'   is attached to the returned number, as attribute data. This attached list
#'   includes data for all special 14 color samples, although the last 6 do not
#'   affect the returned CRI
#' @param tol numeric For the CRI to be meaningful the chromaticities of the
#'   test and reference illuminants must be sufficiently close in the CIE 1960
#'   uniform chromaticity space. If the tolerance is exceeded, the function
#'   returns NA. The default tol=5.4e-3 is the one recommended by the CIE, but
#'   the argument allows the user to override it.
#' @param named logical Whether to set the name attribute of the returned
#'   value to the name of the spectrum passed as argument if possible.
#'
#' @details Please see \code{\link[colorSpec]{computeCRI}} for the details of
#'   the computations and references.
#'   
#' @return A numeric value between zero and 100, or NA if the light is not
#'   white enough.
#'
#' @examples
#'
#' spct_CRI(white_led.source_spct)
#' spct_CRI(sun.spct)
#'
#' @export
#' 
spct_CRI <- function(spct,  
                     adapt = TRUE, 
                     attach = FALSE, 
                     tol = 5.4e-3,
                     named = FALSE) {
  if (!requireNamespace("spacesXYZ", quietly = TRUE)) { 
    warning("Package 'spacesXYZ' needs to be installed to compute CRI.")
    return(NA_real_)
  }
  spct.name <- substitute(spct)
  if (is.name(spct.name)) {
    spct.name <- paste(as.character(spct.name), "_", sep = "")
  } else {
    spct.name <- "NN.spct_"
  }
  if (photobiology::is.source_spct(spct) && 
      photobiology::getMultipleWl(spct) == 1) {
    x <- spct2colorSpec(spct)
  } else {
    warning("Spectrum in an object of class source_spct required!")
    return(NA_real_)
  }
  z <- colorSpec::computeCRI(x = x,
                             adapt = adapt,
                             attach = attach,
                             tol = tol)
  if (!named) {
    unname(z)
  } else {
    names(z) <- gsub("spct_", spct.name, names(z))
    z
  }
}

#' Correlated color temperature
#'
#' Wrapper on function \code{\link[colorSpec]{computeCCT}} from package
#' 'colorSpec' that accepts \code{\link[photobiology]{source_spct}} objects.
#'
#' @param spct source_spct A single light source spectrum.
#' @param isotherms character A vector whose elements match one of the available
#'   isotherm families: 'robertson', 'mccamy', and 'native'. Matching is partial
#'   and case-insensitive. When more than one family is given, a matrix is
#'   returned, see Value. When isotherms = 'native' the isotherms are defined
#'   implicitly as lines perpendicular to the locus, see Details in
#'   \code{\link[spacesXYZ]{CCTfromXYZ}}. The character NA
#'   (\code{NA_character_}) is taken as a synonym for 'native'.
#' @param locus character Valid values are 'robertson' and 'precision', see
#'   above. Matching is partial and case-insensitive.
#' @param strict logical The CIE considers the CCT of a chromaticity uv to be
#'   meaningful only if the distance from uv to the Planckian locus is less than
#'   or equal to 0.05 (in CIE UCS 1960). If strict=FALSE, then this condition is
#'   ignored. Otherwise, the distance is computed along the corresponding
#'   isotherm, and if it exceeds 0.05 the returned CCT is set to NA.
#' @param named logical Whether to set the name attribute of the returned
#'   value to the name of the spectrum passed as argument if possible.
#'
#' @details Please see \code{\link[colorSpec]{computeCCT}} for the details of
#'   the computations and references.
#'   
#' @return A numeric value for "color temperature " in degrees Kelvin. 
#'
#' @examples
#'
#' spct_CCT(white_led.source_spct)
#' spct_CCT(sun.spct)
#'
#' @export
#' 
spct_CCT <- function(spct, 
                     isotherms = 'robertson', 
                     locus = 'robertson', 
                     strict = FALSE,
                     named = FALSE) {
  if (!requireNamespace("spacesXYZ", quietly = TRUE)) { 
    warning("Package 'spacesXYZ' needs to be installed to compute CCT.")
    return(NA_real_)
  }
  spct.name <- substitute(spct)
  if (is.name(spct.name)) {
    spct.name <- paste(as.character(spct.name), "_", sep = "")
  } else {
    spct.name <- "NN.spct_"
  }
  if (photobiology::is.source_spct(spct) && 
      photobiology::getMultipleWl(spct) == 1) {
    x <- spct2colorSpec(spct)
  } else {
    warning("Spectrum in an object of class source_spct required!")
    return(NA_real_)
  }
  z <- colorSpec::computeCCT(x = x,
                             isotherms = isotherms,
                             locus = locus,
                             strict = strict)
  if (!named) {
    unname(z)
  } else {
    names(z) <- gsub("spct_", spct.name, names(z))
    z
  }
}

#' Spectral (color) similarity index
#'
#' Wrapper on function \code{\link[colorSpec]{computeSSI}} from package
#' 'colorSpec' that accepts \code{\link[photobiology]{source_spct}} objects.
#'
#' @param spct,reference.spct source_spct Single light source spectra.
#' @param digits integer The number of digits after the decimal point in the
#'   returned vector. According to Holm the output should be rounded to the
#'   nearest integer, which corresponds to \code{digits = 0}. To return full
#'   precision, set \code{digits = Inf}.
#' @param isotherms character This is only used when reference=NULL. It is
#'   passed to \code{\link[colorSpec]{computeCCT}} in order to compute the CCT
#'   of each test spectrum.
#' @param locus character This is only used when reference=NULL. It is
#'   passed to \code{\link[colorSpec]{computeCCT}} in order to compute the CCT
#'   of each test spectrum.
#' @param named logical Whether to set the name attribute of the returned
#'   value to the name of the spectrum passed as argument if possible.
#'
#' @details Please see \code{\link[colorSpec]{computeSSI}} for the details of
#'   the computations and references.
#'
#' @return A numeric value between zero and 100.
#' 
#' @export
#' 
#' @examples
#'
#' spct_SSI(white_led.source_spct, sun.spct)
#' 
spct_SSI <- function(spct,
                     reference.spct = NULL, 
                     digits = 0, 
                     isotherms = 'mccamy', 
                     locus = 'robertson',
                     named = FALSE) {
  if (!requireNamespace("spacesXYZ", quietly = TRUE)) { 
    warning("Package 'spacesXYZ' needs to be installed to compute SSI.")
    return(NA_real_)
  }
  spct.name <- substitute(spct)
  if (is.name(spct.name)) {
    spct.name <- paste(as.character(spct.name), "_", sep = "")
  } else {
    spct.name <- "NN.spct_"
  }
  if (photobiology::is.source_spct(spct) && 
      photobiology::getMultipleWl(spct) == 1) {
    x <- spct2colorSpec(spct)
  } else {
    warning("Spectrum in an object of class source_spct required!")
    return(NA_real_)
  }
  if (is.null(reference.spct) || nrow(reference.spct) == 0) {
    z <- colorSpec::computeSSI(x = spct2colorSpec(spct),
                               digits = digits,
                               isotherms = isotherms,
                               locus = locus)
  } else {
    if (photobiology::is.source_spct(reference.spct)) {
      reference <- spct2colorSpec(reference.spct)
    } else {
      warning("Reference spectrum in an object of class source_spct required!")
      return(NA_real_)
    }
    z <- colorSpec::computeSSI(x = x,
                               reference = reference,
                               digits = digits,
                               isotherms = isotherms,
                               locus = locus)
  }
  if (!named) {
    unname(z)
  } else {
    names(z) <- gsub("spct_", spct.name, names(z))
    z
  }
}

