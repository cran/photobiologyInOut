# Some "export" functions are not defined as S3 methods for 'hyperSpec'. For
# 'colorSpec' as.colorSpec methods are now defined as the maintainer has updated
# colorSpec after reading this note (thanks Glenn!). Defining these methods as
# S3 generics in this independently developed package would be looking for
# future trouble. The "import" methods are now defined also as "as...." S3
# methods. As a single class "colorSpec" corresponds to several different
# classes from package "photobiology", a function is defined which automatically
# detects the suitable output class. In contrast the S3 methods trigger an error
# if there is a mismatch.
# 

# R matrix ----------------------------------------------------------------
#  Moved to 'photobiology' and expanded adding coercion methods.

# hyperSpec ---------------------------------------------------------------

#' Convert 'hyperSpec::hyperSpec' objects
#' 
#' Convert \code{hyperSpec::hyperSpec} objects containing VIS and UV radiation
#' data into spectral objects (xxxx_spct, xxxx_mspct) as defined in package 
#' 'photobiology' and vice versa, preserving as much information as possible. As
#' \code{hyperSpec} can contain other kinds of spectral data, it does make sense
#' to use these functions only with objects containing data that can be handled
#' by both packages.
#' 
#' @note Objects of class \code{hyperSpec::hyperSpec} contain metadata or class 
#'   data from which the quantity measured and the units of expression can be 
#'   obtained. However, units as included in the objects are not well 
#'   documented making automatic conversion difficult. When using this function
#'   the user may need to use parameter \code{multiplier} to scale the data to
#'   what is expected by the object constructors defined in package
#'   'photobiology' and use parameter \code{spct.data.var} to select the
#'   quantity.
#'   
#'   \code{hyperSpec::hyperSpec} objects may use memory more efficiently
#'   than spectral objects of the classes for collections of spectra defined in
#'   package 'photobiology' as wavelengths are assumed to be the same for all
#'   member spectra, and stored only once while this assumption is not made for
#'   collections of spectra, allowing different wavelengths and lengths for the
#'   component spectra. Wavelengths are stored for each spectrum, but as
#'   spectral classes are derived from 'tbl_df' in many cases no redundant
#'   copies of wavelength data will be made in memory in spite of the more
#'   flexible semantics of the objects.
#'   
#' @section Warning!: Always check the sanity of the imported or exported data
#'   values, as guessing is needed when matching the different classes, and the
#'   functions defined here are NOT guaranteed to return valid data wihtout help
#'   from the user through optional function arguments.
#' 
#' @param x hyperSpec object
#' @param member.class character One of the spectrum classes defined in package
#'   'photobiology'.
#' @param spct.data.var character The name to be used for the 'spc' data when
#'   constructing the spectral objects.
#' @param multiplier numeric A multiplier to be applied to the 'spc' data to do
#'   unit or scale conversion. For example "a.u." units in some examples in
#'   package 'hyperSpec' seem to have scale factors applied.
#' @param ... currently ignored.
#' 
#' @export
#' 
#' @examples 
#' # example run only if 'hyperSpec' is available
#' if (requireNamespace("hyperSpec", quietly = TRUE)) {
#'   library(hyperSpec)
#'   data(laser)
#'   wl(laser) <- 
#'     list(wl = 1e7 / (1/405e-7 - wl (laser)),
#'          label = expression (lambda / nm))
#'   laser.mspct <- hyperSpec2mspct(laser, "source_spct", "s.e.irrad")
#'   class(laser.mspct)
#' }
#' 
hyperSpec2mspct <- function(x, 
                            member.class, 
                            spct.data.var,
                            multiplier = 1,
                            ...) {
  if (requireNamespace("hyperSpec", quietly = TRUE)) {
    stopifnot(inherits(x, "hyperSpec"))
    # spc data (spectra) are stored as rows in a matrix, consequently
    # we transpose the matrix so that each spectrum is in a column
    y <- cbind(hyperSpec::wl(x), t(x$spc) * multiplier) 
    colnames(y) <- c("w.length", paste("spc", 1:nrow(x), sep = ""))
    y <- tibble::as_tibble(y)
    z <- photobiology::split2mspct(x = y, 
                                   member.class = member.class, 
                                   spct.data.var = spct.data.var)
    other.vars <- setdiff(colnames(x), "spc")
    for (r in seq_along(z)) { # nrow(x) is same as length(z)
      for (var in other.vars) {
        z[[r]][[var]] <- rep(x@data[[var]][r], hyperSpec::nwl(x))
      }
    }
    comment(z) <- paste('Converted from "hyperSpec" object\n',
                        'dim: ', 
                        paste(names(dim(x)), dim(x), collapse = " "),
                        '\ncolnames: ', paste(colnames(x), collapse = ", "),
                        sep = "")
  } else {
    warning("Package 'hyperSpec' needs to be installed.")
    z <- call(sub("spct", "mscpt", member.class))
  }
  z
}

#' @rdname hyperSpec2mspct
#'
#' @export
#'
hyperSpec2spct <- function(x, multiplier = 1, ...) {
  y <- hyperSpec2mspct(x, multiplier = multiplier, ...)
  z <- photobiology::rbindspct(y)
  comment(z) <- comment(x[[1]])
  z
}

#' @rdname hyperSpec2mspct
#' 
#' @export
#' 
mspct2hyperSpec <- function(x, 
                            spct.data.var,
                            multiplier = 1,
                            ...) {
  if (requireNamespace("hyperSpec", quietly = TRUE)) {
    stopifnot(is.any_mspct(x))
    spct.names <- names(x)
    spct.selector <- rep(TRUE, length(x))
    for (i in seq_along(x)) {
      temp <- x[[i]]
      s.column <- temp[[spct.data.var]] * multiplier
      wl.current <- temp[["w.length"]]
      if (i == 1L) {
        mat <- s.column
        wl.prev <- wl.current
      } else {
        if (!all(wl.current == wl.prev)) {
          spct.selector[i] <- FALSE
          next()
        }
        mat <- rbind(mat, s.column) # as row!
      }
    }
    methods::new("hyperSpec", 
                 spc = mat, 
                 wavelength = wl.prev, 
                 data = data.frame(spct_name = factor(spct.names[spct.selector])), 
                 labels = list(spct = expression("spct_name"), 
                               spc = expression("I / a.u."),
                               .wavelength = expression("lambda/nm")))
  } else {
    warning("Package 'hyperSpec' needs to be installed.")
    NA
  }
}

#' @rdname hyperSpec2mspct
#' 
#' @export
#' 
spct2hyperSpec <- function(x, 
                           spct.data.var = NULL,
                           multiplier = 1,
                           ...) {
  y <- photobiology::generic_mspct(list(x), class = class(x)[1])
  mspct2hyperSpec(y,
                  spct.data.var,
                  multiplier)
}

# pavo --------------------------------------------------------------------

#' Convert "pavo::rspec" objects
#' 
#' Convert between 'pavo::rspec' objects containing spectral reflectance data 
#' into spectral objects (xxxx_spct, xxxx_mspct) as defined in package
#' 'photobiology'.
#' 
#' @note Objects of class \code{pavo::rspec} do not contain metadata or class 
#'   data from which the quantity measured and the units of expression could be 
#'   obtained. When using this function the user needs to use parameter 
#'   \code{multiplier} to convert the data to what is expected by the object 
#'   constructors defined in package 'photobiology' and use parameter 
#'   \code{spct.data.var} to select the quantity.
#'   
#'   \code{pavo::rspec} objects may use memory more efficiently than spectral 
#'   objects of the classes for collections of spectra defined in package 
#'   'photobiology' as wavelengths are assumed to be the same for all member 
#'   spectra, and stored only once while this assumption is not made for 
#'   collections of spectra, allowing different wavelengths and lengths for the 
#'   component spectra. Wavelengths are stored for each spectrum, but as 
#'   spectral classes are derived from 'tbl_df' in many cases no redundant 
#'   copies of wavelength data will be made in memory in spite of the more 
#'   flexible semantics of the objects.
#'   
#' @section Warning!: Always check the sanity of the imported or exported data 
#'   values, as guessing is needed when matching the different classes, and the 
#'   functions defined here are NOT guaranteed to return valid data wihtout help
#'   from the user through optional function arguments.
#'   
#' @param x rspec object
#' @param member.class character One of the spectrum classes defined in package 
#'   'photobiology'.
#' @param spct.data.var character The name to be used for the 'spc' data when 
#'   constructing the spectral objects.
#' @param multiplier numeric A multiplier to be applied to the 'rspc' data to do
#'   unit or scale conversion.
#' @param ... currently ignored.
#'   
#' @export
#' 
#' @examples 
#' # example run only if 'pavo' is available
#' if (requireNamespace("pavo", quietly = TRUE)) {
#'   library(pavo)
#'   data(sicalis, package = "pavo")
#'   sicalis.mspct <- rspec2mspct(sicalis)
#'   class(sicalis.mspct)
#' 
#'   data(teal, package = "pavo")
#'   teal.spct <- rspec2spct(teal)
#'   class(teal.spct)
#'   levels(teal.spct[["spct.idx"]])
#'   angles <- seq(from = 15, to = 75, by = 5) # from teal's documentation
#'   teal.spct[["angle"]] <- angles[as.numeric(teal.spct[["spct.idx"]])]
#'   teal.spct
#' }
#' 
rspec2mspct <- function(x, 
                        member.class = "reflector_spct", 
                        spct.data.var = "Rpc", 
                        multiplier = 1,
                        ...) {
  if (requireNamespace("pavo", quietly = TRUE)) {
  stopifnot(inherits(x, "rspec"))
  spct.names <- colnames(x)[-1]
  z <- photobiology::split2mspct(x = x, 
                                 member.class = member.class, 
                                 spct.data.var = spct.data.var,
                                 w.length.var = "wl")
  names(z) <- spct.names
  comment(z) <- paste('Converted from "pavo::rspec" object\n',
                      'dim: ', 
                      paste(names(dim(x)), dim(x), collapse = " "),
                      '\ncolnames: ', paste(colnames(x), collapse = ", "),
                      sep = "")
  } else {
    warning("Package 'pavo' needs to be installed.")
    z <- call(sub("spct", "mscpt", member.class))
  }
  z
}

#' @rdname rspec2mspct
#'
#' @export
#'
rspec2spct <- function(x, multiplier = 1, ...) {
  y <- rspec2mspct(x, multiplier = multiplier, ...)
  z <- photobiology::rbindspct(y)
  comment(z) <- comment(x[[1]])
  z
}


# colorSpec ---------------------------------------------------------------

#' Convert 'colorSpec::colorSpec' objects
#' 
#' Convert 'colorSpec::colorSpec' objects into spectral objects (xxxx_spct,
#' xxxx_mspct) as defined in package 'photobiology' and vice versa preserving
#' as much information as possible.
#' 
#' @details Objects of class \code{colorSpec::colorSpec} do not contain metadata
#' or class data from which the units of expression could be obtained. When
#' using function \code{colorSpec2mspct} the user needs to use parameter
#' \code{multiplier} to convert the data to what is expected by the object
#' constructors defined in package 'photobiology' but should only rarely need to
#' use parameter \code{spct.data.var} to select the quantity.
#'
#' \code{colorSpec::colorSpec} objects may use memory more efficiently than
#' spectral objects of the classes for collections of spectra defined in package
#' 'photobiology' as wavelengths are assumed to be the same for all member
#' spectra, and stored only once while this assumption is not made for
#' collections of spectra, allowing different wavelengths and lengths for the
#' component spectra. When using \code{as.colorSpec} methods to convert
#' collections of spectra into \code{colorSpec} objects, if the wavelengths of
#' the individual spectra differe, only the shared range of wavelengths is
#' retained and within the this range, wavelngth values are made consistent by
#' interpolation.
#'   
#' @note In \code{generic_mspct} objects, wavelengths are stored for each
#'   spectrum, individual \code{generic_spct} objects. However, as spectral
#'   classes are derived from 'tbl_df' in many cases no redundant copies of
#'   wavelength data will be made in memory in spite of the more flexible
#'   semantics of the objects.
#' 
#' @section Warning!: Always check the sanity of the imported or exported data
#'   values, as guessing is needed when matching the different classes, and the
#'   functions defined here are NOT guaranteed to return valid data without help
#'   from the user through optional function arguments.
#' 
#' @param x colorSpec object
#' @param multiplier numeric A multiplier to be applied to the 'spc' data to do
#'   unit or scale conversion.
#' @param ... currently ignored.
#' 
#' @export
#' 
#' @examples
#' # example run only if 'colorSpec' is available
#' if (requireNamespace("colorSpec", quietly = TRUE)) {
#'   library(colorSpec)
#'   colorSpec2mspct(Fs.5nm)
#'   colorSpec2spct(Fs.5nm)
#'   colorSpec2mspct(C.5nm)
#'   colorSpec2spct(C.5nm)
#' } 
#' 
colorSpec2mspct <- function(x, multiplier = 1, ...) {
  if (requireNamespace("colorSpec", quietly = TRUE)) {
    stopifnot(inherits(x, "colorSpec"))
    stopifnot(multiplier > 0)
    spct.type <- colorSpec::type(x)
    spct.quantity <- colorSpec::quantity(x)
    spct.metadata <- attr(x, "metadata", exact = TRUE)
    comment.spct <- paste('Converted from "colorSpec::colorSpec" object\n',
                          '"type": ', spct.type, "\n",
                          '"quantity": ', spct.quantity, "\n",
                          '"metadata path": ', spct.metadata[["path"]], "\n",
                          '"metadata header": \n', 
                          paste(spct.metadata[["header"]], collapse = "\n"), sep = "")
    y <- as.data.frame(as.matrix(x) * multiplier)
    y[["w.length"]] <- colorSpec::wavelength(x)
    if (spct.type == "light") {
      if (colorSpec::is.radiometric(x)) {
        z <- photobiology::split2source_mspct(y, 
                                              spct.data.var = "s.e.irrad")
      } else if (colorSpec::is.actinometric(x)) {
        z <- photobiology::split2source_mspct(y, 
                                              spct.data.var = "s.q.irrad")
      } else {
        stop("unkown 'quantity': ", spct.quantity)
      }
    } else if (spct.type == 'responsivity.light') {
      if ( colorSpec::is.radiometric(x)) {
        z <- photobiology::split2response_mspct(y, 
                                                spct.data.var = "s.e.response")
      } else if (colorSpec::is.actinometric(x)) {
        z <- photobiology::split2response_mspct(y, 
                                                spct.data.var = "s.q.response")
      } else {
        stop("unkown 'quantity': ", spct.quantity)
      }
    } else if (spct.type == 'material') {
      if (spct.quantity == 'reflectance') {
        z <- photobiology::split2reflector_mspct(y, 
                                                 spct.data.var = "Rfr")
      } else if (spct.quantity == 'transmittance') {
        z <- photobiology::split2filter_mspct(y, 
                                              spct.data.var = "Tfr")
      } else if (spct.quantity == 'absorbance') {
        z <- photobiology::split2filter_mspct(y, 
                                              spct.data.var = "A")
      } else {
        stop("unkown 'quantity': ", spct.quantity)
      }
    } else {
      return(list())
      # z <- photobiology::split2generic_mspct(y, 
      #                                        spct.data.var = spct.quantity)
    }
    comment(z) <- comment.spct
    scaled <- ifelse(multiplier == 1, FALSE, multiplier)
    photobiology::setScaled(z, scaled)
  } else {
    warning("Package 'colorSpec' needs to be installed.")
    NA
  }
  z
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.source_spct.colorSpec <- function(x, multiplier = 1, ...) {
  if (requireNamespace("colorSpec", quietly = TRUE)) {
    stopifnot(colorSpec::type(x) == "light")
    colorSpec2spct(x = x, 
                   multiplier = multiplier, 
                   ...)
  } else {
    warning("Package 'colorSpec' needs to be installed.")
    photobiology::source_spct()
  }
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.source_mspct.colorSpec <- function(x, multiplier = 1, ...) {
  if (requireNamespace("colorSpec", quietly = TRUE)) {
    stopifnot(colorSpec::type(x) == "light")
  colorSpec2spct(x = x, 
                 multiplier = multiplier, 
                 ...)
  } else {
    warning("Package 'colorSpec' needs to be installed.")
    photobiology::source_mspct()
  }
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.response_spct.colorSpec <- function(x, multiplier = 1, ...) {
  if (requireNamespace("colorSpec", quietly = TRUE)) {
    stopifnot(colorSpec::type(x) == "responsivity.light")
  colorSpec2spct(x = x, 
                 multiplier = multiplier, 
                 ...)
  } else {
    warning("Package 'colorSpec' needs to be installed.")
    photobiology::response_spct()
  }
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.response_mspct.colorSpec <- function(x, multiplier = 1, ...) {
  if (requireNamespace("colorSpec", quietly = TRUE)) {
    stopifnot(colorSpec::type(x) == "responsivity.light")
  colorSpec2mspct(x = x, 
                  multiplier = multiplier, 
                  ...)
  } else {
    warning("Package 'colorSpec' needs to be installed.")
    photobiology::response_mspct()
  }
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.filter_spct.colorSpec <- function(x, multiplier = 1, ...) {
  if (requireNamespace("colorSpec", quietly = TRUE)) {
    stopifnot(colorSpec::type(x) == "material" &&
              colorSpec::quantity(x) %in% c("absorbance", "transmittance"))
  colorSpec2spct(x = x, 
                 multiplier = multiplier, 
                 ...)
  } else {
    warning("Package 'colorSpec' needs to be installed.")
    photobiology::filter_spct()
  }
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.filter_mspct.colorSpec <- function(x, multiplier = 1, ...) {
  if (requireNamespace("colorSpec", quietly = TRUE)) {
    stopifnot(colorSpec::type(x) == "material" &&
                colorSpec::quantity(x) %in% c("absorbance", "transmittance"))
    colorSpec2mspct(x = x, 
                    multiplier = multiplier, 
                    ...)
  } else {
    warning("Package 'colorSpec' needs to be installed.")
    photobiology::filter_mspct()
  }
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.reflector_spct.colorSpec <- function(x, multiplier = 1, ...) {
  if (requireNamespace("colorSpec", quietly = TRUE)) {
    stopifnot(colorSpec::type(x) == "material" &&
                colorSpec::quantity(x) == "reflectance")
    colorSpec2spct(x = x, 
                   multiplier = multiplier, 
                   ...)
  } else {
    warning("Package 'colorSpec' needs to be installed.")
    photobiology::reflector_spct()
  }
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.reflector_mspct.colorSpec <- function(x, multiplier = 1, ...) {
  if (requireNamespace("colorSpec", quietly = TRUE)) {
    stopifnot(colorSpec::type(x) == "material" &&
                colorSpec::quantity(x) == "reflectance")
    colorSpec2mspct(x = x, 
                    multiplier = multiplier, 
                    ...)
  } else {
    warning("Package 'colorSpec' needs to be installed.")
    photobiology::reflector_mspct()
  }
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.chroma_mspct.colorSpec <- function(x, multiplier = 1, ...) {
  if (requireNamespace("colorSpec", quietly = TRUE)) {
    colorSpec2mspct(x = x, 
                    multiplier = multiplier, 
                    ...)
  } else {
    warning("Package 'colorSpec' needs to be installed.")
    photobiology::chroma_mspct()
  }
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
colorSpec2spct <- function(x, multiplier = 1, ...) {
  y <- colorSpec2mspct(x, multiplier = multiplier, ...)
  if (length(y) < 2) {
    z <- y[[1]]
  } else {
    z <- photobiology::rbindspct(y)
  }
  comment(z) <- comment(x[[1]])
  z
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
colorSpec2chroma_spct <- function(x, multiplier = 1, ...) {
  spct.type <- colorSpec::type(x)
  spct.quantity <- colorSpec::quantity(x)
  stopifnot(spct.quantity == 'energy->neural')
  stopifnot(colorSpec::numSpectra(x) == 3)
  stopifnot(sort(tolower(names(x))) == c("x", "y", "z"))
  stopifnot(multiplier > 0)
  spct.metadata <- attr(x, "metadata", exact = TRUE)
  comment.spct <- paste('Converted from "colorSpec::colorSpec" object\n',
                        '"type": ', spct.type, "\n",
                        '"quantity": ', spct.quantity, "\n",
                        '"metadata path": ', spct.metadata[["path"]], "\n",
                        '"metadata header": \n', 
                        paste(spct.metadata[["header"]], collapse = "\n"), sep = "")
  y <- as.data.frame(as.matrix(x) * multiplier)
  y[["w.length"]] <- colorSpec::wavelength(x)
  z <- photobiology::as.chroma_spct(y)
  comment(z) <- comment.spct
  z
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.chroma_spct.colorSpec <- function(x, 
                                     multiplier = 1, 
                                     ...) {
  colorSpec2spct(x = x, 
                 multiplier = multiplier, 
                 ...)
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.chroma_mspct.colorSpec <- function(x, 
                                      multiplier = 1, 
                                      ...) {
  colorSpec2mspct(x = x, 
                  multiplier = multiplier, 
                  ...)
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.generic_spct.colorSpec <- function(x, 
                                      multiplier = 1, 
                                      ...) {
  force(x)
  colorSpec2spct(x = x,
                 multiplier = multiplier,
                 ...)
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
as.generic_mspct.colorSpec <- function(x, 
                                       multiplier = 1, 
                                       ...) {
  force(x)
  colorSpec2mspct(x = x, multiplier = multiplier, ...)
}

#' @rdname as.colorSpec
#'   
#' @param spct.data.var character The name of the variable to read spectral data
#'   from.
#'   
#' @export
#' 
mspct2colorSpec <- function(x, 
                            spct.data.var = NULL,
                            multiplier = 1,
                            ...) {
  stopifnot(photobiology::is.any_mspct(x))
  if (length(x) > 1) {
    # colorSpec expects all spectra to share the same wavelengths
    wl.range <- c(photobiology::wl_range(x))
    wl.stepsize <- photobiology::wl_stepsize(x)
    wl.ranges.consistent <- 
      length(unique(wl.range[["min.wl"]])) == 1 &&
      length(unique(wl.range[["max.wl"]])) == 1
    wl.stepsizes.consistent <- length(unique(wl.stepsize[["min.step.wl"]])) == 1 &&
      length(unique(wl.stepsize[["max.step.wl"]])) == 1
    
    if (!wl.ranges.consistent || !wl.stepsizes.consistent) {
      # overlapping range
    wl.range.out <- c(max(wl.range[["min.wl"]]),
                      min(wl.range[["max.wl"]]))
    wl.stepsize.out <- stats::median(wl.stepsize[["min.step.wl"]]) / 2
    # we try to find a nearby "nice" stepsize
    wl.stepsize.out <- ifelse(wl.stepsize.out >= 1,
                              trunc(wl.stepsize.out),
                              ifelse(wl.stepsize.out >= 0.25,
                              trunc(wl.stepsize.out * 4) / 4,
                              round(wl.stepsize.out, digits = 2)))
    wl.out <- seq(from = wl.range.out[1],
                  to = wl.range.out[2],
                  by = wl.stepsize.out)
    x <- photobiology::interpolate_mspct(x,
                                         w.length.out = wl.out,
                                         fill = NULL)
    message("Spectra interpolated and trimmed as wavelengths differed.")
    }
  }
  #  warning("Deprecated: please use as.colorSpec() instead.")
  class.mspct <- class(x)[1]
  comment.mspct <- comment(x)
  comment.mspct <- 
    paste('Converted from "', class.mspct, '" object\n',
          comment.mspct, sep = "")
  if (class.mspct == "source_mspct") {
    if (is.null(spct.data.var)) {
      x <- photobiology::q2e(x, action = "replace")
      spct.data.var <- "s.e.irrad"
    }
    quantity <- 'energy'
  } else if (class.mspct == "response_mspct") {
    if (is.null(spct.data.var)) {
      x <- photobiology::q2e(x, action = "replace")
      spct.data.var <- "s.e.response"
    }
    quantity <- 'energy->neural'
  } else if (class.mspct == "filter_mspct") {
    if (is.null(spct.data.var)) {
      x <- photobiology::A2T(x, action = "replace")
      spct.data.var <- "Tfr"
    }
    quantity <- 'transmittance'
  } else if (class.mspct == "reflector_mspct") {
    if (is.null(spct.data.var)) {
      spct.data.var <- "Rfr"
    }
    quantity <- 'reflectance'
  } else if (class.mspct == "object_mspct") {
    warning("automatic conversion from 'object_spct' to 'colorSpec' not possible")
    return(colorSpec::colorSpec(data = numeric(), 
                                wavelength = numeric()))
  }
  spct.names <- names(x)
  spct.selector <- rep(TRUE, length.out = length(x))
  for (i in seq_along(x)) {
    temp <- x[[i]]
    s.column <- temp[[spct.data.var]]
    wl.current <- temp[["w.length"]]
    if (i == 1L) {
      mat <- s.column
      wl.prev <- wl.current
    } else {
      if (!all(wl.current == wl.prev)) {
        spct.selector[i] <- FALSE
        next()
      }
      mat <- cbind(mat, s.column) # as column!
    }
  }
  if (!is.null(dim(mat))) {
    colnames(mat) <- spct.names[spct.selector]
    add.name <- FALSE
  } else {
    add.name <- TRUE
  }
  z <- colorSpec::colorSpec(data = mat * multiplier, 
                            wavelength = wl.prev,
                            quantity = quantity
  )
  if (add.name) {
    colorSpec::specnames(z) <- spct.names[spct.selector]
  }
  attr(z, "metadata") <- list(path = "", header = comment.mspct)
  z
}

#' @rdname as.colorSpec
#' 
#' @export
#' 
spct2colorSpec <- function(x, 
                           spct.data.var = NULL,
                           multiplier = 1,
                           ...) {
  y <- photobiology::generic_mspct(list(x), class = class(x)[1])
  mspct2colorSpec(y,
                  spct.data.var,
                  multiplier)
}

#' @rdname as.colorSpec
#' 
#' @export
#' 
chroma_spct2colorSpec <- function(x, 
                                  spct.data.var = NULL,
                                  multiplier = 1,
                                  ...) {
  stopifnot(photobiology::is.chroma_spct(x))
  colorSpec::colorSpec(data = as.matrix(x[ , c("x", "y", "z")]) * multiplier,
                       wavelength = x[["w.length"]],
                       quantity = 'power->neural')
}

#' Convert into 'colorSpec::colorSpec' objects
#' 
#' Convert spectral objects (xxxx_spct, xxxx_mspct) as defined in package 
#' 'photobiology' into colorSpec objects preserving as much information as
#'  possible.
#' 
#' @note Objects of class \code{colorSpec::colorSpec} do not contain metadata or
#'   class data from which the units of expression could be obtained. When using
#'   this function the user needs to use parameter \code{multiplier} to convert 
#'   the data to what is expected by the object constructors defined in package 
#'   'photobiology' but should only rarely need to use parameter
#'   \code{spct.data.var} to select the quantity.
#'   
#'   \code{colorSpec::colorSpec} objects may use memory more efficiently than
#'   spectral objects of the classes for collections of spectra defined in
#'   package 'photobiology' as wavelengths are assumed to be the same for all
#'   member spectra, and stored only once while this assumption is not made for
#'   collections of spectra, allowing different wavelengths and lengths for the
#'   component spectra. Wavelengths are stored for each spectrum, but as
#'   spectral classes are derived from 'tbl_df' in many cases no redundant
#'   copies of wavelength data will be made in memory in spite of the more
#'   flexible semantics of the objects.
#' 
#' @section Warning!: Always check the sanity of the returned data values, as
#'   guessing is needed when matching the different classes, and the functions
#'   defined here are NOT guaranteed to return valid data without help from the
#'   user through optional function arguments.
#' 
#' @param x R object
#' @param spct.data.var character The name of the variable to read spectral data
#'   from.
#' @param multiplier numeric A multiplier to be applied to the 'spc' data to do
#'   unit or scale conversion.
#' @param ... currently ignored.
#' 
#' @name as.colorSpec
#' 
#' @importFrom colorSpec as.colorSpec
#' 
#' @export
#' 
#' @examples 
#' 
#' if (requireNamespace("colorSpec", quietly = TRUE)) {
#'   library(colorSpec)
#'   as.colorSpec(polyester.spct)
#'   as.colorSpec(sun.spct)
#'   as.colorSpec(filter_mspct(list(polyester.spct, yellow_gel.spct)))
#' }
#' 
as.colorSpec.generic_mspct <- function(x, 
                                       spct.data.var = NULL,
                                       multiplier = 1,
                                       ...) {
  mspct2colorSpec(x = x,
                  spct.data.var = spct.data.var,
                  multiplier = multiplier,
                  ...)
}

#' @describeIn as.colorSpec
#' 
#' @export
#' 
as.colorSpec.generic_spct <- function(x, 
                                      spct.data.var = NULL,
                                      multiplier = 1,
                                      ...) {
  spct2colorSpec(x = x,
                 spct.data.var = spct.data.var,
                 multiplier = multiplier,
                 ...)
}


#' @describeIn as.colorSpec
#' 
#' @export
#' 
as.colorSpec.chroma_spct <- function(x, 
                                     spct.data.var = NULL,
                                     multiplier = 1,
                                     ...) {
  colorSpec::colorSpec(data = as.matrix(x[ , c("x", "y", "z")]) * multiplier,
                       wavelength = x[["w.length"]],
                       quantity = 'energy->neural')
}

