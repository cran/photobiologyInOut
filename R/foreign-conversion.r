# The "export" functions are not defined as S3 methods named as.hyperSpec and
# as.colorSpec as would have been natural, because these generics are not
# defined in the respective packages. Defining them as S3 generics in this
# independently developed package would be looking for future trouble. The
# "import" functions could have been defined as "as...." S3 methods, but
# naming symmetry between import and export functions seemed more natural.
# Another reason is that a single class "colorSpec" corresponds to several
# different classes from package "photobiology" and consequently the class
# the object is to be converted to depends on the information stored as
# attributes in a "colorSpec" object rather than on the class it belongs to. 

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
#' @param multiplier numeric A multiplier to be applied to the 'spc' data to
#'   do unit or
#'   scale conversion. For example "a.u." units in some examples in package
#'   'hyperSpec' seem to have scale factors applied.
#' @param ... currently ignored.
#' 
#' @export
#' 
#' @examples 
#' 
#' library(hyperSpec)
#' data(laser)
#' wl(laser) <- 
#' list (wl = 1e7 / (1/405e-7 - wl (laser)),
#'       label = expression (lambda / nm))
#' laser.mspct <- hyperSpec2mspct(laser, "source_spct", "s.e.irrad")
#' class(laser.mspct)
#' 
hyperSpec2mspct <- function(x, 
                            member.class, 
                            spct.data.var,
                            multiplier = 1,
                            ...) {
  stopifnot(inherits(x, "hyperSpec"))
  # spc data (spectra) are stored as rows in a matrix, consequently
  # we transpose the matrix so that each spectrum is in a column
  y <- cbind(hyperSpec::wl(x), t(x$spc) * multiplier) 
  colnames(y) <- c("w.length", paste("spc", 1:nrow(x), sep = ""))
  y <- tibble::as_data_frame(y)
  z <- split2mspct(x = y, 
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
  z
}

#' @rdname hyperSpec2mspct
#'
#' @export
#'
hyperSpec2spct <- function(x, multiplier = 1, ...) {
  y <- hyperSpec2mspct(x, multiplier = multiplier, ...)
  z <- rbindspct(y)
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
}

#' @rdname hyperSpec2mspct
#' 
#' @export
#' 
spct2hyperSpec <- function(x, 
                           spct.data.var = NULL,
                           multiplier = 1,
                           ...) {
  y <- generic_mspct(list(x), class = class(x)[1])
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
#' 
#' library(pavo)
#' data(sicalis)
#' sicalis.mspct <- rspec2mspct(sicalis)
#' class(sicalis.mspct)
#' 
#' data(teal)
#' teal.spct <- rspec2spct(teal)
#' class(teal.spct)
#' levels(teal.spct[["spct.idx"]])
#' angles <- seq(from = 15, to = 75, by = 5) # from teal's documentation
#' teal.spct[["angle"]] <- angles[as.numeric(teal.spct[["spct.idx"]])]
#' teal.spct
#' 
rspec2mspct <- function(x, 
                        member.class = "reflector_spct", 
                        spct.data.var = "Rpc", 
                        multiplier = 1,
                        ...) {
  stopifnot(pavo::is.rspec(x))
  spct.names <- colnames(x)[-1]
  z <- split2mspct(x = x, 
                   member.class = member.class, 
                   spct.data.var = spct.data.var,
                   w.length.var = "wl")
  names(z) <- spct.names
  comment(z) <- paste('Converted from "pavo::rspec" object\n',
                      'dim: ', 
                      paste(names(dim(x)), dim(x), collapse = " "),
                      '\ncolnames: ', paste(colnames(x), collapse = ", "),
                      sep = "")
  z
}

#' @rdname rspec2mspct
#'
#' @export
#'
rspec2spct <- function(x, multiplier = 1, ...) {
  y <- rspec2mspct(x, multiplier = multiplier, ...)
  z <- rbindspct(y)
  comment(z) <- comment(x[[1]])
  z
}



# colorSpec ---------------------------------------------------------------

#' Convert 'colorSpec::colorSpec' objects
#' 
#' Convert 'colorSpec::colorSpec' objects into spectral objects (xxxx_spct,
#' xxxx_mspct) as defined in package 'photobiology' and vice veersa preserving
#' as much information as possible.
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
#' @section Warning!: Always check the sanity of the imported or exported data
#'   values, as guessing is needed when matching the different classes, and the
#'   functions defined here are NOT guaranteed to return valid data wihtout help
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
#' 
#' library(colorSpec)
#' colorSpec2mspct(Fs.5nm)
#' colorSpec2spct(Fs.5nm)
#' colorSpec2mspct(C.5nm)
#' colorSpec2spct(C.5nm)
#' 
colorSpec2mspct <- function(x, multiplier = 1, ...) {
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
    if (spct.quantity == 'power') {
      z <- photobiology::split2source_mspct(y, 
                                            spct.data.var = "s.e.irrad")
    } else if (spct.quantity == 'photons/sec') {
      z <- photobiology::split2source_mspct(y, 
                                            spct.data.var = "s.q.irrad")
    } else {
      stop("unkown 'quantity': ", spct.quantity)
    }
  } else if (spct.type == 'responsivity.light') {
    if (spct.quantity %in% c('power->electrical', 'power->neural', 'power->action')) {
      z <- photobiology::split2response_mspct(y, 
                                              spct.data.var = "s.e.response")
    } else if (spct.quantity %in% c('photons->electrical', 'photons->neural', 'photons->action')) {
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
  setScaled(z, scaled)
  z
}

#' @rdname colorSpec2mspct
#'
#' @export
#'
colorSpec2spct <- function(x, multiplier = 1, ...) {
  y <- colorSpec2mspct(x, multiplier = multiplier, ...)
  z <- rbindspct(y)
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
  stopifnot(spct.quantity == 'power->neural')
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
#' @param spct.data.var character The name of the variable to read spectral data
#'   from.
#'   
#' @export
#' 
mspct2colorSpec <- function(x, 
                            spct.data.var = NULL,
                            multiplier = 1,
                            ...) {
  stopifnot(is.any_mspct(x))
  class.mspct <- class(x)[1]
  comment.mspct <- comment(x)
  comment.mspct <- 
    paste('Converted from "', class.mspct, '" object\n',
          comment.mspct, sep = "")
  if (class.mspct == "source_mspct") {
    if (is.null(spct.data.var)) {
      x <- q2e(x, action = "replace")
      spct.data.var <- "s.e.irrad"
    }
    quantity <- 'power'
  } else if (class.mspct == "response_mspct") {
    if (is.null(spct.data.var)) {
      x <- q2e(x, action = "replace")
      spct.data.var <- "s.e.response"
    }
    quantity <- 'power->action'
  } else if (class.mspct == "filter_mspct") {
    if (is.null(spct.data.var)) {
      x <- A2T(x, action = "replace")
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

#' @rdname colorSpec2mspct
#' 
#' @export
#' 
spct2colorSpec <- function(x, 
                           spct.data.var = NULL,
                           multiplier = 1,
                           ...) {
  y <- generic_mspct(list(x), class = class(x)[1])
  mspct2colorSpec(y,
                  spct.data.var,
                  multiplier)
}

#' @rdname colorSpec2mspct
#' 
#' @export
#' 
chroma_spct2colorSpec <- function(x, 
                                  spct.data.var = NULL,
                                  multiplier = 1,
                                  ...) {
  stopifnot(is.chroma_spct(x))
  colorSpec::colorSpec(data = as.matrix(x[ , c("x", "y", "z")]) * multiplier,
                       wavelength = x[["w.length"]],
                       quantity = 'power->neural')
}
