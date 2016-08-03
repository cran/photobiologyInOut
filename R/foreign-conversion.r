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

# R matrix ----------------------------------------------------------------

#' Convert a collection of spectra into a matrix
#' 
#' Convert an object of class \code{generic_mspct} or a derived class into an R
#' matrix with wavelengths saved as an attribute and spectral data in rows
#' or columns.
#' 
#' @note Only collections of spectra containing spectra with exactly the same
#' \code{w.length} values can by converted. If needed, the spectra can be
#' re-expressed before attempting the conversion to a matrix.
#' 
#' @param x generic_mspct object.
#' @param spct.data.var character The name of the variable containing the spectral data.
#' @param byrow logical. If FALSE (the default) the matrix is filled with the
#'   spectra stored by columns, otherwise the matrix is filled by rows.
#' @param ... currently ignored.
#' 
#' @section Warning!: This conversion preserves the spectral data but discards
#'   almost all the metadata contained in the spectral objects. In other words a
#'   matrix created with this function cannot be used to recreate the original
#'   object unless the same metadata is explicitly supplied when converting the
#'   matrix into new collection of spectra.
#' 
#' @export
#' 
mspct2mat <- function(x, 
                      spct.data.var,
                      byrow = attr(x, "mspct.byrow"),
                      ...) {
  stopifnot(is.any_mspct(x))
  if (length(x) == 0L) {
    return(matrix(numeric()))
  }
  spct.names <- names(x)
  spct.selector <- rep(TRUE, length(x))
  mat <- numeric()
  for (i in seq_along(x)) {
    temp <- x[[i]]
    s.column <- temp[[spct.data.var]]
    wl.current <- temp[["w.length"]]
    if (i == 1L) {
      wl.prev <- wl.current
    } 
    if (!all(wl.current == wl.prev) || length(s.column) == 0L) {
      spct.selector[i] <- FALSE
      next()
    }
    mat <- c(mat, s.column) # one long numeric vector
  }
  if (any(!spct.selector)) {
    warning("Spectra dropped: ", sum(!spct.selector), " out of ", length(spct.selector), ".")
  }
  if (byrow) {
    z <- matrix(mat, nrow = sum(spct.selector), byrow = byrow,
                dimnames = list(spct = c(spct.names[spct.selector]), 
                                w.length = wl.prev))
  } else {
    z <- matrix(mat, ncol = sum(spct.selector), byrow = byrow,
                dimnames = list(w.length = wl.prev,
                                spct = c(spct.names[spct.selector])))
  }
  attr(z, "w.length") <- wl.prev
  comment(z) <- comment(x)
  z
}

#' Convert a matrix into a collection of spectra
#' 
#' Convert an R object of class matrix into a \code{generic_mspct} or a derived 
#' class.
#' 
#' @note Only \code{matrix} objects that have rows or columns of the same length
#'   as the numeric vector of walengths supplied can be converted. The resulting
#'   spectra will be built using the constructors and subjected to the same
#'   checks as if built individually. Only collections with all members of the
#'   same class can be built with this function. Additional named arguments can
#'   be supplied to set the same metadata attributes to all the member spectra.
#'   In the case of square matrices, an explicit argument is needed for
#'   \code{byrow} making it good practice for scripts and package code to not
#'   rely on the automatic default.
#'   
#' @param x matrix object.
#' @param w.length numeric A vector of walength values sorted in strictly ascending
#'   order (nm).
#' @param spct.data.var character The name of the variable that will contain the 
#'   spectral data. This indicates what physical quantity is stored in the matrix 
#'   and the units of expression used.
#' @param member.class character The name of the class of the individual spectra
#'   to be constructed.
#' @param multiplier numeric A multiplier to be applied to the values in \code{x} to do
#'   unit or scale conversion.
#' @param byrow logical Flag indicating whether each spectrum is stored in a row or
#'   a column of the matrix. By default this value is set based on the length of
#'   \code{w.length} and the dimensions of \code{x}.
#' @param spct.names character Vector of names to be assigned to collection members,
#'   either of length 1, or with length equal to the number of spectra.
#' @param ... other arguments passed to the constructor of collection members.
#' 
#' @export
#' 
#' @examples 
#' 
#' x <- matrix(1:100, ncol = 2)
#' wl <- (301:350)
#' z <- mat2mspct(x, wl, "filter_spct", "Tpc")
#' 
#' x <- matrix(1:100, nrow = 2, byrow = TRUE)
#' wl <- (301:350)
#' z <- mat2mspct(x, wl, "filter_spct", "Tpc", byrow = TRUE, spct.name = c("A", "B"))
#' 
mat2mspct <- function(x,
                      w.length,
                      member.class,
                      spct.data.var,
                      multiplier = 1,
                      byrow = NULL,
                      spct.names = "spct_",
                      ...) {
  stopifnot(is.matrix(x))
  if (length(spct.names) == 0) {
    spct.names = "spct"
  }
  if (is.null(byrow)) {
    if (nrow(x) == ncol(x)) {
      stop("For square matrices an argument for 'byrow' is mandatory")
    } else if (nrow(x) == length(w.length)) {
      byrow <- FALSE
    } else if (ncol(x) == length(w.length)) {
      byrow <- TRUE
    } else {
      stop("Length of 'w.length' vector is different to that of spectral data.")
    }
  }
  # spc data (spectra) can be stored as rows or as colums in a matrix, 
  # consequently if stored by rows we transpose the matrix.
  if (byrow) {
    x <- t(x)
  }
  ncol <- ncol(x)
  stopifnot(nrow(x) == length(w.length))
  y <- cbind(w.length, x * multiplier)
  y <- tibble::as_data_frame(y)
  if (length(spct.names) == ncol) {
    colnames(y) <- c("w.length", spct.names)
  } else {
    colnames(y) <- c("w.length", paste(spct.names[1], 1:ncol, sep = ""))
  }
  z <- split2mspct(x = y, 
                   member.class = member.class, 
                   spct.data.var = spct.data.var,
                   ncol = ncol,
                   ...)
  comment(z) <- paste('Converted from an R "matrix" object\n',
                      'with ', length(z), ' spectra stored ',
                      ifelse(byrow, "in rows.", "in columns."),
                      sep = "")
  z
}

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
