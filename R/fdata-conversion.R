# FDA ---------------------------------------------------------------------

#' Convert spectra into 'fda.usc::fdata' objects
#' 
#' Convert spectral objects (xxxx_spct, xxxx_mspct) as defined in package
#' 'photobiology' into \code{fda.usc::fdata} objects, preserving as much
#' information as possible. As \code{fdata} objects can contain other kinds
#' of data, the reverse conversion is supported (experimentally) and mainly
#' for 'fdata' objects returned by the functional data analysis methods from
#' package \code{fda.usc} to spectral data previosuly exported in the
#' opposite direction.
#' 
#' @section Warning!: When converting multiple spectra, all the spectra to be
#' included in the \code{fdata} object must share the same wavelength values.
#' Spectra that do not fulfil this condition will be skipped. The data variable
#' needs also to be present in all individual spectra as no conversions are 
#' applied automatically by this function. If a different name, indicating a
#' different quantity or a different base of expression is encountered, the
#' affected spectrum is skipped with a warning.
#' 
#' @param x generic_mspct or generic_spct object or an object belonging to a
#'   derived class, or an object of class 'fdata' depending on the function.
#' @param spct.data.var character The name of the column containing data to
#'   export. If \code{NULL} the first spectral data column found is used.
#' @param multiplier numeric A multiplier to be applied to the 'spc' data to do
#'   unit or scale conversion.
#' @param ... possibly additional named arguments passed to object constructors.
#' 
#' @export
#' 
#' @examples 
#' if (requireNamespace("fda.usc", quietly = TRUE)) {
#' # from spectra to fdata
#'   sun.fdata <- spct2fdata(sun.spct)
#'   str(sun.fdata)
#'   polyester.fdata <- spct2fdata(polyester.spct)
#'   str(polyester.fdata)
#' # from fdata to spectra
#'   fdata2spct(sun.fdata)
#'   fdata2spct(sun.fdata, drop.idx = TRUE)
#'   fdata2spct(polyester.fdata, drop.idx = TRUE)
#' }
#' 
#' @export
#' 
mspct2fdata <- function(x, 
                        spct.data.var = NULL,
                        multiplier = 1,
                        ...) {
  if (requireNamespace("fda.usc", quietly = TRUE)) {
    stopifnot(photobiology::is.any_mspct(x))
    if (is.null(spct.data.var )) {
      spct.data.var <- setdiff(colnames(x[[1]]), "w.length")[1]
    }
    w.length <- x[[1]][["w.length"]]
    spct.names <- names(x)
    spct.selector <- rep(TRUE, length(x))
    vec <- numeric()
    for (i in seq_along(x)) {
      temp <- x[[i]]
      if (!spct.data.var %in% colnames(temp)) {
        message("Missing column, skipping: ", spct.names[i])
        spct.selector[i] <- FALSE
        next()
      }
      if (!isTRUE(all.equal(w.length, x[[i]][["w.length"]]))) {
        message("Inconsistent 'w.length', skipping: ", spct.names[i])
        spct.selector[i] <- FALSE
        next()
      }
      s.column <- temp[[spct.data.var]] * multiplier
      vec <- c(vec, s.column)
    }
    
    mat <- matrix(vec, nrow = sum(spct.selector), dimnames = list(spct.names[spct.selector]), byrow = TRUE)
    z <- fda.usc::fdata(mdata = mat,
                        argvals = w.length,
                        names = list(main = "", 
                                     xlab = "w.length", 
                                     ylab = spct.data.var),
                        fdata2d = FALSE,
                        ...)
    comment(z) <- paste('Converted from object of class "', class(x)[1], '\"', sep = "")
    #    message("Created 'fdata' object containing ", sum(spct.selector), " spectra.")
    z
  } else {
    warning("Package 'fda.usc' needs to be installed.")
    NA
  }
}

#' @rdname mspct2fdata
#' 
#' @export
#' 
spct2fdata <- function(x, 
                       spct.data.var = NULL,
                       multiplier = 1,
                       ...) {
  stopifnot(is.any_spct(x))
  if (is.null(spct.data.var )) {
    spct.data.var <- setdiff(colnames(x), "w.length")[1]
  }
  stopifnot(spct.data.var %in% colnames(x))
  if (getMultipleWl(x) == 1L) {
    y <- photobiology::generic_mspct(list(x), class = class(x)[1])
  } else {
    y <- subset2mspct(x)
  }
  mspct2fdata(y, 
              spct.data.var = spct.data.var,
              multiplier = multiplier,
              ...)
}

#' @rdname mspct2fdata
#' 
#' @param member.class character Name of the class of the spectrum or of the
#'    members of the collection of spectra.
#' 
#' @export
#' 
fdata2spct <- function(x, 
                       multiplier = 1,
                       member.class = NULL,
                       drop.idx = FALSE,
                       ...) {
  if (requireNamespace("fda.usc", quietly = TRUE)) {
    stopifnot(inherits(x, "fdata"))
    if (!is.null(member.class) &&
        !member.class %in% photobiology::spct_classes()) {
      stop("Bad argument to 'member.class': ", member.class)
    }
    
    spct.data.var <- x[["names"]][["ylab"]]
    if (is.null(member.class)) {
      if (grepl("irrad", spct.data.var)) {
        member.class <- "source_spct"
      } else if (grepl("response", spct.data.var)) {
        member.class <- "response_spct"
      } else if (grepl("Tfr|Tpc|Afr|Apc", spct.data.var)) {
        member.class <- "filter_spct"
      } else if (grepl("^A$", spct.data.var)) {
        member.class <- "filter_spct"
      } else if (grepl("Rfr|Rpc", spct.data.var)) {
        member.class <- "reflector_spct"
      }
    } else if (!is.null(member.class) && 
               !member.class %in% photobiology::spct_classes()) {
      stop("Bad argument to 'member.class': ", member.class)
    }
    
    wavelength.var <- x[["names"]][["xlab"]]
    if (!x[["names"]][["xlab"]] == "w.length") {
      warning("Found ", wavelength.var, " instead of 'w.length' in 'fdata' object!")
    }
    
    spct.names <- rownames(x[["data"]])
    if (is.null(spct.names)) {
      if (!is.null(x[["names"]][["main"]])) {
        name.root <- x[["names"]][["main"]]
      } else {
        name.root <- "row"
      }
      spct.names <- paste(name.root, seq(1, nrow(x[["data"]])), sep = "_")
    }
    if (length(spct.names) == 1L && drop.idx) {
      args <- list(w.length = rep(x[["argvals"]], times = length(spct.names)),
                   as.vector(x[["data"]]),
                   multiple.wl = 1L,
                   ...)
    } else {
      args <- list(w.length = rep(x[["argvals"]], times = length(spct.names)),
                   as.vector(x[["data"]]),
                   spct.idx = rep(spct.names, each = length(x[["argvals"]])),
                   multiple.wl = length(spct.names),
                   idfactor = "spct.idx",
                   ...)
    }
    names(args)[2] <- spct.data.var
    z <- do.call(member.class, args = args)
    photobiology::what_measured(z) <- x[["names"]][["main"]]
    comment(z) <- paste("Converted from class 'fdata' into '", class(z), "'", sep = "")
    z
  } else {
    warning("Package 'fda.usc' needs to be installed.")
    NA
  }
}

#' @rdname mspct2fdata
#' 
#' @param drop.idx logical Flag indicating whether to drop or keep
#'   \code{idx.var} in the collection members.
#' 
#' @export
#' 
fdata2mspct <- function(x, 
                       multiplier = 1,
                       member.class = NULL,
                       drop.idx = FALSE,
                       ...) {
  y <- fdata2spct(x = x,
                  multiplier = multiplier,
                  member.class = member.class,
                  drop.idx = FALSE,
                  ...)
  photobiology::subset2mspct(y, drop.idx = drop.idx)
}
  
