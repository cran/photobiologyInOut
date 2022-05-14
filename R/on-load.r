.onLoad <- function(libname, pkgname) {
  options(colorSpec.logformat = "%t{%H:%M:%OS3} %l %n::%f(). %m",
          colorSpec.loglevel = "WARN",
          colorSpec.stoponerror = FALSE)
}