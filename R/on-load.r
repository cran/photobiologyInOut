.onLoad <- function(libname, pkgname) {
  # depends on 'colorspec' version and on whether it has been attached or not
  if (!any(grepl("colorspec", names(.Options))) && 
      utils::packageVersion("colorSpec") < "1.5-0") {
    options(colorSpec.logformat = "%t{%H:%M:%OS3} %l %n::%f(). %m",
            colorSpec.loglevel = "WARN",
            colorSpec.stoponerror = TRUE)
  }
}