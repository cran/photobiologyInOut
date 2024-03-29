% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read-tuv-file.r
\name{read_qtuv_txt}
\alias{read_qtuv_txt}
\title{Read Quick TUV output file.}
\usage{
read_qtuv_txt(
  file,
  ozone.du = NULL,
  label = NULL,
  tz = NULL,
  locale = readr::default_locale()
)
}
\arguments{
\item{file}{character string with the name of a text file.}

\item{ozone.du}{numeric Ozone column in Dobson units.}

\item{label}{character string, but if \code{NULL} the value of \code{file} is
used, and if \code{NA} the "what.measured" attribute is not set.}

\item{tz}{character Time zone is by default read from the file.}

\item{locale}{The locale controls defaults that vary from place to place. The
default locale is US-centric (like R), but you can use
\code{\link[readr]{locale}} to create your own locale that controls things
like the default time zone, encoding, decimal mark, big mark, and day/month
names.}
}
\value{
a source_spct object obtained by finding the center of wavelength
  intervals in the Quick TUV output file, and adding variables \code{zenith.angle} and
  \code{date}.
}
\description{
Reads and parses the header of a text file output by the Quick TUV on-line
web front-end at UCAR to extract the header and spectral data. The time field
is converted to a date.
}
\note{
The ozone column value used in the simulation cannot be retrieved from
the file. Tested only with Quick TUV versison 5.2 on 2018-07-30. This 
function can be expected to be robust to variations in the position of lines
in the imported file and resistent to the presence of extraneous text or
even summaries.
}
\references{
\url{https://www.acom.ucar.edu/Models/TUV/Interactive_TUV/}
}
