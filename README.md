
# photobiologyInOut <img src="man/figures/logo.png" align="right" width="120" />

<!-- badges: start -->

[![CRAN
version](https://www.r-pkg.org/badges/version-last-release/photobiologyInOut)](https://cran.r-project.org/package=photobiologyInOut)
[![cran
checks](https://badges.cranchecks.info/worst/photobiologyInOut.svg)](https://cran.r-project.org/web/checks/check_results_photobiologyInOut.html)
[![photobiologyInOut status
badge](https://aphalo.r-universe.dev/badges/photobiologyInOut)](https://aphalo.r-universe.dev/photobiologyInOut)
[![R-CMD-check](https://github.com/aphalo/photobiologyInOut/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aphalo/photobiologyInOut/actions/workflows/R-CMD-check.yaml)
[![Documentation](https://img.shields.io/badge/documentation-photobiologyInOut-informational.svg)](https://docs.r4photobiology.info/photobiologyInOut/)
[![doi](https://img.shields.io/badge/doi-10.32614/CRAN.package.photobiologyInOut-blue.svg)](https://doi.org/10.32614/CRAN.package.photobiologyInOut)
<!-- badges: end -->

Package ‘**photobiologyInOut**’ provides functions for importing
spectral data from diverse sources including instrument-specific files,
spectral data output by solar-radiation simulation models, and files
from on-line data repositories. It also includes functions for
conversion of spectral data from/to classes defined in other R packages.
Package ‘**photobiologyInOut**’ complements other packages in the
‘**r4photobiology suite**’ by allowing reading and writing “foreign”
spectral data as well as reading non-spectral data acquired with data
loggers. The functions transfer metadata from file headers or companion
metadata files to attributes of the returned objects, in most cases both
the headers as text and parsed metadata extracted from them.

Data files from **spectrometers** from the following suppliers are
currently supported: CID Bio-Science, Wasatch Photonics, Avantes,
LI-COR, Macam Photonics, and Ocean Optics/Ocean Insight.

Data files from **data loggers** from the following suppliers are
currently supported: Campbell Scientific and YoctoPuce.

Data files from output by **radiation transfer models**: libRadtran and
TUV, including running simulations with the Quick TUV calculator.

Data objects of classes defined in **R packages**: ‘hyperSpec’ (2-way),
‘colorSpec’, ‘pavo’ and ‘fda.usc’ (2-way).

Data files downloaded from **repositories of spectral data**: ASTER
(NASA’s ECOSTRESS Spectral Library), FReD (Floral Reflectance Database)
and CIE’s (International Commision of Illumination) data tables,
corresponding to CIE standards.

Developing a package like this is a never-ending task as I have only a
limited sample of output files for testing and formats are quite
variable. The functions may not work with different software or firmware
versions used for acquiring spectral data from instruments. Even the
format of files can depend on the current locale and operating system.

This package is part of a suite of R packages for photobiological
calculations described at the
[r4photobiology](https://www.r4photobiology.info) web site.

## Warning

**The functions in this package work with the example files I have had
access to for testing, but they may not work with your own files as file
formats vary and are subject to revision.**

**PLEASE, BE VERY CAREFUL WHEN USING THIS PACKAGE. DO CHECK THAT UNITS
USED IN THE IMPORTED FILE ARE THOSE EXPECTED BY THESE FUNCTIONS AND THAT
THE VALUES IN THE RETRIEVED DATA ARE THOSE EXPECTED!**

*If the functions in this package do not work with your files, these
functions hopefully will be useful as examples for developing your own
functions. If you develop new functions or improve the existing ones,
please, do contribute them back to this project as a pull request at
GitHub or contact the maintainer. Sharing with the maintainer files that
fail to be read correctly is also of great help.*

## Installation

Installation of the most recent stable version from CRAN:

``` r
install.packages("photobiologyInOut")
```

Installation of the current unstable version from R-Universe CRAN-like
repository:

``` r
install.packages('photobiologyInOut', 
                 repos = c('https://aphalo.r-universe.dev', 
                           'https://cloud.r-project.org'))
```

Installation of the current unstable version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("aphalo/photobiologyInOut")
```

## Documentation

HTML documentation is available at
(<https://docs.r4photobiology.info/photobiologyInOut/>), including a
*User Guide*.

News on updates to the different packages of the ‘r4photobiology’ suite
are regularly posted at (<https://www.r4photobiology.info/>).

Two articles introduce the basic ideas behind the design of the suite
and its use: Aphalo P. J. (2015)
(<https://doi.org/10.19232/uv4pb.2015.1.14>) and Aphalo P. J. (2016)
(<https://doi.org/10.19232/uv4pb.2016.1.15>).

A book is under preparation, and the draft is currently available at
(<https://leanpub.com/r4photobiology/>).

A handbook written before the suite was developed contains useful
information on the quantification and manipulation of ultraviolet and
visible radiation: Aphalo, P. J., Albert, A., Björn, L. O., McLeod, A.
R., Robson, T. M., & Rosenqvist, E. (Eds.) (2012) Beyond the Visible: A
handbook of best practice in plant UV photobiology (1st ed., p. xxx +
174). Helsinki: University of Helsinki, Department of Biosciences,
Division of Plant Biology. ISBN 978-952-10-8363-1 (PDF),
978-952-10-8362-4 (paperback). PDF file available from
(<https://doi.org/10.31885/9789521083631>).

## Contributing

Pull requests, bug reports, and feature requests are welcome at
(<https://github.com/aphalo/photobiologyInOut>). Contribution of example
data files that could be supported in future versions will be very much
appreciated.

## Citation

If you use this package to produce scientific or commercial
publications, please cite according to:

``` r
citation("photobiologyInOut")
#> To cite package ‘photobiologyInOut’ in publications use:
#> 
#>   Aphalo, Pedro J. (2015) The r4photobiology suite. UV4Plants Bulletin,
#>   2015:1, 21-29. DOI:10.19232/uv4pb.2015.1.14
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     author = {Pedro J. Aphalo},
#>     title = {The r4photobiology suite},
#>     journal = {UV4Plants Bulletin},
#>     volume = {2015},
#>     number = {1},
#>     pages = {21-29},
#>     year = {2015},
#>     doi = {10.19232/uv4pb.2015.1.14},
#>   }
```

## License

© 2015-2025 Pedro J. Aphalo (<pedro.aphalo@helsinki.fi>). Released under
the GPL, version 2 or greater. This software carries no warranty of any
kind.
