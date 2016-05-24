# README for phootbiologyInOut #

`photobiologyInOut` is an R package providing functions for importing spectral data from diverse sources including instrument-specific files, spectral data output by solar-radiation simulation models, and from spectral objects of classes defined by other R packages. It complements other packages in the  R suite for _photobiology_ by allowing reading and writing "foreign" spectral data.

This package is at a rather early stage of development, and I have only a limited sample of output files for testing. The functions may not work with different software or firmware versions used for generating spectral data. Even the format of files can depend on the current locale and operating system. The locale used for reading in data can be set trough a function argument as the files being read-in may have been created under a different locale setting.

Please, see the web site [R4Photobiology](http://www.r4photobiology.info) for details on other packages available as part of the suite, and on how to install them.

_The functions in this package work with the example files I have used for testing, but they may not work with your own files as file formats may vary._

__PLEASE, BE VERY CAREFUL WHEN USING THIS PACKAGE. DO CHECK THAT UNITS USED IN THE IMPORTED FILE ARE THOSE EXPECTED BY THESE FUNCTIONS AND THAT THE VALUES ARE AS EXPECTED!__

_If they do not work with your files, they hopefully are useful as examples for developing your own functions. If you develop new functions or improve the existing ones, please, do contribute them back to this project._
