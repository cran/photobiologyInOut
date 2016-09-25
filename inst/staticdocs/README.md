# photobiologyInOut #

Package '**photobiologyInOut**' provides functions for importing spectral data from diverse sources including instrument-specific files, as well as spectral data output by solar-radiation simulation models. It also includes functions for exchanging spectral data with other R packages. Package '**photobiologyInOut**' complements other packages in the '**r4photobiology suite**' by allowing reading and writing _foreign_ spectral data as well as reading data saved from data loggers.

Developing a package like this is a never-ending task as I have only a limited sample of output files for testing and formats are quite variable. The functions may not work with different software or firmware versions used for acquiring spectral data from instruments. Even the format of files can depend on the current locale and operating system.

Please, see the [r4photobiology](http://www.r4photobiology.info) web site for details on the suite.

## Warning ##

**The functions in this package work with the example files I have access to for testing, but they may not work with your own files as file formats vary.**

**PLEASE, BE VERY CAREFUL WHEN USING THIS PACKAGE. DO CHECK THAT UNITS USED IN THE IMPORTED FILE ARE THOSE EXPECTED BY THESE FUNCTIONS AND THAT THE VALUES ARE AS EXPECTED!**

_If they do not work with your files, they hopefully will be useful as examples for developing your own functions. If you develop new functions or improve the existing ones, please, do contribute them back to this project._
