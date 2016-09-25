## ----setup, include=FALSE, cache=FALSE---------------
library(knitr)
opts_chunk$set(fig.path='figure/pos-', fig.align='center', fig.show='hold',
               fig.width=7, fig.height=6, size="footnotesize")
options(replace.assign = TRUE, width = 55,
        warnPartialMatchAttr = FALSE,
        warnPartialMatchDollar = FALSE,
        warnPartialMatchArgs = FALSE)

## ----example-0-hiden, eval=TRUE, include=FALSE-------
# this may be needed in some geographic locations as some Windows TZ strings are
# not recognized by all versions of R
Sys.setenv(TZ = 'UTC')
library(photobiology)
library(photobiologyWavebands)
library(photobiologyInOut)
library(lubridate)
library(ggplot2)
library(ggmap)
library(ggspectra)
library(hyperSpec)
library(colorSpec)
library(pavo)
library(readr)

## ----own-set-up, echo=FALSE, include=FALSE-----------
my_version <- packageVersion("photobiologyInOut")

## ----example-0-hiden, eval=FALSE, include=TRUE-------
#  # this may be needed in some geographic locations as some Windows TZ strings are
#  # not recognized by all versions of R
#  Sys.setenv(TZ = 'UTC')
#  library(photobiology)
#  library(photobiologyWavebands)
#  library(photobiologyInOut)
#  library(lubridate)
#  library(ggplot2)
#  library(ggmap)
#  library(ggspectra)
#  library(hyperSpec)
#  library(colorSpec)
#  library(pavo)
#  library(readr)

## ----------------------------------------------------
options(tibble.print_max = 5,
        tibble.print_min = 3,
        photobiology.strict.range = NA_integer_)

## ----------------------------------------------------
jazraw.spct <- read_oo_jazdata(file = "data-vignettes/spectrum.jaz")
jazraw.spct <- trim_wl(jazraw.spct, range = c(250, 900))

## ----------------------------------------------------
plot(jazraw.spct)

## ----------------------------------------------------
getWhenMeasured(jazraw.spct)

## ----------------------------------------------------
getInstrDesc(jazraw.spct)

## ----------------------------------------------------
getInstrSettings(jazraw.spct)

## ----------------------------------------------------
jaz.spct <- read_oo_jazirrad(file = "data-vignettes/spectrum.JazIrrad")
jaz0.spct <- jaz.spct
jaz.spct <- trim_wl(jaz.spct, range = c(290, 800))

## ----------------------------------------------------
plot(jaz.spct)

## ----------------------------------------------------
jaz.spct <- fshift(jaz0.spct, range = c(255, 290), f = "mean")
jaz.spct <- trim_wl(jaz.spct, range = c(290, 800))
plot(jaz.spct)

## ----------------------------------------------------
jaz.spct <- smooth_spct(jaz.spct)
plot(jaz.spct)

## ----------------------------------------------------
e_irrad(jaz.spct, PAR())       # W m-2

## ----------------------------------------------------
plot(read_oo_jazirrad(file = "data-vignettes/spectrum.JazIrrad"))

## ----------------------------------------------------
plot(read_oo_jazirrad(file = "data-vignettes/spectrum.JazIrrad"),
     range = c(250,850))

## ----------------------------------------------------
plot(smooth_spct(read_oo_jazirrad(file = "data-vignettes/spectrum.JazIrrad")),
     range = c(250,850))

## ----------------------------------------------------
plot(read_oo_ssirrad(file = "data-vignettes/spectrum.SSIrrad"))

## ----------------------------------------------------
plot(read_avaspec_csv(file = "data-vignettes/spectrum-avaspec.csv"),
     range = c(280, 900), unit.out = "photon")

## ----------------------------------------------------
plot(read_macam_dta(file = "data-vignettes/spectrum.DTA"))

## ----------------------------------------------------
licor.spct <- read_licor_prn(file = "data-vignettes/spectrum.PRN")

## ----------------------------------------------------
licor.spct
cat(comment(licor.spct))
plot(licor.spct, unit.out = "photon")

## ----------------------------------------------------
day.dat <- read_csi_dat(file = "data-vignettes/cr6-day.dat")
day.dat

## ----------------------------------------------------
hour.dat <- read_csi_dat(file = "data-vignettes/cr6-hour.dat")
ggplot(hour.dat, aes(TIMESTAMP, PAR_Den_Avg)) + geom_line()

## ----------------------------------------------------
tuv.spct <- read_tuv_usrout(file = "data-vignettes/usrout.txt",
                            date = ymd("2014-03-21"))
summary(subset(tuv.spct, spct.idx == "A"))
tuv.spct

## ----fig.height=10-----------------------------------
plot(tuv.spct, annotations = c("colour.guide")) +
  facet_wrap(~date, ncol = 2)

## ----------------------------------------------------
tuv.mspct <- subset2mspct(tuv.spct)
tuv.mspct

## ----------------------------------------------------
tuv_nd.spct <- read_tuv_usrout(file = "data-vignettes/usrout.txt")
tuv_nd.spct

## ----------------------------------------------------
lrt.df <- read.table(file = "data-vignettes/libradtran-plain-2col.dat",
                     col.names = c("w.length", "s.e.irrad"))
summary(lrt.df)
libradtran.spct <- source_spct(w.length = lrt.df$w.length,
                               s.e.irrad = lrt.df$s.e.irrad * 1e-3)
plot(libradtran.spct, range = c(250, 2500), unit.out = "photon")

## ----------------------------------------------------
lbr.multi.spct <- read_libradtran_vesa("data-vignettes/libradtran-multi.dat")
print(lbr.multi.spct, n = 5)

## ----------------------------------------------------
z.spct <- read_fmi_cum("data-vignettes/2014-08-21_cum.hel")
class_spct(z.spct)
getWhenMeasured(z.spct)
z.spct

## ----------------------------------------------------
z.mspct <- read_m_fmi_cum(c("data-vignettes/2014-08-21_cum.hel",
                            "data-vignettes/2014-08-22_cum.hel"))
class(z.mspct)
getWhenMeasured(z.mspct)
z.mspct

## ----------------------------------------------------
files <- list.files("./data-vignettes/", "*cum.hel")
files <- paste("./data-vignettes/", files, sep = "")
z1.mspct <- read_m_fmi_cum(files)
class(z1.mspct)
getWhenMeasured(z1.mspct)
z1.mspct

## ----message=FALSE-----------------------------------
z2.mspct <-
  read_m_fmi_cum(files,
                 geocode = geocode("Kumpula, Helsinki, Finland",
                                   source = "google"))
class(z2.mspct)
getWhenMeasured(z2.mspct)
getWhereMeasured(z2.mspct)
z2.mspct

## ----------------------------------------------------
x <- matrix(1:100, ncol = 2)
wl <- 501:550 # in nanometres
z <- mat2mspct(x, wl, "filter_spct", "Tpc")
z

## ----------------------------------------------------
z <- mat2mspct(x, wl, "filter_spct", "Tpc", spct.names = c("A", "B"))
z

## ----------------------------------------------------
xrow <- matrix(1:100, nrow = 2, byrow = TRUE)
z1 <- mat2mspct(xrow, wl, "filter_spct", "Tpc")
z1

## ----------------------------------------------------
z2c.mat <- mspct2mat(z2.mspct, "s.e.irrad")
class(z2c.mat)
dim(z2c.mat)
head(dimnames(z2c.mat)$spct)
head(dimnames(z2c.mat)$w.length)
head(attr(z2c.mat, "w.length"))

## ----------------------------------------------------
z2r.mat <- mspct2mat(z2.mspct, "s.e.irrad", byrow = TRUE)
class(z2r.mat)
dim(z2r.mat)
head(dimnames(z2r.mat)$spct)
head(dimnames(z2r.mat)$w.length)
head(attr(z2r.mat, "w.length"))

## ----------------------------------------------------
z2.hspct <- mspct2hyperSpec(z2.mspct, "s.e.irrad")
class(z2.hspct)
plot(z2.hspct)

## ----------------------------------------------------
class(laser)
laser
plot(laser)

## ----------------------------------------------------
wl(laser) <- list (
  wl = 1e7 / (1/405e-7 - wl (laser)),
  label = expression (lambda / nm)
)
laser
plot(laser)
laser.mspct <-
  hyperSpec2mspct(laser, "source_spct", "s.e.irrad", multiplier = 1e-3)
ggplot(laser.mspct[[1]]) +
  geom_line() +
  stat_peaks(geom = "text", vjust = -1, label.fmt = "%.6g nm", color = "red")

## ----------------------------------------------------
fluorescent.mspct <- colorSpec2mspct(Fs.5nm)
print(fluorescent.mspct, n = 3, n.members = 3)

## ----------------------------------------------------
colorSpec2mspct(Hoya)

## ----------------------------------------------------
fluorescent.spct <- colorSpec2spct(Fs.5nm)
plot(fluorescent.spct) + aes(linetype = spct.idx)

## ----------------------------------------------------
colorSpec2chroma_spct(xyz1931.5nm)

## ----------------------------------------------------
sun.cspec <- spct2colorSpec(sun.spct)
plot(sun.cspec, col = "blue")

## ----------------------------------------------------
spct2colorSpec(yellow_gel.spct)

## ----------------------------------------------------
chroma_spct2colorSpec(beesxyzCMF.spct)

## ----------------------------------------------------
data(sicalis)
class(sicalis)
names(sicalis)

## ----------------------------------------------------
sicalis.mspct <- rspec2mspct(sicalis, "reflector_spct", "Rpc")
summary(sicalis.mspct[[1]])
summary(sicalis.mspct[[2]])
summary(sicalis.mspct[[3]])

## ----------------------------------------------------
ggplot(rbindspct(sicalis.mspct[1:3])) +
  aes(linetype = spct.idx) +
  ylim(0,0.3) +
  geom_line()

## ----------------------------------------------------
print(sicalis.mspct[c(TRUE, FALSE, FALSE)])
ggplot(rbindspct(sicalis.mspct[c(TRUE, FALSE, FALSE)])) +
  aes(linetype = spct.idx) +
  ylim(0,0.15) +
  geom_line() +
  ggtitle("'crown' reflectance spectra")

## ----------------------------------------------------
refl.by.band <- reflectance(sicalis.mspct, w.band = list(Red(), Green(), Blue(), UVA()))
refl.by.band$body.part <- c("crown", "throat", "breast")

## ----------------------------------------------------
refl.red <- reflectance(sicalis.mspct, w.band = Red())
names(refl.red)[2] <- "red.reflectance"
refl.red$body.part <- c("crown", "throat", "breast")
ggplot(refl.red, aes(x = body.part, y = red.reflectance)) +
  stat_summary(fun.data = "mean_se", color = "red") +
  geom_point(alpha = 0.5)

## ----------------------------------------------------
my.locale <- locale(decimal_mark = ",", tz = "EET")
read_oo_jazirrad(file = "data-vignettes/spectrum-comma.JazIrrad",
                 locale = my.locale)

## ----warning=FALSE-----------------------------------
jaz01.spct <- read_oo_jazirrad(file = "data-vignettes/spectrum.JazIrrad",
                               date = NULL)
getWhenMeasured(jaz01.spct)

## ----warning=FALSE-----------------------------------
jaz02.spct <- read_oo_jazirrad(file = "data-vignettes/spectrum.JazIrrad",
                               date = ymd_hms("2015-11-15 12:00:00"))
getWhenMeasured(jaz02.spct)

## ----warning=FALSE-----------------------------------
jaz03.spct <- read_oo_jazirrad(file = "data-vignettes/spectrum.JazIrrad",
                               date = now())
getWhenMeasured(jaz03.spct)

## ----message=FALSE,warning=FALSE---------------------
jaz04.spct <- read_oo_jazirrad(file = "data-vignettes/spectrum.JazIrrad",
                               geocode = geocode("Vikki, 00790 Helsinki, Finland",
                                                 source = "google"))
jaz04.spct
getWhereMeasured(jaz04.spct)

