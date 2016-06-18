library("photobiology")
library("photobiologyInOut")
library("lubridate")
library("readr")

context("read Ocean Optics")

test_that("jazz", {

  # warnings triggered by negative irradiance values in file
  suppressWarnings(jaz.spct <- 
                     read_oo_jazirrad(file = "data-test/spectrum.JazIrrad",
                                      tz = "EET"))
  
  expect_equal(nrow(jaz.spct), 2048)
  expect_equal(ncol(jaz.spct), 2)
  expect_equal(jaz.spct[1, 1], 188.8252, tolerance = 0.0001)
  expect_equal(jaz.spct[2048, 1], 1033.148, tolerance = 0.0001)
  expect_is(jaz.spct[[1]], "numeric")
  expect_equal(sum(is.na(jaz.spct[[1]])), 0)
  expect_true(all(sign(jaz.spct[[1]]) > 0))
  expect_is(jaz.spct[[2]], "numeric")
  expect_equal(sum(is.na(jaz.spct[[2]])), 0)
  expect_is(jaz.spct, "source_spct")
  expect_named(jaz.spct, c("w.length", "s.e.irrad"))
  expect_equal(as.numeric(getWhenMeasured(jaz.spct), tz = "EET"),
               as.numeric(ymd_hms("2015-02-03 09:44:41", tz = "EET"), tz = "EET"))
  expect_equal(getWhereMeasured(jaz.spct), 
               data.frame(lon = NA_real_, lat = NA_real_))
  expect_equal(getWhatMeasured(jaz.spct), "File: data-test/spectrum.JazIrrad")
  expect_equal(getTimeUnit(jaz.spct), "second")
  expect_gt(length(comment(jaz.spct)), 0)
})

test_that("jazz_raw", {
  
  jaz.spct <- read_oo_jazdata(file = "data-test/spectrum.jaz", 
                              tz = "EET")
  
  expect_equal(nrow(jaz.spct), 2048)
  expect_equal(ncol(jaz.spct), 2)
  expect_equal(jaz.spct[1, 1], 190.313904, tolerance = 0.000001)
  expect_equal(jaz.spct[2048, 1], 892.611511, tolerance = 0.000001)
  expect_is(jaz.spct[[1]], "numeric")
  expect_equal(sum(is.na(jaz.spct[[1]])), 0)
  expect_true(all(sign(jaz.spct[[1]]) > 0))
  expect_is(jaz.spct[[2]], "numeric")
  expect_equal(sum(is.na(jaz.spct[[2]])), 0)
  expect_is(jaz.spct, "raw_spct")
  expect_named(jaz.spct, c("w.length", "counts"))
  expect_equal(as.numeric(getWhenMeasured(jaz.spct), tz = "EET"),
               as.numeric(ymd_hms("2016-04-25 12:49:02", tz = "EET"), tz = "EET"))
  expect_equal(getWhereMeasured(jaz.spct), 
               data.frame(lon = NA_real_, lat = NA_real_))
  expect_equal(getWhatMeasured(jaz.spct), "File: data-test/spectrum.jaz")
  expect_gt(length(comment(jaz.spct)), 0)
})


test_that("SpectraSuite", {
  
  ss.spct <- read_oo_ssirrad(file = "data-test/spectrum.SSIrrad", 
                             tz = "CET")
  
  expect_equal(nrow(ss.spct), 1044)
  expect_equal(ncol(ss.spct), 2)
  expect_equal(ss.spct[1, 1], 199.08)
  expect_equal(ss.spct[1044, 1], 998.61)
  expect_is(ss.spct[[1]], "numeric")
  expect_equal(sum(is.na(ss.spct[[1]])), 0)
  expect_true(all(sign(ss.spct[[1]]) > 0))
  expect_is(ss.spct[[2]], "numeric")
  expect_equal(sum(is.na(ss.spct[[2]])), 0)
  expect_is(ss.spct, "source_spct")
  expect_named(ss.spct, c("w.length", "s.e.irrad"))
  expect_equal(as.numeric(getWhenMeasured(ss.spct), tz = "CET"),
               as.numeric(ymd_hms("2013-05-06 15:13:40", tz = "CET"), tz = "CET"))
  expect_equal(getWhereMeasured(ss.spct), 
               data.frame(lon = NA_real_, lat = NA_real_))
  expect_equal(getWhatMeasured(ss.spct), NA)
  expect_equal(getTimeUnit(ss.spct), "second")
  expect_gt(length(comment(ss.spct)), 0)
})

test_that("jazz-comma", {
  
  my.locale <- readr::locale("en", decimal_mark = ",", 
                             tz = "EET")
  # warnings triggered by negative irradiance values in file
  suppressWarnings(jaz.spct <- 
                     read_oo_jazirrad(file = "data-test/spectrum-comma.JazIrrad", 
                                      locale = my.locale))
  
  expect_equal(nrow(jaz.spct), 2048)
  expect_equal(ncol(jaz.spct), 2)
  expect_equal(jaz.spct[1, 1], 188.8252, tolerance = 0.0001)
  expect_equal(jaz.spct[2048, 1], 1033.148, tolerance = 0.0001)
  expect_is(jaz.spct[[1]], "numeric")
  expect_equal(sum(is.na(jaz.spct[[1]])), 0)
  expect_true(all(sign(jaz.spct[[1]]) > 0))
  expect_is(jaz.spct[[2]], "numeric")
  expect_equal(sum(is.na(jaz.spct[[2]])), 0)
  expect_is(jaz.spct, "source_spct")
  expect_named(jaz.spct, c("w.length", "s.e.irrad"))
  expect_equal(as.numeric(getWhenMeasured(jaz.spct), tz = "EET"),
               as.numeric(ymd_hms("2015-02-03 09:44:41", tz = "EET"), tz = "EET"))
  expect_equal(getWhereMeasured(jaz.spct), 
               data.frame(lon = NA_real_, lat = NA_real_))
  expect_equal(getWhatMeasured(jaz.spct), "File: data-test/spectrum-comma.JazIrrad")
  expect_equal(getTimeUnit(jaz.spct), "second")
  expect_gt(length(comment(jaz.spct)), 0)
})

test_that("jazz_raw", {
  
  my.locale <- readr::locale("en", decimal_mark = ",", tz = "EET")
  jaz.spct <- read_oo_jazdata(file = "data-test/spectrum-comma.jaz", 
                              locale = my.locale)
  
  expect_equal(nrow(jaz.spct), 2048)
  expect_equal(ncol(jaz.spct), 2)
  expect_equal(jaz.spct[1, 1], 190.313904, tolerance = 0.000001)
  expect_equal(jaz.spct[2048, 1], 892.611511, tolerance = 0.000001)
  expect_is(jaz.spct[[1]], "numeric")
  expect_equal(sum(is.na(jaz.spct[[1]])), 0)
  expect_true(all(sign(jaz.spct[[1]]) > 0))
  expect_is(jaz.spct[[2]], "numeric")
  expect_equal(sum(is.na(jaz.spct[[2]])), 0)
  expect_is(jaz.spct, "raw_spct")
  expect_named(jaz.spct, c("w.length", "counts"))
  expect_equal(as.numeric(getWhenMeasured(jaz.spct), tz = "EET"),
               as.numeric(ymd_hms("2016-04-25 12:49:02", tz = "EET"), tz = "EET"))
  expect_equal(getWhereMeasured(jaz.spct), 
               data.frame(lon = NA_real_, lat = NA_real_))
  expect_equal(getWhatMeasured(jaz.spct), "File: data-test/spectrum-comma.jaz")
  expect_gt(length(comment(jaz.spct)), 0)
})


test_that("SpectraSuite", {
  
  my.locale <- readr::locale("en", decimal_mark = ",", tz = "CET")

  ss.spct <- read_oo_ssirrad(file = "data-test/spectrum-comma.SSIrrad", 
                             locale = my.locale)
  
  expect_equal(nrow(ss.spct), 1044)
  expect_equal(ncol(ss.spct), 2)
  expect_equal(ss.spct[1, 1], 199.08)
  expect_equal(ss.spct[1044, 1], 998.61)
  expect_is(ss.spct[[1]], "numeric")
  expect_equal(sum(is.na(ss.spct[[1]])), 0)
  expect_true(all(sign(ss.spct[[1]]) > 0))
  expect_is(ss.spct[[2]], "numeric")
  expect_equal(sum(is.na(ss.spct[[2]])), 0)
  expect_is(ss.spct, "source_spct")
  expect_named(ss.spct, c("w.length", "s.e.irrad"))
  expect_equal(as.numeric(getWhenMeasured(ss.spct), tz = "CET"), 
               as.numeric(ymd_hms("2013-05-06 15:13:40", tz = "CET"), tz = "CET"))
  expect_equal(getWhereMeasured(ss.spct), 
               data.frame(lon = NA_real_, lat = NA_real_))
  expect_equal(getWhatMeasured(ss.spct), NA)
  expect_equal(getTimeUnit(ss.spct), "second")
  expect_gt(length(comment(ss.spct)), 0)
})

test_that("pi_raw", {
  
  my.date <- now(tz = "EET")
  pi.spct <- read_oo_pidata(file = "data-test/spectrum.pi", 
                            date = my.date,
                            npixels = 2048)
  
  expect_equal(nrow(pi.spct), 2048)
  expect_equal(ncol(pi.spct), 2)
  expect_equal(pi.spct[1, 1], 188.41408000, tolerance = 0.000001)
  expect_equal(pi.spct[2048, 1], 1035.61297366, tolerance = 0.000001)
  expect_is(pi.spct[[1]], "numeric")
  expect_equal(sum(is.na(pi.spct[[1]])), 0)
  expect_true(all(sign(pi.spct[[1]]) > 0))
  expect_is(pi.spct[[2]], "numeric")
  expect_equal(sum(is.na(pi.spct[[2]])), 0)
  expect_is(pi.spct, "raw_spct")
  expect_named(pi.spct, c("w.length", "counts"))
  expect_equal(as.numeric(getWhenMeasured(pi.spct), tz = "EET"),
               as.numeric(my.date, tz = "EET"))
  expect_equal(getWhereMeasured(pi.spct), 
               data.frame(lon = NA_real_, lat = NA_real_))
  expect_equal(getWhatMeasured(pi.spct), "File: data-test/spectrum.pi")
  expect_gt(length(comment(pi.spct)), 0)
})
