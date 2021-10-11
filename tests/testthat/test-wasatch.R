library("photobiology")
library("photobiologyInOut")
library("lubridate")

context("read ENLIGHTEN .csv file)")

test_that("single spectrum (energy)", {

  file.name <- 
    system.file("extdata", "enlighten-wasatch-raman.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  wasatch.source.spct <- 
    read_wasatch_csv(file = file.name, s.qty = c(Processed = "s.e.irrad"))
  
  expect_equal(nrow(wasatch.source.spct), 1024)
  expect_equal(ncol(wasatch.source.spct), 2)
  expect_equal(wasatch.source.spct[1, 1],  792.33, tolerance = 0.01)
  expect_equal(wasatch.source.spct[1024, 1], 1067.8, tolerance = 0.01)
  expect_is(wasatch.source.spct[[1]], "numeric")
  expect_equal(sum(is.na(wasatch.source.spct[[1]])), 0)
  expect_true(all(sign(wasatch.source.spct[[1]]) > 0))
  expect_is(wasatch.source.spct[[2]], "numeric")
  expect_equal(sum(is.na(wasatch.source.spct[[2]])), 0)
  expect_is(wasatch.source.spct, "source_spct")
  expect_named(wasatch.source.spct, c("w.length", "s.e.irrad"))
  expect_is(getWhenMeasured(wasatch.source.spct), "POSIXct")
  expect_equivalent(getWhereMeasured(wasatch.source.spct), 
                    data.frame(lon = NA_real_, lat = NA_real_, 
                               address = NA_character_, 
                               stringsAsFactors = FALSE))
  expect_equal(length(getWhatMeasured(wasatch.source.spct)), 1)
  expect_equal(length(comment(wasatch.source.spct)), 1)
})

test_that("single spectrum (raw counts)", {
  
  file.name <- 
    system.file("extdata", "enlighten-wasatch-scope.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  wasatch.raw.spct <- read_wasatch_csv(file = file.name)
  
  expect_equal(nrow(wasatch.raw.spct), 1024)
  expect_equal(ncol(wasatch.raw.spct), 3)
  expect_equal(wasatch.raw.spct[1, 2],  247.94, tolerance = 0.01)
  expect_equal(wasatch.raw.spct[1024, 2], 709.67, tolerance = 0.01)
  expect_is(wasatch.raw.spct[[1]], "integer")
  expect_is(wasatch.raw.spct[[2]], "numeric")
  expect_equal(sum(is.na(wasatch.raw.spct[[1]])), 0)
#  expect_true(all(sign(wasatch.spct[[1]]) > 0))
  expect_is(wasatch.raw.spct[[3]], "numeric")
  expect_equal(sum(is.na(wasatch.raw.spct[[3]])), 0)
  expect_is(wasatch.raw.spct, "raw_spct")
  expect_named(wasatch.raw.spct, c("Pixel", "w.length", "counts"))
  expect_is(getWhenMeasured(wasatch.raw.spct), "POSIXct")
  expect_equivalent(getWhereMeasured(wasatch.raw.spct), 
                    data.frame(lon = NA_real_, lat = NA_real_, 
                               address = NA_character_, 
                               stringsAsFactors = FALSE))
  expect_equal(length(getWhatMeasured(wasatch.raw.spct)), 1)
  expect_equal(length(comment(wasatch.raw.spct)), 1)
})

test_that("single spectrum (absorbance with keep)", {
  
  file.name <- 
    system.file("extdata", "enlighten-wasatch-absorbance.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  wasatch.abs.spct <- read_wasatch_csv(file = file.name)
  
  expect_equal(nrow(wasatch.abs.spct), 1024)
  expect_equal(ncol(wasatch.abs.spct), 6)
  expect_equal(wasatch.abs.spct[1, 2],  247.94, tolerance = 0.01)
  expect_equal(wasatch.abs.spct[1024, 2], 709.67, tolerance = 0.01)
  expect_is(wasatch.abs.spct[[1]], "integer")
  expect_is(wasatch.abs.spct[[2]], "numeric")
  expect_equal(sum(is.na(wasatch.abs.spct[[1]])), 0)
  #  expect_true(all(sign(wasatch.spct[[1]]) > 0))
  expect_is(wasatch.abs.spct[[3]], "numeric")
  expect_is(wasatch.abs.spct[[4]], "numeric")
  expect_is(wasatch.abs.spct[[5]], "numeric")
  expect_is(wasatch.abs.spct[[6]], "numeric")
  expect_equal(sum(is.na(wasatch.abs.spct[[3]])), 0)
  expect_is(wasatch.abs.spct, "filter_spct")
  expect_named(wasatch.abs.spct, c("Pixel", "w.length", "A", "Raw", "Dark", "Reference"))
  expect_is(getWhenMeasured(wasatch.abs.spct), "POSIXct")
  expect_equivalent(getWhereMeasured(wasatch.abs.spct), 
                    data.frame(lon = NA_real_, lat = NA_real_, 
                               address = NA_character_, 
                               stringsAsFactors = FALSE))
  expect_equal(length(getWhatMeasured(wasatch.abs.spct)), 1)
  expect_equal(length(comment(wasatch.abs.spct)), 1)
})

test_that("single spectrum (absorbance with drop)", {
  
  file.name <- 
    system.file("extdata", "enlighten-wasatch-absorbance.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  wasatch.abs1.spct <- read_wasatch_csv(file = file.name, extra.cols = "drop")
  
  expect_equal(nrow(wasatch.abs1.spct), 1024)
  expect_equal(ncol(wasatch.abs1.spct), 2)
  expect_equal(wasatch.abs1.spct[1, 1],  247.94, tolerance = 0.01)
  expect_equal(wasatch.abs1.spct[1024, 1], 709.67, tolerance = 0.01)
  expect_is(wasatch.abs1.spct[[1]], "numeric")
  expect_equal(sum(is.na(wasatch.abs1.spct[[1]])), 0)
  #  expect_true(all(sign(wasatch.spct[[1]]) > 0))
  expect_is(wasatch.abs1.spct[[2]], "numeric")
  expect_equal(sum(is.na(wasatch.abs1.spct[[2]])), 0)
  expect_is(wasatch.abs1.spct, "filter_spct")
  expect_named(wasatch.abs1.spct, c("w.length", "A"))
  expect_is(getWhenMeasured(wasatch.abs1.spct), "POSIXct")
  expect_equivalent(getWhereMeasured(wasatch.abs1.spct), 
                    data.frame(lon = NA_real_, lat = NA_real_, 
                               address = NA_character_, 
                               stringsAsFactors = FALSE))
  expect_equal(length(getWhatMeasured(wasatch.abs1.spct)), 1)
  expect_equal(length(comment(wasatch.abs1.spct)), 1)
})

test_that("single spectrum (absorbance dark)", {
  
  file.name <- 
    system.file("extdata", "enlighten-wasatch-absorbance.csv", 
                package = "photobiologyInOut", mustWork = TRUE)

  wasatch.abs2.spct <- read_wasatch_csv(file = file.name, 
                                        s.qty = c(Dark = "counts"),
                                        extra.cols = "drop")
  
  expect_equal(nrow(wasatch.abs2.spct), 1024)
  expect_equal(ncol(wasatch.abs2.spct), 2)
  expect_equal(wasatch.abs2.spct[1, 1],  247.94, tolerance = 0.01)
  expect_equal(wasatch.abs2.spct[1024, 1], 709.67, tolerance = 0.01)
  expect_is(wasatch.abs2.spct[[1]], "numeric")
  expect_equal(sum(is.na(wasatch.abs2.spct[[1]])), 0)
  #  expect_true(all(sign(wasatch.spct[[1]]) > 0))
  expect_is(wasatch.abs2.spct[[2]], "numeric")
  expect_equal(sum(is.na(wasatch.abs2.spct[[2]])), 0)
  expect_is(wasatch.abs2.spct, "raw_spct")
  expect_named(wasatch.abs2.spct, c("w.length", "counts"))
  expect_is(getWhenMeasured(wasatch.abs2.spct), "POSIXct")
  expect_equivalent(getWhereMeasured(wasatch.abs2.spct), 
                    data.frame(lon = NA_real_, lat = NA_real_, 
                               address = NA_character_, 
                               stringsAsFactors = FALSE))
  expect_equal(length(getWhatMeasured(wasatch.abs2.spct)), 1)
  expect_equal(length(comment(wasatch.abs2.spct)), 1)
})

test_that("single spectrum (transmission with drop)", {
  
  file.name <- 
    system.file("extdata", "enlighten-wasatch-transmission.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  suppressWarnings( # warning caused by negative values in file
    wasatch.tr.spct <- read_wasatch_csv(file = file.name, 
                                      extra.cols = "drop")
  )
  
  expect_equal(nrow(wasatch.tr.spct), 1024)
  expect_equal(ncol(wasatch.tr.spct), 2)
  expect_equal(wasatch.tr.spct[1, 1],  247.94, tolerance = 0.01)
  expect_equal(wasatch.tr.spct[1024, 1], 709.67, tolerance = 0.01)
  expect_is(wasatch.tr.spct[[1]], "numeric")
  expect_equal(sum(is.na(wasatch.tr.spct[[1]])), 0)
  #  expect_true(all(sign(wasatch.spct[[1]]) > 0))
  expect_is(wasatch.tr.spct[[2]], "numeric")
  expect_equal(sum(is.na(wasatch.tr.spct[[2]])), 0)
  expect_is(wasatch.tr.spct, "filter_spct")
  expect_named(wasatch.tr.spct, c("w.length", "Tfr"))
  expect_is(getWhenMeasured(wasatch.tr.spct), "POSIXct")
  expect_equivalent(getWhereMeasured(wasatch.tr.spct), 
                    data.frame(lon = NA_real_, lat = NA_real_, 
                               address = NA_character_, 
                               stringsAsFactors = FALSE))
  expect_equal(length(getWhatMeasured(wasatch.tr.spct)), 1)
  expect_equal(length(comment(wasatch.tr.spct)), 1)
})

test_that("single spectrum (transmission with split)", {
  
  file.name <- 
    system.file("extdata", "enlighten-wasatch-transmission.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  suppressWarnings( # warning caused by negative values in file
    wasatch.tr.mspct <- read_wasatch_csv(file = file.name, 
                                      extra.cols = "split")
  )
  
  expect_equal(length(wasatch.tr.mspct), 4)
  expect_is(wasatch.tr.mspct, "generic_mspct")
  expect_is(wasatch.tr.mspct[["Processed"]], "filter_spct")
  expect_is(wasatch.tr.mspct[["Raw"]], "raw_spct")
  expect_is(wasatch.tr.mspct[["Dark"]], "raw_spct")
  expect_is(wasatch.tr.mspct[["Reference"]], "raw_spct")
})
