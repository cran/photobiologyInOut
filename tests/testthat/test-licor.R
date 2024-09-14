library("photobiology")
library("photobiologyInOut")
library("lubridate")

context("read LI-1800 PRN file)")

test_that("single spectrum (quantum)", {

  file.name <- 
    system.file("extdata", "spectrum.PRN", 
                package = "photobiologyInOut", mustWork = TRUE)
  licor.spct <- read_licor_prn(file = file.name)
  
  expect_equal(nrow(licor.spct), 601)
  expect_equal(ncol(licor.spct), 2)
  expect_equal(licor.spct[1, 1], 300)
  expect_equal(min(licor.spct), 300)
  expect_equal(licor.spct[601, 1], 900)
  expect_equal(max(licor.spct), 900)
  expect_is(licor.spct[[1]], "numeric")
  expect_equal(sum(is.na(licor.spct[[1]])), 0)
  expect_true(all(sign(licor.spct[[1]]) > 0))
  expect_is(licor.spct[[2]], "numeric")
  expect_true(all(sign(licor.spct[[1]]) >= 0))
  expect_equal(sum(is.na(licor.spct[[2]])), 0)
  expect_is(licor.spct, "source_spct")
  expect_named(licor.spct, c("w.length", "s.q.irrad"))
  when.measured <- when_measured(licor.spct)
  expect_equal(month(when.measured), 8)
  expect_equal(day(when.measured), 23)
  expect_equal(hour(when.measured), 16)
  expect_equal(minute(when.measured), 32)
  expect_equal(second(when.measured), 0)
  expect_equal(getWhereMeasured(licor.spct), 
               tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
  expect_equal(getWhatMeasured(licor.spct), "File: spectrum.PRN")
  #  expect_gt(length(getWhatMeasured(licor.spct)), 0)
  expect_gt(length(comment(licor.spct)), 0)
  
  licor.spct <- 
    read_licor_prn(file = file.name, label = "test-label")
  expect_equal(getWhatMeasured(licor.spct), 
               "File: spectrum.PRN\ntest-label")
  expect_equal(nrow(licor.spct), 601)
  expect_equal(ncol(licor.spct), 2)
  
  file.name <- 
    system.file("extdata", "spectrum-licor-long.PRN", 
                package = "photobiologyInOut", mustWork = TRUE)
  licor.spct <- read_licor_prn(file = file.name)
  
  expect_equal(nrow(licor.spct), 401)
  expect_equal(ncol(licor.spct), 2)
  expect_equal(licor.spct[1, 1], 300)
  expect_equal(min(licor.spct), 300)
  expect_equal(licor.spct[401, 1], 1100)
  expect_equal(max(licor.spct), 1100)
  expect_is(licor.spct[[1]], "numeric")
  expect_equal(sum(is.na(licor.spct[[1]])), 0)
  expect_true(all(sign(licor.spct[[1]]) > 0))
  expect_is(licor.spct[[2]], "numeric")
  expect_equal(what_measured(licor.spct), 
               "File: spectrum-licor-long.PRN")
  
})

test_that("two spectra (quantum)", {
  
  file.name1 <- 
    system.file("extdata", "spectrum.PRN", 
                package = "photobiologyInOut", mustWork = TRUE)
  file.name2 <- 
    system.file("extdata", "spectrum-licor-long.PRN", 
                package = "photobiologyInOut", mustWork = TRUE)
  file.names <- c(file.name1, file.name2)
  licor.mspct <- read_m_licor_prn(files = file.names, tz = NULL)
  
  expect_equal(nrow(licor.mspct[[1]]), 601)
  expect_equal(ncol(licor.mspct[[1]]), 2)
  expect_equal(licor.mspct[[1]][1, 1], 300)
  expect_equal(min(licor.mspct[[1]]), 300)
  expect_equal(licor.mspct[[1]][601, 1], 900)
  expect_equal(max(licor.mspct[[1]]), 900)
  expect_is(licor.mspct[[1]][[1]], "numeric")
  expect_equal(sum(is.na(licor.mspct[[1]][[1]])), 0)
  expect_true(all(sign(licor.mspct[[1]][[1]]) > 0))
  expect_is(licor.mspct[[1]][[2]], "numeric")
  expect_true(all(sign(licor.mspct[[1]][[1]]) >= 0))
  expect_equal(sum(is.na(licor.mspct[[1]][[2]])), 0)
  expect_is(licor.mspct[[1]], "source_spct")
  expect_named(licor.mspct[[1]], c("w.length", "s.q.irrad"))
  when.measured <- when_measured(licor.mspct[[1]])
  expect_equal(month(when.measured), 8)
  expect_equal(day(when.measured), 23)
  expect_equal(hour(when.measured), 16)
  expect_equal(minute(when.measured), 32)
  expect_equal(second(when.measured), 0)
  expect_equal(getWhereMeasured(licor.mspct[[1]]), 
               tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
  expect_equal(getWhatMeasured(licor.mspct[[1]]), "File: spectrum.PRN")
  #  expect_gt(length(getWhatMeasured(licor.mspct[[1]])), 0)
  expect_gt(length(comment(licor.mspct[[1]])), 0)
  
  expect_equal(nrow(licor.mspct[[2]]), 401)
  expect_equal(ncol(licor.mspct[[2]]), 2)
  expect_equal(licor.mspct[[2]][1, 1], 300)
  expect_equal(min(licor.mspct[[2]]), 300)
  expect_equal(licor.mspct[[2]][401, 1], 1100)
  expect_equal(max(licor.mspct[[2]]), 1100)
  expect_is(licor.mspct[[2]][[1]], "numeric")
  expect_equal(sum(is.na(licor.mspct[[2]][[1]])), 0)
  expect_true(all(sign(licor.mspct[[2]][[1]]) > 0))
  expect_is(licor.mspct[[2]][[2]], "numeric")
  expect_equal(what_measured(licor.mspct[[2]]), 
               "File: spectrum-licor-long.PRN")
  
})

test_that("single spectrum Tfr", {
  # for the time being we use a file with reflectance
  file.name <- 
    system.file("extdata", "reflectance.PRN", 
                package = "photobiologyInOut", mustWork = TRUE)
  licor.spct <- read_licor_prn(file = file.name, s.qty = "Tfr")
  
  expect_equal(nrow(licor.spct), 226)
  expect_equal(ncol(licor.spct), 2)
  expect_equal(licor.spct[1, 1], 350)
  expect_equal(min(licor.spct), 350)
  expect_equal(licor.spct[226, 1], 800)
  expect_equal(max(licor.spct), 800)
  expect_is(licor.spct[[1]], "numeric")
  expect_equal(sum(is.na(licor.spct[[1]])), 0)
  expect_true(all(sign(licor.spct[[1]]) > 0))
  expect_is(licor.spct[[2]], "numeric")
  #  expect_true(all(sign(licor.spct[[1]]) >= 0))
  expect_equal(sum(is.na(licor.spct[[2]])), 0)
  expect_is(licor.spct, "filter_spct")
  expect_named(licor.spct, c("w.length", "Tfr"))
  # expect_equal(as.numeric(getWhenMeasured(licor.spct), tz = "UTC"), 
  #              as.numeric(ymd_hms("0000-08-23 14:52:11", tz = "UTC"), tz = "UTC"))
  expect_equal(getWhereMeasured(licor.spct), 
               tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
  expect_gt(length(getWhatMeasured(licor.spct)), 0)
  expect_gt(length(comment(licor.spct)), 0)
  
})

test_that("single spectrum Rfr", {
  
  file.name <- 
    system.file("extdata", "reflectance.PRN", 
                package = "photobiologyInOut", mustWork = TRUE)
  licor.spct <- read_licor_prn(file = file.name, s.qty = "Rfr")
  
  expect_equal(nrow(licor.spct), 226)
  expect_equal(ncol(licor.spct), 2)
  expect_equal(licor.spct[1, 1], 350)
  expect_equal(min(licor.spct), 350)
  expect_equal(licor.spct[226, 1], 800)
  expect_equal(max(licor.spct), 800)
  expect_is(licor.spct[[1]], "numeric")
  expect_equal(sum(is.na(licor.spct[[1]])), 0)
  expect_true(all(sign(licor.spct[[1]]) > 0))
  expect_is(licor.spct[[2]], "numeric")
  #  expect_true(all(sign(licor.spct[[1]]) >= 0))
  expect_equal(sum(is.na(licor.spct[[2]])), 0)
  expect_is(licor.spct, "reflector_spct")
  expect_named(licor.spct, c("w.length", "Rfr"))
  # expect_equal(as.numeric(getWhenMeasured(licor.spct), tz = "UTC"), 
  #              as.numeric(ymd_hms("0000-08-23 14:52:11", tz = "UTC"), tz = "UTC"))
  expect_equal(getWhereMeasured(licor.spct), 
               tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
  expect_gt(length(getWhatMeasured(licor.spct)), 0)
  expect_gt(length(comment(licor.spct)), 0)
  
})

