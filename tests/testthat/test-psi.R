library("photobiology")
library("photobiologyInOut")
library("lubridate")

context("read PSI SpectraPen CSV file)")

test_that("multiple spectra, energy, file top", {

  file.name <- 
    system.file("extdata", "spectrum-psi-spectrapen-SP.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  psi.mspct <- read_spectrapen_csv(file = file.name,
                                  tz = "UTC")
  expect_is(psi.mspct, "source_mspct")
  expect_named(psi.mspct, c("spct.13", "spct.14", "spct.15", "spct.16", "spct.17", "spct.18"))
  expect_equal(as.numeric(getWhenMeasured(psi.mspct[[1]]), tz = "UTC"),
               as.numeric(ymd_hms("2022-10-11 19:08:01 UTC", tz = "UTC"), tz = "UTC"))
  for (spct in psi.mspct) {
    psi.spct <- spct
    expect_equal(nrow(psi.spct), 256)
    expect_equal(ncol(psi.spct), 2)
    expect_equal(psi.spct[1, 1], 327.1)
    expect_equal(psi.spct[256, 1], 793)
    expect_is(psi.spct[[1]], "numeric")
    expect_equal(sum(is.na(psi.spct[[1]])), 0)
    expect_true(all(sign(psi.spct[[1]]) > 0))
    expect_is(psi.spct[[2]], "numeric")
    #  expect_true(all(sign(psi.spct[[2]]) >= 0))
    expect_equal(sum(is.na(psi.spct[[2]])), 0)
    expect_is(psi.spct, "source_spct")
    expect_named(psi.spct, c("w.length", "s.e.irrad"))
    expect_equal(getWhereMeasured(psi.spct), 
                 tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
    expect_gt(length(getWhatMeasured(psi.spct)), 0)
    expect_gt(length(comment(psi.spct)), 0)
  }
})

test_that("multiple spectra, photon, file bottom", {
  
  file.name <- 
    system.file("extdata", "spectrum-psi-spectrapen-SP.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  psi.mspct <- read_spectrapen_csv(file = file.name, tz = "UTC", start.row = 266)
  expect_is(psi.mspct, "source_mspct")
  expect_named(psi.mspct, c("spct.13", "spct.14", "spct.15", "spct.16", "spct.17", "spct.18"))
  expect_equal(as.numeric(getWhenMeasured(psi.mspct[[1]]), tz = "UTC"),
               as.numeric(ymd_hms("2022-10-11 19:08:01 UTC", tz = "UTC"), tz = "UTC"))
  for (spct in psi.mspct) {
    psi.spct <- spct
    expect_equal(nrow(psi.spct), 256)
    expect_equal(ncol(psi.spct), 2)
    expect_equal(psi.spct[1, 1], 327.1)
    expect_equal(psi.spct[256, 1], 793)
    expect_is(psi.spct[[1]], "numeric")
    expect_equal(sum(is.na(psi.spct[[1]])), 0)
    expect_true(all(sign(psi.spct[[1]]) > 0))
    expect_is(psi.spct[[2]], "numeric")
    #  expect_true(all(sign(psi.spct[[2]]) >= 0))
    expect_equal(sum(is.na(psi.spct[[2]])), 0)
    expect_is(psi.spct, "source_spct")
    expect_named(psi.spct, c("w.length", "s.e.irrad"))
    expect_equal(getWhereMeasured(psi.spct), 
                 tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
    expect_gt(length(getWhatMeasured(psi.spct)), 0)
    expect_gt(length(comment(psi.spct)), 0)
  }
})

