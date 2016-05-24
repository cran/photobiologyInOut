library("photobiology")
library("photobiologyInOut")
library("lubridate")

context("read LI-1800 PRN file)")

test_that("single spectrum (quantum)", {

  licor.spct <- read_licor_prn(file = "data-test/spectrum.PRN")
  
  expect_equal(nrow(licor.spct), 601)
  expect_equal(ncol(licor.spct), 2)
  expect_equal(licor.spct[1, 1], 300)
  expect_equal(licor.spct[601, 1], 900)
  expect_is(licor.spct[[1]], "numeric")
  expect_equal(sum(is.na(licor.spct[[1]])), 0)
  expect_true(all(sign(licor.spct[[1]]) > 0))
  expect_is(licor.spct[[2]], "numeric")
#  expect_true(all(sign(licor.spct[[1]]) >= 0))
  expect_equal(sum(is.na(licor.spct[[2]])), 0)
  expect_is(licor.spct, "source_spct")
  expect_named(licor.spct, c("w.length", "s.q.irrad"))
  # expect_equal(as.numeric(getWhenMeasured(licor.spct), tz = "UTC"), 
  #              as.numeric(ymd_hms("0000-08-23 14:52:11", tz = "UTC"), tz = "UTC"))
  expect_equal(getWhereMeasured(licor.spct), 
               data.frame(lon = NA_real_, lat = NA_real_))
  expect_gt(length(getWhatMeasured(licor.spct)), 0)
  expect_gt(length(comment(licor.spct)), 0)
})

