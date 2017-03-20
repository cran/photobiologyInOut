library("photobiology")
library("photobiologyInOut")
library("lubridate")

context("read FMI (Anders)")

test_that("single spectrum", {

  file.name <- 
    system.file("extdata", "2014-08-21_cum.hel", 
                package = "photobiologyInOut", mustWork = TRUE)
  fmi.spct <- read_fmi_cum(file = file.name)
  
  expect_equal(nrow(fmi.spct), 511)
  expect_equal(ncol(fmi.spct), 2)
  expect_equal(fmi.spct[1, 1], 290)
  expect_equal(fmi.spct[511, 1], 800)
  expect_is(fmi.spct[[1]], "numeric")
  expect_equal(sum(is.na(fmi.spct[[1]])), 0)
  expect_true(all(sign(fmi.spct[[1]]) > 0))
  expect_is(fmi.spct[[2]], "numeric")
  expect_true(all(sign(fmi.spct[[1]]) >= 0))
  expect_equal(sum(is.na(fmi.spct[[2]])), 0)
  expect_is(fmi.spct, "source_spct")
  expect_named(fmi.spct, c("w.length", "s.e.irrad"))
  # expect_equal(getWhenMeasured(fmi.spct), 
  #              ymd("2014-08-21", tz = "UTC"))
#  expect_gt(length(comment(fmi.spct)), 0)
})

test_that("multiple spectrum", {
  
  my.files <- 
    system.file("extdata", c("2014-08-21_cum.hel", "2014-08-22_cum.hel"), 
                package = "photobiologyInOut", mustWork = TRUE)
  fmi.mspct <- read_m_fmi_cum(files = my.files)
  
  expect_is(fmi.mspct, "source_mspct")
  expect_equal(length(fmi.mspct), 2)
  expect_equal(dim(fmi.mspct), c(2, 1))
  expect_false(attr(fmi.mspct, "mspct.byrow", exact = TRUE))
  expect_named(fmi.mspct, c( "2014_08_21_cum.hel", "2014_08_22_cum.hel"))
  expect_equal(nrow(fmi.mspct[[1]]), 511)
  expect_equal(ncol(fmi.mspct[[1]]), 2)
  expect_equal(fmi.mspct[[1]][1, 1], 290)
  expect_equal(fmi.mspct[[1]][511, 1], 800)
  expect_is(fmi.mspct[[1]][[1]], "numeric")
  expect_equal(sum(is.na(fmi.mspct[[1]][[1]])), 0)
  expect_true(all(sign(fmi.mspct[[1]][[1]]) > 0))
  expect_is(fmi.mspct[[1]][[2]], "numeric")
  expect_true(all(sign(fmi.mspct[[1]][[1]]) >= 0))
  expect_equal(sum(is.na(fmi.mspct[[1]][[2]])), 0)
  expect_is(fmi.mspct[[1]], "source_spct")
  expect_named(fmi.mspct[[1]], c("w.length", "s.e.irrad"))
  # expect_equal(getWhenMeasured(fmi.mspct[[1]]), 
  #              ymd("2014-08-21", tz = "UTC"))
  #  expect_gt(length(comment(fmi.spct)), 0)
})

