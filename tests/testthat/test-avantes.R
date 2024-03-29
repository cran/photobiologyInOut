library("photobiology")
library("photobiologyInOut")
library("lubridate")

context("read Avaspec .csv file)")

test_that("single spectrum (quantum)", {

  file.name <- 
    system.file("extdata", "spectrum-avaspec.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  avaspec.spct <- read_avaspec_csv(file = file.name)
  
  expect_equal(nrow(avaspec.spct), 1604)
  expect_equal(ncol(avaspec.spct), 2)
  expect_equal(avaspec.spct[1, 1],  172.485, tolerance = 0.001)
  expect_equal(avaspec.spct[1604, 1], 1100.222, tolerance = 0.001)
  expect_is(avaspec.spct[[1]], "numeric")
  expect_equal(sum(is.na(avaspec.spct[[1]])), 0)
  expect_true(all(sign(avaspec.spct[[1]]) > 0))
  expect_is(avaspec.spct[[2]], "numeric")
#  expect_true(all(sign(avaspec.spct[[1]]) >= 0))
  expect_equal(sum(is.na(avaspec.spct[[2]])), 0)
  expect_is(avaspec.spct, "source_spct")
  expect_named(avaspec.spct, c("w.length", "s.e.irrad"))
  expect_true(is.na(getWhenMeasured(avaspec.spct)))
  expect_equivalent(getWhereMeasured(avaspec.spct), 
                    data.frame(lon = NA_real_, lat = NA_real_, address = NA_character_, 
                               stringsAsFactors = FALSE))
  expect_gt(length(getWhatMeasured(avaspec.spct)), 0)
  expect_gt(length(comment(avaspec.spct)), 0)
})

