library("photobiology")
library("photobiologyInOut")
library("lubridate")

context("read CID CI-710s (.CSV file)")

test_that("Absorbance Mode", {

  file.name <- 
    system.file("extdata", "cid-spectravue-Abs-Measurements.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  cid_Abs.spct <- read_cid_spectravue_csv(file = file.name)
  
  expect_is(cid_Abs.spct, "object_spct")
  expect_equal(nrow(cid_Abs.spct), 1235)
  expect_equal(ncol(cid_Abs.spct), 3)
  expect_equal(cid_Abs.spct[1, 1], 380.3094, tolerance = 1e-4)
  expect_equal(wl_min(cid_Abs.spct), 380.3094, tolerance = 1e-4)
  expect_equal(cid_Abs.spct[1235, 1], 1099.885, tolerance = 1e-4)
  expect_equal(wl_max(cid_Abs.spct), 1099.885, tolerance = 1e-4)
  expect_is(cid_Abs.spct[[1]], "numeric")
  expect_equal(sum(is.na(cid_Abs.spct[[1]])), 0)
  expect_true(all(sign(cid_Abs.spct[[1]]) > 0))
  expect_is(cid_Abs.spct[[2]], "numeric")
  expect_true(all(sign(cid_Abs.spct[[2]]) >= 0))
  expect_equal(sum(is.na(cid_Abs.spct[[2]])), 0)
  expect_is(cid_Abs.spct[[3]], "numeric")
  expect_true(all(sign(cid_Abs.spct[[3]]) >= 0))
  expect_equal(sum(is.na(cid_Abs.spct[[3]])), 0)
  expect_named(cid_Abs.spct, c("w.length", "Rfr", "Tfr"))
  expect_equal(as.numeric(getWhenMeasured(cid_Abs.spct), tz = "UTC"), 
               as.numeric(ymd_hms("2022-03-11 17:18:08 UTC", tz = "UTC"), tz = "UTC"))
  expect_equal(getWhereMeasured(cid_Abs.spct), 
               tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
  expect_gt(length(getWhatMeasured(cid_Abs.spct)), 0)
  expect_gt(length(comment(cid_Abs.spct)), 0)
  
  file.name <- 
    system.file("extdata", "cid-spectravue-Abs-Measurements.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  cid_A.spct <- read_cid_spectravue_csv(file = file.name, absorbance.to = "A")
  
  expect_is(cid_A.spct, "filter_spct")
  expect_equal(nrow(cid_A.spct), 1235)
  expect_equal(ncol(cid_A.spct), 2)
  expect_equal(cid_A.spct[1, 1], 380.3094, tolerance = 1e-4)
  expect_equal(wl_min(cid_A.spct), 380.3094, tolerance = 1e-4)
  expect_equal(cid_A.spct[1235, 1], 1099.885, tolerance = 1e-4)
  expect_equal(wl_max(cid_A.spct), 1099.885, tolerance = 1e-4)
  expect_is(cid_A.spct[[1]], "numeric")
  expect_equal(sum(is.na(cid_A.spct[[1]])), 0)
  expect_true(all(sign(cid_A.spct[[1]]) > 0))
  expect_is(cid_A.spct[[2]], "numeric")
  expect_true(all(sign(cid_A.spct[[2]]) >= 0))
  expect_equal(sum(is.na(cid_A.spct[[2]])), 0)
  expect_named(cid_A.spct, c("w.length", "A"))
  expect_equal(as.numeric(getWhenMeasured(cid_Abs.spct), tz = "UTC"), 
               as.numeric(ymd_hms("2022-03-11 17:18:08 UTC", tz = "UTC"), tz = "UTC"))
  expect_equal(getWhereMeasured(cid_Abs.spct), 
               tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
  expect_gt(length(getWhatMeasured(cid_A.spct)), 0)
  expect_gt(length(comment(cid_A.spct)), 0)
  
})

test_that("Reflectance Mode", {

  file.name <- 
    system.file("extdata", "cid-spectravue-Rpc-Measurements.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  cid_Rpc.spct <- read_cid_spectravue_csv(file = file.name)

  expect_is(cid_Rpc.spct, "reflector_spct")
  expect_equal(nrow(cid_Rpc.spct), 1235)
  expect_equal(ncol(cid_Rpc.spct), 2)
  expect_equal(cid_Rpc.spct[1, 1], 380.3094, tolerance = 1e-4)
  expect_equal(wl_min(cid_Rpc.spct), 380.3094, tolerance = 1e-4)
  expect_equal(cid_Rpc.spct[1235, 1], 1099.885, tolerance = 1e-4)
  expect_equal(wl_max(cid_Rpc.spct), 1099.885, tolerance = 1e-4)
  expect_is(cid_Rpc.spct[[1]], "numeric")
  expect_equal(sum(is.na(cid_Rpc.spct[[1]])), 0)
  expect_true(all(sign(cid_Rpc.spct[[1]]) > 0))
  expect_is(cid_Rpc.spct[[2]], "numeric")
  expect_true(all(sign(cid_Rpc.spct[[2]]) >= 0))
  expect_equal(sum(is.na(cid_Rpc.spct[[2]])), 0)
  expect_named(cid_Rpc.spct, c("w.length", "Rfr"))
  expect_equal(as.numeric(getWhenMeasured(cid_Rpc.spct), tz = "UTC"), 
               as.numeric(ymd_hms("2022-03-08 15:02:47 UTC", tz = "UTC"), tz = "UTC"))
  expect_equal(getWhereMeasured(cid_Rpc.spct), 
               tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
  expect_gt(length(getWhatMeasured(cid_Rpc.spct)), 0)
  expect_gt(length(comment(cid_Rpc.spct)), 0)
  
})

test_that("Multiple Modes", {
  
  file.name <- 
    system.file("extdata", "cid-spectravue-multi-Measurements.csv", 
                package = "photobiologyInOut", mustWork = TRUE)
  cid_multi.mspct <- read_cid_spectravue_csv(file = file.name)
  
  expect_is(cid_multi.mspct, "generic_mspct")
  cid_multi1.spct <- cid_multi.mspct[[1]]
  cid_multi2.spct <- cid_multi.mspct[[2]]
  
  expect_is(cid_multi1.spct, "filter_spct")
  expect_is(cid_multi2.spct, "reflector_spct")
  
  expect_equal(nrow(cid_multi1.spct), 1235)
  expect_equal(nrow(cid_multi2.spct), 1235)
  expect_equal(ncol(cid_multi1.spct), 2)
  expect_equal(ncol(cid_multi2.spct), 2)
  expect_equal(cid_multi1.spct[1, 1], 380.3094, tolerance = 1e-4)
  expect_equal(cid_multi2.spct[1, 1], 380.3094, tolerance = 1e-4)
  expect_equal(wl_min(cid_multi1.spct), 380.3094, tolerance = 1e-4)
  expect_equal(wl_min(cid_multi2.spct), 380.3094, tolerance = 1e-4)
  expect_equal(cid_multi1.spct[1235, 1], 1099.885, tolerance = 1e-4)
  expect_equal(cid_multi2.spct[1235, 1], 1099.885, tolerance = 1e-4)
  expect_equal(wl_max(cid_multi1.spct), 1099.885, tolerance = 1e-4)
  expect_equal(wl_max(cid_multi2.spct), 1099.885, tolerance = 1e-4)
  expect_is(cid_multi1.spct[[1]], "numeric")
  expect_is(cid_multi2.spct[[1]], "numeric")
  expect_equal(sum(is.na(cid_multi1.spct[[1]])), 0)
  expect_equal(sum(is.na(cid_multi2.spct[[1]])), 0)
  expect_true(all(sign(cid_multi1.spct[[1]]) > 0))
  expect_true(all(sign(cid_multi2.spct[[1]]) > 0))
  expect_is(cid_multi1.spct[[2]], "numeric")
  expect_is(cid_multi2.spct[[2]], "numeric")
  expect_true(all(sign(cid_multi1.spct[[2]]) >= 0))
  expect_true(all(sign(cid_multi2.spct[[2]]) >= 0))
  expect_equal(sum(is.na(cid_multi1.spct[[2]])), 0)
  expect_equal(sum(is.na(cid_multi2.spct[[2]])), 0)
  expect_named(cid_multi1.spct, c("w.length", "Tfr"))
  expect_named(cid_multi2.spct, c("w.length", "Rfr"))
  expect_equal(as.numeric(getWhenMeasured(cid_multi1.spct), tz = "UTC"), 
               as.numeric(ymd_hms("2022-03-11 17:19:21 UTC", tz = "UTC"), tz = "UTC"))
  expect_equal(as.numeric(getWhenMeasured(cid_multi2.spct), tz = "UTC"), 
               as.numeric(ymd_hms("2022-03-11 17:19:21 UTC", tz = "UTC"), tz = "UTC"))
  expect_equal(getWhereMeasured(cid_multi1.spct), 
               tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
  expect_equal(getWhereMeasured(cid_multi2.spct), 
               tibble::tibble(lon = NA_real_, lat = NA_real_, address = NA_character_))
  expect_gt(length(getWhatMeasured(cid_multi1.spct)), 0)
  expect_gt(length(getWhatMeasured(cid_multi2.spct)), 0)
  expect_gt(length(comment(cid_multi1.spct)), 0)
  expect_gt(length(comment(cid_multi1.spct)), 0)
  
})

