context("Parser")

test_that("parser works", {
  markup <- "
  > Neurons
  expressed: SYT1
  reference: this should not be parsed

  > N2
  expressed uniq: MBP

  > N3
  expressed: OOO
  not expressed: XXX
  "

  marker.info <- parseMarkerFile(markup, is.text=T)

  # General parsing
  expect_true(all(names(marker.info) == c("Neurons", "N2", "N3")))
  expect_equal(length(marker.info), 3)
  expect_equal(length(marker.info[[1]]), 3)

  # Parent
  for (n in c("Neurons", "N2", "N3")) {
    expect_equal(marker.info[[n]]$parent, "root")
  }

  # Expressed
  expect_equal(marker.info[[1]]$expressed, "SYT1")
  expect_equal(marker.info[[2]]$expressed, "MBP")
  expect_equal(marker.info[[3]]$expressed, "OOO")

  # Not expressed
  expect_equal(length(marker.info[[1]]$not_expressed), 1)
  expect_equal(length(marker.info[[2]]$not_expressed), 0)
  expect_equal(length(marker.info[[3]]$not_expressed), 2)

  expect_equal(marker.info[[1]]$not_expressed, "MBP")
  expect_true(all(c("MBP", "XXX") %in% marker.info[[3]]$not_expressed))
})

test_that("expressed field is presented", {
  markup <- "
  > Neurons
  expressed: SYT1

  > N2
  test: MBP
  "

  expect_error(parseMarkerFile(markup, is.text=T))

  markup <- "
  > Neurons
  expressed: SYT1

  > N2
  expressed:
  "

  expect_error(parseMarkerFile(markup, is.text=T))

  markup <- "
  > Neurons
  expressed: SYT1

  > N2
  expressed uniq: SSS
  "

  parseMarkerFile(markup, is.text=T)
})
