context("Annotation")

test_that("normalization works", {
  cm <- matrix(runif(200), nrow=25)
  cm.tfidf <- normalizeTfIdfWithFeatures(cm)

  expect_true(all(dim(cm.tfidf) == dim(cm)))

  cm.tc <- normalizeTCWithFeatures(cm)
  expect_true(all(dim(cm.tc) == dim(cm)))
  expect_true(all(abs(apply(cm.tc, 1, max) - 1) < 1e-10))
})
