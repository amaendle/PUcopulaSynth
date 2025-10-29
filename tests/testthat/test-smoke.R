test_that("pipeline runs end-to-end on a tiny toy dataset", {
  skip_if_not_installed("PUcopula")
  set.seed(1)

  toy <- data.frame(
    x = rnorm(50),
    y = rpois(50, 2),
    z = factor(sample(c("a","b","c"), 50, TRUE)),
    w = ordered(sample(1:4, 50, TRUE))
  )

  pre  <- preprocessData(toy)
  cop  <- fitPUcopula(pre$data, driver_strength_factor = 0.3, bin_size = 3, family = "binom")
  marg <- estimateMarginals(pre$data, k = 3)

  syn  <- generateSynthetic(25, cop, marg,
                            original_levels = pre$original_levels,
                            original_varnames = names(toy),
                            original_classes  = sapply(toy, class))
  expect_s3_class(syn, "data.frame")
  expect_equal(nrow(syn), 25)
  expect_setequal(names(syn), names(toy))
})
