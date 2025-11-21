#' Draw from a fitted PUcopula
#' @param model A PUcopula model from [fitPUcopula()]
#' @param n Number of rows
#' @return Matrix of U(0,1) draws
#' @export
simulateCopula <- function(model, n) {
  model@rand(n)
}

#' Backtransform copula samples using fitted marginals
#'
#' Takes a matrix of copula samples (on [0,1]) and a list of marginal models
#' and returns a data.frame on the original marginal scales.
#'
#' @param u Numeric matrix of size n x p with copula samples
#' @param marginals List of marginals from [estimateMarginals()]
#'
#' @return data.frame with n rows and length(marginals) columns
#' @keywords internal
backtransform_marginals <- function(u, marginals) {
  n <- nrow(u)
  df <- as.data.frame(matrix(nrow = n, ncol = length(marginals)))
  names(df) <- names(marginals)

  for (j in seq_along(marginals)) {
    m <- marginals[[j]]

    if (inherits(m, "ecdf")) {
      xs <- environment(m)$x
      ps <- m(xs)
      df[[j]] <- stats::approx(ps, xs, xout = u[, j], ties = "ordered")$y

    } else if (inherits(m, "logspline")) {
      df[[j]] <- logspline::qlogspline(u[, j], m)

    } else if (inherits(m, "oldlogspline")) {
      df[[j]] <- logspline::qoldlogspline(u[, j], m)

    } else if (inherits(m, "table")) {
      vbreaks <- cumsum(m)
      if (length(vbreaks) <= 2) {
        df[[j]] <- (findInterval(u[, j], vbreaks) + 1) |> as.factor()
        levels(df[[j]]) <- names(m)
      } else {
        df[[j]] <- transform_u(u[, j], as.numeric(names(vbreaks)), vbreaks)
      }

    } else {
      df[[j]] <- u[, j]
    }
  }

  df
}

#' Generate synthetic data
#'
#' Combines a fitted PUcopula and marginal models to produce a synthetic
#' data.frame. Optionally restores original factor structure, names and classes.
#'
#' @param n Integer, number of rows to generate
#' @param copula A PUcopula model
#' @param marginals List of marginals from [estimateMarginals()]
#' @param original_levels Optional `preprocessData()$original_levels`
#' @param original_varnames Optional vector of original column names
#' @param original_classes Optional named vector of original classes
#' @return Synthetic data.frame
#' @export
#' @importFrom logspline qlogspline qoldlogspline
#' @importFrom stats approx
generateSynthetic <- function(
  n,
  copula,
  marginals,
  original_levels = NULL,
  original_varnames = NULL,
  original_classes = NULL
) {
  if (is.character(n)) stop("'n' must be numeric.")
  if (n == 0) return(NULL)

  # 1. simulate from copula
  u <- simulateCopula(copula, n)

  # df <- as.data.frame(matrix(nrow = n, ncol = length(marginals)))
  # names(df) <- names(marginals)
  #
  # for (j in seq_along(marginals)) {
  #   m <- marginals[[j]]
  #   if (inherits(m, "ecdf")) {
  #     xs <- environment(m)$x
  #     ps <- m(xs)
  #     df[[j]] <- stats::approx(ps, xs, xout = u[, j], ties = "ordered")$y
  #  # } else if (inherits(m$qfun, "logspline")) { #???????????????????
  #   } else if (inherits(m, "logspline")) { #!!!!!!!!!!!!!!!!!!!!!
  #    # df[[j]] <- logspline::qlogspline(u[, j], m$qfun)#???
  #     df[[j]] <- logspline::qlogspline(u[, j], m)#!!
  #  # } else if (inherits(m$qfun, "oldlogspline")) { #????????????????????????
  #   } else if (inherits(m, "oldlogspline")) { #!!!!!!!!!!!!!!!!!!!!!!!!!!
  #     ##df[[j]] <- logspline::qoldlogspline(u[, j], m$qfun) #??
  #     df[[j]] <- logspline::qoldlogspline(u[, j], m) #!!
  #  # } else if (inherits(m$qfun, "table")) { #??????????????????????
  #   } else if (inherits(m, "table")) { #!!!!!!!!!!!!!!!!!!!
  #     #vbreaks <- cumsum(m$qfun) #??
  #     vbreaks <- cumsum(m) #!!
  #     if (length(vbreaks) <= 2) {
  #       df[[j]] <- (findInterval(u[, j], vbreaks) + 1) |> as.factor()
  #       #levels(df[[j]]) <- names(m$qfun) ##???
  #       levels(df[[j]]) <- names(m) #!!
  #     } else {
  #       df[[j]] <- transform_u(u[, j], as.numeric(names(vbreaks)), vbreaks)
  #     }
  #   } else {
  #     df[[j]] <- u[, j]
  #   }
  # }

  # 2. backtransform marginals
  df <- backtransform_marginals(u, marginals)

  # 3. restore original structure
  if (!is.null(original_levels)) {
    df_raw <- df
    df <- postprocessData(df_raw, original_levels$dummies)
  }
  if (!is.null(original_varnames)) {
    df <- df[original_varnames]
  }
  if (!is.null(original_classes)) {
    int_cols <- names(original_classes)[original_classes == "integer"]
    if (length(int_cols)) {
      df[int_cols] <- lapply(df[int_cols], function(x) as.integer(round(x)))
    }
  }
  df
}
