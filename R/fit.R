#' Fit a PUcopula model
#'
#' Fits a PUcopula model on a (preprocessed) data matrix with optional
#' rank-binning and jitter for numeric variables.
#'
#' @param data Preprocessed data.frame (e.g., `preprocessData()$data`)
#' @param driver_strength_factor Numeric scalar or vector in (0,1] used to scale rows per variable
#' @param bin_size Numeric scalar, vector, or named list with bin sizes
#' @param jitter FALSE, numeric (single) or named list mapping variables to jitter factors
#' @param family PUcopula family, e.g. "binom" or "nbinom"
#' @return A `PUcopula::PUCopula` model
#' @export
#' @import PUcopula
fitPUcopula <- function(
  data, driver_strength_factor = 0.5, bin_size = 3, jitter = FALSE, family = "binom"
) {
  if (!requireNamespace("PUcopula", quietly = TRUE)) {
    stop("PUcopula package is required", call. = FALSE)
  }
  dataTable <- data

  driver_strength <- lapply(as.list(driver_strength_factor), function(x) {
    max(1, floor(round(nrow(dataTable) * x)))
  })

  varnames <- colnames(dataTable)
  varnames_clean <- sub("(\\.cat\\.[0-9]+|\\.oriname)$", "", varnames)
  varnames_clean_unique <- unique(varnames_clean)

  num_jitter <- function(x, factor) base::jitter(as.numeric(x), factor = factor)
  if (isFALSE(jitter)) {
    # no-op
  } else if (is.numeric(jitter) && length(jitter) == 1) {
    dataTable <- lapply(dataTable, num_jitter, factor = jitter)
  } else if (is.list(jitter)) {
    for (iname in names(jitter)) {
      if (iname %in% varnames_clean) {
        idx <- which(varnames_clean == iname)
        for (id in idx) dataTable[[id]] <- num_jitter(dataTable[[id]], factor = jitter[[iname]])
      }
    }
  }

  if (is.numeric(bin_size) && length(bin_size) == 1) {
    dataTable <- lapply(dataTable, rank_bin_smooth, k = bin_size) |>
      as.data.frame() |> as.matrix()
    bin_size_list <- NULL
  } else if (is.numeric(bin_size) && length(bin_size) == length(varnames_clean_unique)) {
    bin_size_list <- as.list(bin_size); names(bin_size_list) <- varnames_clean_unique
  } else if (is.list(bin_size)) {
    bin_size_list <- bin_size
  } else {
    stop("bin_size must be a single numeric, a numeric vector of variable length, or a named list.")
  }

  if (!is.null(bin_size_list)) {
    fill_value <- if (length(bin_size_list) == 1 && is.null(names(bin_size_list))) bin_size_list[[1]] else NA
    new_bin_size <- bin_size_list[varnames_clean] |>
      sapply(function(x) if (is.null(x)) fill_value else x) |>
      as.numeric()
    dataTable <- mapply(rank_bin_smooth, dataTable, new_bin_size, SIMPLIFY = FALSE) |>
      as.data.frame() |> as.matrix()
  }

  driver_strength_list <- driver_strength
  fill_value <- if (length(driver_strength_list) == 1 && is.null(names(driver_strength_list))) driver_strength_list[[1]] else NA
  driver_strength_vec <- driver_strength_list[varnames_clean] |>
    sapply(function(x) if (is.null(x)) fill_value else x) |>
    as.numeric()

  PUcopula::PUCopula(
    family = family,
    pars.a = driver_strength_vec,
    patch  = "rook",
    data   = as.matrix(dataTable)
  )
}
