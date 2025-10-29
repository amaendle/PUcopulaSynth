#' Estimate marginal models
#'
#' Fits logspline marginals for numeric/ordered variables, and empirical
#' probability tables for binary or trivial variables. Optional k-NN smoothing
#' is applied to numeric columns.
#'
#' @param data Preprocessed data.frame
#' @param method Character; for numeric and ordered factors (currently "spline")
#' @param k Numeric scalar, vector, or named list for k-NN smoothing
#' @return A named list of marginal models (each element has `qfun`)
#' @export
#' @importFrom logspline logspline
estimateMarginals <- function(data, method = "spline", k = NULL) {
  if (!requireNamespace("logspline", quietly = TRUE)) stop("logspline package is required")

  method <- rep_len(method, 2)
  dataTable <- data

  varnames <- colnames(dataTable)
  varnames_clean <- sub("(\\.cat\\.[0-9]+|\\.oriname)$", "", varnames)  # what about .lev.???
  varnames_clean_unique <- unique(varnames_clean)

  if (is.numeric(k) && length(k) == 1) {
    dataTable <- lapply(
      dataTable,
      function(x) if (is.numeric(x)) knnsmoother(x, k = k) else x
    ) |> as.data.frame()
    k_list <- NULL
  } else if (is.numeric(k) && length(k) == length(varnames_clean_unique)) {
    k_list <- as.list(k); names(k_list) <- varnames_clean_unique
  } else if (is.list(k)) {
    k_list <- k
  } else if (is.numeric(k) && length(k) != length(varnames_clean_unique) && length(k) != 1) {
    stop("k must be length 1, length = #variables, or a named list.")
  } else if (is.null(k)) {
    k_list <- NULL
  }

  if (!is.null(k_list)) {
    fill_value <- if (length(k_list) == 1 && is.null(names(k_list))) k_list[[1]] else NA
    k_vec <- k_list[varnames_clean] |>
      sapply(function(x) if (is.null(x)) fill_value else x) |>
      as.numeric()
    dataTable <- mapply(
      function(x, k_i) if (is.numeric(x)) knnsmoother(x, k = k_i) else x,
      dataTable, k_vec, SIMPLIFY = FALSE
    ) |> as.data.frame()
  }

  marginals <- mapply(function(col, varname) {
    if (method[1] == "spline") {
      mod <-
        if (check_if_trivial(col)) {
          prop.table(table(col))
        } else if (check_if_binary(col)) {
          prop.table(table(col))
        } else if (is.factor(col)) {
          if (method[2] == "spline") {
            tryCatch(
              logspline::logspline(jitter(as.numeric(as.character(col_jittered))), amount=0.5),
              error = function(e) stop(sprintf("logspline failed for '%s': %s", varname, conditionMessage(e)), call. = FALSE)
            )
          } else {
            prop.table(table(col))
          }
        } else {
          if (is.integer(col)) {
            # jitter to make estimation numerically safer
            col <- jitter(col, amount=0.5)
          }
          tryCatch(
            logspline::logspline(col),
            error = function(e) stop(sprintf("logspline failed for '%s': %s", varname, conditionMessage(e)), call. = FALSE)
          )
        }
      list(qfun = mod)
    } else {
      stop("Only 'spline' supported for numeric at the moment.")
    }
  }, dataTable, names(dataTable))

  names(marginals) <- names(dataTable)
  marginals
}
