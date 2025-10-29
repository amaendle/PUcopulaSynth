#' @keywords internal
#' @noRd
rank_bin_smooth <- function(x, k) {
  n <- length(x)
  if (sum(!is.na(x)) < k) stop("Not enough non-NA values for chosen k")
  ranks <- rank(x, ties.method = "average", na.last = "keep")
  non_na_idx <- which(!is.na(ranks))
  sorted_non_na_idx <- order(ranks[non_na_idx])
  sorted_indices <- non_na_idx[sorted_non_na_idx]
  num_bins <- floor(length(sorted_indices) / k)
  breaks <- floor(seq(1, length(sorted_indices) + 1, length.out = num_bins + 1))
  smoothed_ranks <- rep(NA, n)
  for (i in seq_len(num_bins)) {
    start_idx <- breaks[i]; end_idx <- breaks[i + 1] - 1
    bin_indices <- sorted_indices[start_idx:end_idx]
    smoothed_ranks[bin_indices] <- mean(ranks[bin_indices])
  }
  smoothed_ranks
}

#' @keywords internal
#' @noRd
check_if_integer <- function(x) is.numeric(x) && all(x %% 1 == 0)

#' @keywords internal
#' @noRd
check_if_binary <- function(x, ignoreNA = TRUE) {
  if (ignoreNA) length(unique(stats::na.omit(x))) == 2 else length(unique(x)) == 2
}

#' @keywords internal
#' @noRd
check_if_trivial <- function(x) length(unique(x)) < 2

#' k-NN smoother for numeric vectors (no DataSHIELD thresholds)
#'
#' @param x Numeric vector (NAs allowed).
#' @param k Integer neighbours in [1, N-1] for non-missing values.
#' @return Numeric vector of same length with smoothed non-missing entries.
#' @export
#' @importFrom RANN nn2
#' @importFrom stats sd
knnsmoother <- function(x, k = 3) {
  if (!is.numeric(x)) stop("x must be numeric")
  out <- x
  idx <- which(!is.na(x))
  if (length(idx) < 2L) return(out)
  if (k < 1 || k > (length(idx) - 1)) {
    stop(sprintf("k must be in [1, %d] for non-missing entries.", length(idx) - 1), call. = FALSE)
  }
  xv <- x[idx]
  if (all(xv == xv[1])) return(out)

  xs <- (xv - mean(xv)) / stats::sd(xv)
  nearest <- RANN::nn2(xs, k = k)
  x.centroid <- vapply(seq_along(xs), function(i) mean(xs[nearest$nn.idx[i, 1:k]]), numeric(1))
  scaling <- stats::sd(xs) / stats::sd(x.centroid)
  out[idx] <- (x.centroid * scaling) * stats::sd(xv) + mean(xv)
  out
}

#' @keywords internal
#' @noRd
find_closest_cat <- function(..., contrast_tab) {
  vec <- c(...)
  dists <- apply(contrast_tab, 1, function(row) sum((vec - row)^2))
  rs <- names(which.min(dists))
  if (is.null(rs)) rs <- NA
  rs
}

#' Piecewise-linear transform used for discrete tables
#' @keywords internal
#' @noRd
#' @importFrom stats approx
transform_u <- function(u, x, p) {
  stopifnot(length(x) >= 2, length(p) == length(x))
  if (is.unsorted(x, strictly = TRUE)) stop("x must be strictly increasing")
  if (any(p < 0 | p > 1) || is.unsorted(p, strictly = TRUE)) stop("p must be strictly increasing in [0,1]")
  gaps <- diff(x)
  if (any(gaps <= 0)) stop("x must be strictly increasing")
  pk <- c(0, p[1:(length(x)-1)], 1)
  mids <- x[-length(x)] + gaps/2
  yk <- c(x[1] - 0.5 * gaps[1], mids, x[length(x)] + 0.5 * gaps[length(x)-1])
  stats::approx(pk, yk, xout = u, method = "linear", ties = "ordered", rule = 2)$y
}
