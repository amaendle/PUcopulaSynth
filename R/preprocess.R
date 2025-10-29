#' Save original names and classes
#' @param data A data.frame
#' @return Character vector (names) / Named character vector (classes)
#' @export
save_original_varnames <- function(data) names(data)

#' @rdname save_original_varnames
#' @export
save_original_classes <- function(data) sapply(data, class)

#' Preprocess data before copula fitting
#'
#' Converts multi-level unordered factors to dummy variables and tags remaining
#' factor columns with '.oriname' while storing levels.
#'
#' @param data A data.frame
#' @return A list with `data` (processed) and `original_levels` (dummies + oriname)
#' @export
#' @importFrom dplyr select rename_with mutate across bind_cols
#' @importFrom tidyselect where matches
#' @importFrom caret dummyVars
preprocessData <- function(data) {
  dataTable <- data
  options(contrasts = c("contr.sum","contr.poly"))

  df.unordered <- dplyr::select(
    dataTable,
    tidyselect::where(~ is.factor(.) && !is.ordered(.) && nlevels(.) > 2)
  )
  unordered.names <- names(df.unordered)
  df.unordered <- dplyr::rename_with(df.unordered, ~ paste0(.x, rep(".cat.", length(.x))))
  df.unordered <- dplyr::mutate(df.unordered, dplyr::across(
    .cols = tidyselect::where(is.factor) & tidyselect::matches("\\.cat\\.$"),
    .fns = droplevels
  ))

  df.ordered <- dplyr::select(
    dataTable,
    tidyselect::where(~ is.ordered(.) && nlevels(.) > 2)
  )
  ordered.names <- names(df.ordered)
  df.ordered <- dplyr::rename_with(df.ordered, ~ paste0(.x, rep(".lev.", length(.x))))
  df.ordered <- dplyr::mutate(df.ordered, dplyr::across(
    .cols = tidyselect::where(is.ordered) & tidyselect::matches("\\.lev\\.$"),
    .fns = droplevels
  ))

  ####df.catlev <- dplyr::bind_cols(lapply(df.unordered, function(x) as.factor(as.integer(x))), df.ordered)
  df.catlev <- dplyr::bind_cols(df.unordered, df.ordered)

  if (ncol(df.catlev)==0)
    return(
      list(
        data = dataTable,
        original_levels = list(
          dummies = NULL,
          oriname = NULL
        )
      )
    )

  dmy <- caret::dummyVars(
    ~ .,
    data = dplyr::select(df.catlev, tidyselect::where(~ !(is.factor(.) && nlevels(.) <= 2))),
    sep = "",
    fullRank = TRUE
  )
  df.dmy <- data.frame(predict(dmy, newdata = df.catlev))
  #df.dmy[] <- lapply(df.dmy, function(x) as.factor(as.integer(x))) # why do we do this anyway
  ########################################################################################
  # Only convert columns ending with ".cat." followed by a number
  cols_cat <- grep("\\.cat\\.\\d+$", names(df.dmy), value = TRUE)
  df.dmy[cols_cat] <- lapply(df.dmy[cols_cat], function(x) as.factor(as.integer(x)))
  ########################################################################################

  dataTable <- dataTable[, !(names(dataTable) %in% c(unordered.names, ordered.names))]

  oriname_levels <- list()
  dataTable <- dplyr::mutate(
    dataTable,
    dplyr::across(
      .cols = tidyselect::where(is.factor),
      .fns = ~ {
        oriname_levels[[dplyr::cur_column()]] <<- levels(.)
        .
      }
    )
  )
  dataTable <- dplyr::rename_with(dataTable, ~ paste0(.x, ".oriname"))

  list(
    data = cbind(dataTable, df.dmy),
    original_levels = list(
      dummies = dmy$lvls,
      oriname = oriname_levels
    )
  )
}

#' Restore original factor structure after synthesis
#'
#' @param data A data.frame containing dummy/encoded columns
#' @param cat_dummy_levels The `dummies` element returned by [preprocessData()]
#' @return A data.frame with factors restored and columns ordered like input
#' @export
#' @importFrom stats contr.poly contr.sum
#' @importFrom dplyr select bind_cols
#' @importFrom tidyr all_of
postprocessData <- function(data, cat_dummy_levels) {
  dataTable <- data
  dmy_lvls <- cat_dummy_levels

  dummy_cols <- lapply(names(dmy_lvls), function(x) names(dataTable)[startsWith(names(dataTable), x)])
  dummy_cols <- stats::setNames(dummy_cols, names(dmy_lvls))
  n_dummies <- sapply(dummy_cols, function(x) as.integer(length(x) + 1))

  if (length(dummy_cols)==0) {
    return(data)
  }

  cat_levels <- dmy_lvls
  ordinal <- !endsWith(names(dummy_cols), ".cat.")

  contrast_cat_tab <- mapply(function(n_dum, is_ord) {
    if (is_ord) stats::contr.poly(n_dum) else stats::contr.sum(n_dum)
  }, n_dummies, as.list(ordinal), SIMPLIFY = FALSE)
  contrast_cat_tab <- mapply(function(tab, nam) { rownames(tab) <- nam; tab },
                             contrast_cat_tab, cat_levels, SIMPLIFY = FALSE)

  new_cols <- mapply(function(cols, contrasts) {
    apply(dataTable[, cols, drop = FALSE], 1, function(row) {
      find_closest_cat(as.numeric(row), contrast_tab = contrasts)
    })
  }, dummy_cols, contrast_cat_tab, SIMPLIFY = FALSE) |>
    as.data.frame() |>
    stats::setNames(names(dmy_lvls))

  new_cols <- mapply(factor, new_cols, cat_levels, SIMPLIFY = FALSE) |> as.data.frame()
  newdf <- dplyr::select(dataTable, -tidyr::all_of(unlist(dummy_cols)))
  newdf <- dplyr::bind_cols(newdf, new_cols)

  names(newdf) <- gsub("\\.(oriname|cat\\.|lev\\.)$", "", names(newdf))
  ordered_names <- unique(gsub("\\.(oriname|cat\\.(\\.L|\\.Q|\\.C|\\d+)|lev\\.(\\.L|\\.Q|\\.C|\\d+))$", "", colnames(dataTable)))
  ## ^^ \\. auch vor \\d+????????????????????
  newdf <- newdf[, match(ordered_names, names(newdf))]
  newdf
}
