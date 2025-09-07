
#' @export BaCoN_pipeline
#' @returns A convenience wrapper that first computes a correlation matrix
#' followed by a BaCoN matrix and stores these objects in a cache directory.
#' The BaCoN computation is executed in an RStudio job intance, which allows
#' parallel computation of multiple BaCoN matrices on.

## ---- BaCoN_pipeline function ----

BaCoN_pipeline <- function(expression_matrix,
                           effect_matrix,
                           cor_method = "pearson",
                           cache_path,
                           bacon_correction_factor = 0.05,
                           pairs_to_remove = NULL,
                           cache_file_suffix = "no_proximity_pairs_1e7_bp",
                           verbose = T) {

  stopifnot("Cache file path required" = !missing(cache_path))

  if (!missing(cache_path)) {dir.create(cache_path, recursive = T, showWarnings = F)}

  .attr <- list(nrow = ncol(expression_matrix),
                ncol = ncol(effect_matrix),
                genespace = prod(ncol(expression_matrix), ncol(effect_matrix)),
                rownames = colnames(expression_matrix),
                colnames = colnames(effect_matrix),
                correlation_method = cor_method,
                intersecting_cell_lines = intersect(rownames(effect_matrix),
                                                    rownames(expression_matrix)),
                bacon_correction_factor = bacon_correction_factor)


  stopifnot("Make sure the cell lines of the expression and effect matrix match!" = (
    all(rownames(expression_matrix) == .attr$intersecting_cell_lines) &
      all(rownames(effect_matrix) == .attr$intersecting_cell_lines))
  )

  if (verbose) {
    message(stringr::str_glue("Correlation method: {cor_method}. {length(.attr$intersecting_cell_lines)} cell lines."))
    message(stringr::str_glue("Matrix dimensions: {.attr$nrow} x {.attr$ncol} -> {.attr$genespace} gene pairs."))
  }

  .files <- c(cormat = "correlation_matrix", bacon = "bacon_matrix")


  if (!is.null(pairs_to_remove)) {
    .files <- c(.files,
                cormat_no_neighbors = stringr::str_glue("correlation_matrix_{cache_file_suffix}"),
                bacon_no_neighbors = stringr::str_glue("bacon_matrix_{cache_file_suffix}"))
  }

  .files <- purrr::map_chr(.files, ~ file.path(cache_path, stringr::str_glue("{.x}.rds")))

  if (!all(purrr::map_lgl(.files, file.exists))) {
    job::job({

      # 1

      if (file.exists(.files[["cormat"]])) {
        if (verbose) {message("Import correlation matrix from cache.")}
        .cormat <- readRDS(.files[["cormat"]])
      } else {
        if (verbose) {message("Computing correlation matrix...")}
        .cormat <- cor(expression_matrix,
                       effect_matrix,
                       use = "pairwise.complete.obs",
                       method = cor_method)
        saveRDS(.cormat, .files[["cormat"]])
      }

      # 2
      if (!is.null(pairs_to_remove)) {
        if (file.exists(.files[["cormat_no_neighbors"]])) {
          if (verbose) {message("Import proximity-adjusted correlation matrix from cache.")}
          .cormat_no_neighbors <- readRDS(.files[["cormat_no_neighbors"]])
        } else {
          if (verbose) {message("Computing proximity-adjusted correlation matrix...")}
          .cormat_no_neighbors <- .cormat

          .i <- sort_gene_pairs(rep(rownames(.cormat_no_neighbors), ncol(.cormat_no_neighbors)),
                                rep(colnames(.cormat_no_neighbors), each = nrow(.cormat_no_neighbors)))

          .cormat_no_neighbors[which(.i %in% pairs_to_remove)] <- NA

          saveRDS(.cormat_no_neighbors, .files[["cormat_no_neighbors"]])
        }
      }

      # 3

      if (file.exists(.files[["bacon"]])) {
        if (verbose) {message("Import BaCoN matrix from cache.")}
        .bacon <- list(bacon_matrix = readRDS(.files[["bacon"]]))
      } else {
        if (verbose) {message("Computing BaCoN matrix...")}
        .bacon <- BaCoN_v2(input_matrix = .cormat,
                           cf = bacon_correction_factor,
                           verbose = verbose)

        saveRDS(.bacon, .files[["bacon"]])
      }

      # 4

      if (!is.null(pairs_to_remove)) {
        if (file.exists(.files[["bacon_no_neighbors"]])) {
          if (verbose) {message("Import proximity-adjusted BaCoN matrix from cache.")}
          .bacon_no_neighbors <- list(bacon_matrix = readRDS(.files$bacon_no_neighbors))
        } else {
          if (verbose) {message("Computing proximity-adjusted BaCoN matrix...")}

          .bacon_no_neighbors <- BaCoN_v2(input_matrix = .cormat_no_neighbors,
                                          cf = bacon_correction_factor,
                                          verbose = verbose)

          saveRDS(.bacon_no_neighbors, .files[["bacon_no_neighbors"]])
        }
      }

      job::export("none")

    }, import = c("BaCoN_v2", "baconize", "sort_gene_pairs",
                  "effect_matrix", "expression_matrix",
                  ".files",
                  "pairs_to_remove",
                  "bacon_correction_factor", "cor_method", "verbose", "cache_path"),
    packages = c("data.table"),
    title = stringr::str_glue("BaCoN ({basename(cache_path)})")
    )

    output <- NULL # if the bacon matrix has not been computed yet, return NULL instead

  } else {
    output <- .files
    if (is.null(pairs_to_remove)) {
      output <- purrr::set_names(output, c("cor", "BaCoN"))
    } else {
      output <- purrr::set_names(output, c("cor", "BaCoN", "cor_no_neighbors", "BaCoN_no_neighbors"))
    }

    output <- purrr::map(output, readRDS) |> abind::abind(along = 3)
    attr(output, "BaCoN_attributes") <- .attr

  }

  return(output)

}
