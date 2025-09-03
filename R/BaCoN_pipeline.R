
#' @export BaCoN_pipeline
#' @returns A convenience wrapper that first computes a correlation matrix
#' followed by a BaCoN matrix and stores these objects in a cache directory.
#' The BaCoN computation is executed in an RStudio job intance, which allows
#' parallel computation of multiple BaCoN matrices on.

## ---- BaCoN_pipeline function ----

BaCoN_pipeline <- \(expression_matrix,
                    effect_matrix,
                    cor_method = "pearson",
                    cache_path = NA,
                    bacon_correction_factor = 0.05,
                    steps = c("correlation", "BaCoN"),
                    show_progress = T,
                    verbose = T) {




  .nrow <- ncol(expression_matrix)
  .ncol <- ncol(effect_matrix)
  .genespace <- prod(.nrow, .ncol)
  .rownames <- colnames(expression_matrix)
  .colnames <- colnames(effect_matrix)
  .cl_intersect <- intersect(rownames(effect_matrix),
                             rownames(expression_matrix))

  if (!(all(rownames(expression_matrix) == .cl_intersect) &
        all(rownames(effect_matrix) == .cl_intersect))) {
    warning("Make sure the cell lines of the expression and effect matrix match!")
  }

  if (verbose) {
    message(paste0("Correlation method: ",
                         cor_method, ". ", length(.cl_intersect), " cell lines."))
    message(paste0("Matrix dimensions: ",
                         .nrow, " x ", .ncol, " -> ", .genespace, " gene pairs."))
  }

  if (!is.na(cache_path)) {
    if (!dir.exists(cache_path)) {
      dir.create(cache_path, recursive = T, showWarnings = F)}
  }

  .cormat_file <- file.path(cache_path, "correlation_matrix.rds")
  .bacon_file <- file.path(cache_path, "bacon_matrix.rds")

  if ("correlation" %in% steps) {
    if (file.exists(.cormat_file)) {
      if (verbose) {message("Import correlation matrix from cache.")}
      .cormat <- readRDS(.cormat_file)
    } else {
      if (verbose) {message("Computing correlation matrix...")}
      .cormat <- cor(expression_matrix,
                     effect_matrix,
                     use = "pairwise.complete.obs",
                     method = cor_method)
      saveRDS(.cormat, .cormat_file)
    }
  }

  if ("BaCoN" %in% steps) {
    if (file.exists(.bacon_file)) {
      if (verbose) {message("Import BaCoN matrix from cache.")}
      .bacon <- list(bacon_matrix = readRDS(.bacon_file))
    } else {
      if (verbose) {message("Starting BaCoN job...")}

      job::job({
        .bacon <- BaCoN(input_matrix = .cormat,
                        cf = bacon_correction_factor,
                        verbose = verbose,
                        show_progress = show_progress,
                        n_threads = 1, detailed_output = T)

        for (. in c("rowwise_progress", "colwise_progress", "bacon_matrix")) {
          saveRDS(.bacon[[.]], file.path(cache_path, paste0(., ".rds")))}

        job::export("none")

      }, import = c("BaCoN", "baconize", ".cormat",
                    "bacon_correction_factor", "cache_path",
                    "verbose", "show_progress"),
      packages = c("data.table", "progress"),
      title = paste0("BaCoN (", basename(cache_path), ")"))

      # if the bacon matrix has not been computed yet, return NULL instead.
      .bacon <- list(bacon_matrix = NULL)
    }
  }

  output <- list(nrow = .nrow,
                 ncol = .ncol,
                 genespace = .genespace,
                 rownames = .rownames,
                 colnames = .colnames,
                 correlation_method = cor_method,
                 intersecting_cell_lines = .cl_intersect,
                 bacon_correction_factor = bacon_correction_factor)

  if ("correlation" %in% steps) {
    output$correlation_matrix <- .cormat
  }
  if ("BaCoN" %in% steps) {
    output$bacon_matrix <- .bacon$bacon_matrix
  }

  return(output)

}
