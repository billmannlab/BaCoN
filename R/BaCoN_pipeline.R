#' @import data.table
#' @import job
#' @importFrom progress progress_bar
#' @importFrom stats cor

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
                    show_progress = T,
                    verbose = T) {

  if (!is.na(cache_path)) {
    if (!base::dir.exists(cache_path)) {
      base::dir.create(cache_path, recursive = T, showWarnings = F)}
  }


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
    #message(base::paste0("[", base::format(base::Sys.time(), "%X"), "]"))
    #if (is.na(cache_path)) {message("Not saving cache.")}
    message(
      ifelse(
        is.na(cache_path),
        str_c("Running pipeline without cache.",
              "It is recommended to provide a cache filepath for large matrices.",
              sep = "\n"),
        base::paste0("Cache directory: ", cache_path)))

    message(base::paste0("Correlation method: ",
                         cor_method, ". ", length(.cl_intersect), " cell lines."))
    message(base::paste0("Matrix dimensions: ",
                         .nrow, " x ", .ncol, " -> ", .genespace, " gene pairs."))
  }

  if (is.na(cache_path)) { # no cache path provided. Do not save the files.
    if (verbose) {message("Computing correlation matrix...")}

    .cormat <- cor(expression_matrix,effect_matrix,
                   use = "pairwise.complete.obs",
                   method = cor_method)

    if (verbose) {message("Computing BaCoN matrix...")}

    .bacon <- BaCoN(input_matrix = .cormat,
                    cf = bacon_correction_factor,
                    verbose = verbose,
                    show_progress = show_progress,
                    n_threads = 1, detailed_output = T)}

  # cache path provided. Import matrices if available, otherwise compute.
  if (!is.na(cache_path)) {

    .cormat_file <- file.path(cache_path, "correlation_matrix.rds")
    .bacon_file <- file.path(cache_path, "bacon_matrix.rds")


    if (file.exists(.cormat_file)) {
      if (verbose) {message("Import correlation matrix from cache.")}
      .cormat <- readRDS(.cormat_file)
    } else {
      .cormat <- cor(expression_matrix[.cl_intersect,],
                     effect_matrix[.cl_intersect,],
                     use = "pairwise.complete.obs",
                     method = cor_method)
      saveRDS(.cormat, .cormat_file)
    }

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
      title = base::paste0("BaCoN (", cache_path, ")"))

      # if the bacon matrix has not been computed yet, return NULL instead.
      .bacon <- list(bacon_matrix = NULL)


      if (verbose) {message("Returning NULL instead of the BaCoN matrix.")}
    }}

  return(list(
    nrow = .nrow,
    ncol = .ncol,
    genespace = .genespace,
    rownames = .rownames,
    colnames = .colnames,
    correlation_method = cor_method,
    intersecting_cell_lines = .cl_intersect,
    bacon_correction_factor = bacon_correction_factor,
    correlation_matrix = .cormat,
    bacon_matrix = .bacon$bacon_matrix))

}
