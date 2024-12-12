#' @import data.table
#' @import job
#' @importFrom progress progress_bar
#' @importFrom stats cor

#' @export BaCoN_pipeline
#' @returns xxx



BaCoN_pipeline <- \(expression_matrix,
                    effect_matrix,
                    cor_method = "pearson",
                    cache_path = NA,
                    bacon_correction_factor = 0.05,
                    show_progress = T,
                    verbose = T) {

  if (!is.na(cache_path)) {
    if (!base::dir.exists(cache_path)) {base::dir.create(cache_path, recursive = T, showWarnings = F)}
  }


  .nrow <- ncol(expression_matrix)
  .ncol <- ncol(effect_matrix)
  .genespace <- prod(.nrow, .ncol)
  .rownames <- colnames(expression_matrix)
  .colnames <- colnames(effect_matrix)
  .cl_intersect <- intersect(rownames(effect_matrix), rownames(expression_matrix))


  if (verbose) {
    message(base::paste0("[", base::format(base::Sys.time(), "%X"), "]"))
    if (is.na(cache_path)) {message("Not saving cache.")}
    if (!is.na(cache_path)) {message(base::paste0("Cache directory: ", cache_path))}
    message(base::paste0("Matrix dimensions: ", .nrow, " x ", .ncol, " -> ", .genespace, " gene pairs."))
    message(base::paste("A", cor_method, "correlation matrix will be computed based on", length(.cl_intersect), "cell lines.", sep = " "))
  }


  #message("\n")
  # message(paste0("PCC cache exists: ", file.exists(file.path(cache_path, str_c(pcc_fname, ".rds")))))
  #message("\n")

  if (is.na(cache_path)) { # no cache path provided. Do not save the files.
    .cormat <- cor(expression_matrix[.cl_intersect,],
                   effect_matrix[.cl_intersect,],
                   use = "pairwise.complete.obs",
                   method = cor_method)

    .bacon <- BaCoN(input_matrix = .cormat,
                    cf = bacon_correction_factor,
                    verbose = verbose,
                    show_progress = show_progress,
                    n_threads = 1, detailed_output = T)}

  if (!is.na(cache_path)) { # cache path provided. Import matrices if available, otherwise compute.

    .cormat_file <- file.path(cache_path, "correlation_matrix.rds")
    .bacon_file <- file.path(cache_path, "bacon_matrix.rds")


    if (file.exists(.cormat_file)) {
      .cormat <- readRDS(.cormat_file)
    } else {
      .cormat <- cor(expression_matrix[.cl_intersect,],
                     effect_matrix[.cl_intersect,],
                     use = "pairwise.complete.obs",
                     method = cor_method)
      saveRDS(.cormat, .cormat_file)
    }

    if (file.exists(.bacon_file)) {
      .bacon <- list(bacon_matrix = readRDS(.bacon_file))
    } else {
      .bacon <- BaCoN(input_matrix = .cormat,
                      cf = bacon_correction_factor,
                      verbose = verbose,
                      show_progress = show_progress,
                      n_threads = 1, detailed_output = T)

      for (. in c("rowwise_progress", "colwise_progress", "bacon_matrix")) {
        saveRDS(.bacon[[.]], file.path(cache_path, paste0(., ".rds")))}
    }
  }

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
