#' @import data.table
#' @importFrom progress progress_bar
#'
#' @export BaCoN
#' @returns A BaCoN-matrix of the input correlation matrix.

## ---- BaCoN function ----

BaCoN <- \(input_matrix,
           cf = 0.05,
           verbose = T,
           show_progress = T,
           n_threads = 1,
           detailed_output = F) {

  .pbformat <- "[:bar] :percent (:current/:total, :tick_rate), elapsed: :elapsedfull, ETA: :eta"

  setDTthreads(n_threads)

  .genespace <- length(input_matrix)
  .nrow <- nrow(input_matrix)
  .ncol <- ncol(input_matrix)
  .rownames <- rownames(input_matrix)
  .colnames <- colnames(input_matrix)
  .y <- .nrow + .ncol
  .rowwise_progress <- rep(NA, .nrow)
  .colwise_progress <- rep(NA, .ncol)

  .data <- list(ID = seq_along(1:.genespace),
                gene1 = rep(.rownames, .ncol),
                gene2 = rep(.colnames, each = .nrow),
                PCC = as.vector(input_matrix))


  setDT(.data)

  setkey(.data, gene1, physical = T)

  if (verbose) {message("Entering phase 1/2 (rowwise computation)...")}
  if (show_progress) {.pb <- progress::progress_bar$new(format = .pbformat,
                                                        total = .nrow,
                                                        width = 75, force = T)}

  .start <- base::Sys.time()
  .data[, bacon_rowwise := {
    if (detailed_output) {
      .grp <- .GRP
      .rowwise_progress[[.grp]] <<- difftime(Sys.time(), .start, units = "secs")
    }
    if (show_progress) {.pb$tick()}
    baconize(PCC, cf)}, by = gene1]

  setkey(.data, ID, physical = T)
  setkey(.data, gene2, physical = T)

  if (verbose) {message("Entering phase 2/2 (columnwise computation)...")}
  if (show_progress) {.pb <- progress::progress_bar$new(format = .pbformat,
                                                        total = .ncol,
                                                        width = 75, force = T)}

  .start <- base::Sys.time()
  .data[, bacon_colwise := {
    if (detailed_output) {
      .grp <- .GRP
      .colwise_progress[[.grp]] <<- difftime(Sys.time(), .start, units = "secs")
    }
    if (show_progress) {.pb$tick()}
    baconize(PCC, cf)
  }, by = gene2]

  setkey(.data, ID, physical = T)

  ###
  if(!all(.data[, ID] == seq_along(1:length(input_matrix)))) {
    warning("Check order of the elements!")
  }
  ###

  .data[, BaCoN := data.table::fcase(PCC >= 0, 1 - (bacon_rowwise + bacon_colwise) / .y,
                                     PCC < 0, -1 + (bacon_rowwise + bacon_colwise) / .y)]

  .matrix <- base::array(data = .data[, get("BaCoN")],
                         dim = c(.nrow, .ncol),
                         dimnames = list(.rownames, .colnames))

  if (detailed_output) {
    return(list(correlation_matrix = input_matrix,
                correction_factor = cf,
                genespace = .genespace,
                nrow = .nrow,
                ncol = .ncol,
                rownames = .rownames,
                colnames = .colnames,
                rowwise_progress = .rowwise_progress,
                colwise_progress = .colwise_progress,
                data = .data,
                bacon_matrix = .matrix))
  } else {
    return(.matrix)

  }}
