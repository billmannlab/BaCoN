#' @import data.table
#' @export BaCoN_datatable
#' @returns A BaCoN-matrix of the input correlation matrix.

## ---- BaCoN_datatable function ----

BaCoN_datatable <- function(input_matrix, cf = 0.05, .th = NA, verbose = T) {
  .mean <- base::mean(input_matrix, na.rm = T)
  .sd <- stats::sd(input_matrix, na.rm = T)

  .d <- data.table::data.table(expand.grid(gene1 = rownames(input_matrix),
                                           gene2 = colnames(input_matrix),
                                           stringsAsFactors = T),
                               PCC = as.vector(input_matrix))

  baconize <- \(.vec, .cf) {
    .out <- base::rep(NA, base::length(.vec))
    i <- .vec >= 0 & !is.na(.vec)
    .out[i] <- base::sapply(.vec[i], \(x) base::sum(.vec > x - .cf, na.rm = T), simplify = T)
    if (base::sum(!i) > 0) {.out[!i] <- base::sapply(.vec[!i], \(x) base::sum(.vec < x + .cf, na.rm = T), simplify = T)}
    .out}

  if (is.na(.th)) {.i <- .d[, .I]}
  if (!is.na(.th)) {.i <- .d[, .I[abs(PCC) > (.mean + .th *.sd - cf)]]}

  if (verbose) {
    message("BaCoN_datatable started.")
    message("[", format(base::Sys.time(), "%X"), "] Entering phase 1/2...")
  }

  .grpn <- data.table::uniqueN(.d[, gene1])

  # Initialize start time
  .start_time <- Sys.time()

  # Perform the operation with progress and time estimation
  .d[.i, b1 := {
    .current <- Sys.time()
    .progress <- round(.GRP / .grpn * 100, 2) # Calculate the percentage of completion
    .elapsed <- difftime(.current, .start_time, units = "secs") # Calculate the elapsed time
    .estimated_total <- as.numeric(.elapsed) / (.GRP / .grpn) # Estimate total time based on the progress
    .eta <- .estimated_total - as.numeric(.elapsed) # Calculate the estimated remaining time
    cat("\r", .progress, "%. ETA: ", round(.eta, 0), " seconds", sep = "") # Print progress and estimated remaining time

    baconize(PCC, cf)}, by = gene1]

  cat("\r \r")

  #.d[.i, b1 := {cat("\r", round(.GRP/grpn*100, 2), "%"); baconize(PCC, cf)}, by = gene1]


  if (verbose) {
    message("[", format(base::Sys.time(), "%X"), "] Entering phase 2/2...")
  }

  .grpn <- data.table::uniqueN(.d[, gene2])
  .d[.i, b2 := {
    .current <- Sys.time()
    .progress <- round(.GRP / .grpn * 100, 2) # Calculate the percentage of completion
    .elapsed <- difftime(.current, .start_time, units = "secs") # Calculate the elapsed time
    .estimated_total <- as.numeric(.elapsed) / (.GRP / .grpn) # Estimate total time based on the progress
    .eta <- .estimated_total - as.numeric(.elapsed) # Calculate the estimated remaining time
    cat("\r", .progress, "%. ETA: ", round(.eta, 0), " seconds", sep = "") # Print progress and estimated remaining time


    baconize(PCC, cf)}, by = gene2]

  cat("\r \r")

  if (verbose) {
    message("[", format(base::Sys.time(), "%X"), "] Done.")
  }

  y <- base::sum(dim(input_matrix))
  .d[, BaCoN := data.table::fcase(PCC >= 0, 1 - (b1 + b2) / y,
                                  PCC < 0, -1 + (b1 + b2) / y)]
  if (!is.na(.th)) {.d[abs(PCC) < .mean + .th *.sd, BaCoN := NA]}


  out <- base::array(.d[, BaCoN], dim = dim(input_matrix), dimnames = dimnames(input_matrix))
  return(out)}
