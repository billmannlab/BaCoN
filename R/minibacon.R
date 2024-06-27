## ---- MiniBaCoN function ----

#' @import data.table
#'

#' @returns A BaCoN-matrix of the input correlation matrix.

miniBaCoN <- function(.matrix, cf = 0.05, .th = NA, verbose = T) {
  .mean <- base::mean(.matrix, na.rm = T)
  .sd <- stats::sd(.matrix, na.rm = T)

  .d <- data.table::data.table(expand.grid(gene1 = rownames(.matrix),
                               gene2 = colnames(.matrix),
                               stringsAsFactors = T),
                   PCC = as.vector(.matrix))

  baconize <- \(.vec, .cf) {
    .out <- base::rep(NA, base::length(.vec))
    i <- .vec >= 0 & !is.na(.vec)
    .out[i] <- base::sapply(.vec[i], \(x) sum(.vec > x - .cf, na.rm = T), simplify = T)
    if (sum(!i) > 0) {.out[!i] <- base::sapply(.vec[!i], \(x) base::sum(.vec < x + .cf, na.rm = T), simplify = T)}
    .out}

  if (is.na(.th)) {.i <- .d[, .I]}
  if (!is.na(.th)) {.i <- .d[, .I[abs(PCC) > (.mean + .th *.sd - cf)]]}

  if (verbose) {
    message("miniBaCoN started... (", format(base::Sys.time(), "%X"), ").")
  }


  .d[.i, b1 := baconize(PCC, cf), by = gene1]

  if (verbose) {
    message("Halfway done... (", format(base::Sys.time(), "%X"), ").")
  }

  .d[.i, b2 := baconize(PCC, cf), by = gene2]

  if (verbose) {
    message("Done. (", format(base::Sys.time(), "%X"), ").")
    }

  y <- sum(dim(.matrix))
  .d[, BaCoN := data.table::fcase(PCC >= 0, 1 - (b1 + b2) / y,
                      PCC < 0, -1 + (b1 + b2) / y)]
  if (!is.na(.th)) {.d[abs(PCC) < .mean + .th *.sd, BaCoN := NA]}


  out <- base::array(.d[, BaCoN], dim = dim(.matrix), dimnames = dimnames(.matrix))
  return(out)}
