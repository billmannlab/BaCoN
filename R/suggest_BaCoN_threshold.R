#' @importFrom stats sd
#' @importFrom utils capture.output
#' @export suggest_BaCoN_threshold
#' @returns a function that suggests a threshold to compute BaCoN for high PCC values.

suggest_BaCoN_threshold <- function(.matrix, verbose = T) {

  .mean <- base::mean(.matrix, na.rm = T)
  .sd <- stats::sd(.matrix, na.rm = T)

  .d <- matrix(NA, nrow = 2, ncol = 4, dimnames = list(c("numeric threshold", "values to compute [%]"), c("none", 1:3)))

  .d["numeric threshold",c("none", 1:3)] <- c(0, .mean + 1:3 * .sd)
  .d["values to compute [%]",c("none", 1:3)] <- sapply(.d["numeric threshold",c("none", 1:3)],
                                                       \(.) sum(base::abs(.matrix) > as.numeric(.), na.rm = T)*100 / base::length(.matrix))
  .z <- colnames(.d)[last(base::which(.d["values to compute [%]",c("none", 1:3)] > 1))]

  if (verbose) {
    message("Possible thresholds:\n")
    message(base::paste0(utils::capture.output(round(.d, 2)), collapse = "\n"))}

  if (verbose) {message(base::paste0("\nSuggesting threshold ", .z, ".\n"))}

  return(.z)}
