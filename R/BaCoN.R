#' @export BaCoN
#' @returns A BaCoN-matrix of the input correlation matrix.

## ---- BaCoN function ----

BaCoN <- function(
    input_matrix,
    cf = 0.05,
    verbose = T) {

  .y <- nrow(input_matrix) + ncol(input_matrix)

  output <- list(merged = array(data = NA,
                                dim = dim(input_matrix),
                                dimnames = dimnames(input_matrix)))

  output$rowwise <- rownames(input_matrix) |>
    purrr::map(~ input_matrix[.x,]) |>
    purrr::map(baconize, .cf = cf, .progress = ifelse(verbose, "1/2...", F))
  output$rowwise <- do.call(rbind, output$rowwise)

  output$colwise <- colnames(input_matrix) |>
    purrr::map(~ input_matrix[,.x]) |>
    purrr::map(baconize, .cf = cf, .progress = ifelse(verbose, "2/2...", F))

    output$colwise <- do.call(cbind, output$colwise)

  .i <- !is.na(input_matrix) & input_matrix >= 0
  output$merged[.i] <- 1 - (output$rowwise[.i] + output$colwise[.i]) / .y


  .i <- !is.na(input_matrix) & input_matrix < 0
  output$merged[.i] <- -1 + (output$rowwise[.i] + output$colwise[.i]) / .y

  output <- output$merged

  attr(output, "BaCoN_attributes") <- list(y = .y)

  return(output)
}
