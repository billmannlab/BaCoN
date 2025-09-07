#' @export BaCoN_v2
#' @returns A BaCoN-matrix of the input correlation matrix.

## ---- BaCoN_v2 function ----

BaCoN_v2 <- function(
    input_matrix,
    cf = 0.05,
    verbose = T) {

  .attr <- list(genespace = length(input_matrix),
                nrow = nrow(input_matrix),
                ncol = ncol(input_matrix),
                rownames = rownames(input_matrix),
                colnames = colnames(input_matrix))

  .attr$y <- .attr$nrow + .attr$ncol

  output <- purrr::imap(purrr::set_names(c("rowwise", "colwise", "merged")),
                 ~ array(dim = dim(input_matrix), dimnames = dimnames(input_matrix)))

  purrr::walk(rownames(input_matrix), ~ {output$rowwise[.x,] <<- baconize(input_matrix[.x,], .cf = cf)}, .progress = ifelse(verbose, "Rowwise (1/2)...", F))
  purrr::walk(colnames(input_matrix), ~ {output$colwise[,.x] <<- baconize(input_matrix[,.x], .cf = cf)}, .progress = ifelse(verbose, "Column-wise (2/2)...", F))

  .i <- !is.na(input_matrix) & input_matrix >= 0
  output$merged[.i] <- 1 - (output$rowwise[.i] + output$colwise[.i]) / .attr$y


  .i <- !is.na(input_matrix) & input_matrix < 0
  output$merged[.i] <- -1 + (output$rowwise[.i] + output$colwise[.i]) / .attr$y

  output <- output$merged

  attr(output, "BaCoN_attributes") <- .attr

  return(output)
  }
