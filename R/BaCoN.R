## ---- BaCoN function ----

#' @importFrom progress progress_bar
#' @importFrom future.apply future_apply
#' @importFrom future plan
#' @importFrom future multisession

#' @importFrom stats sd
#' @export BaCoN
#' @returns A BaCoN-matrix of the input correlation matrix.


BaCoN <- function(input_matrix, corr_f = 0.05,
                  threshold = "none",
                  #negative_th = "none",
                  #parallel = F,
                  n_cores = 1,
                  verbose = T,
                  showProgress = F,
                  simulate = F) {

  if (verbose & simulate) {
    message("Simulating BaCoN run.")
  }

  if (verbose & n_cores > 1) {
    message(paste0("BaCoN is set up to run with ", n_cores, " cores."))
    }

  if (verbose) {
    message(paste0("\nChosen threshold: ", threshold, ", chosen correction factor: ", corr_f, "."))
  }

  threshold <- ifelse(threshold == "none", NA, as.numeric(threshold))

  if (!is.na(threshold)) {
    threshold <- mean(input_matrix, na.rm = T) + threshold * stats::sd(input_matrix, na.rm = T)
    }

  if (simulate) {
    message("BaCoN simulation complete.")
    return(NULL)
  }

  if (!simulate) {

    appl_func <- ifelse(n_cores > 1,
                        future.apply::future_apply,
                        base::apply)

    if (n_cores > 1) {
      options(future.globals.maxSize = +Inf) # 1000*1024^2
      future::plan(future::multisession, workers = n_cores)
    }

    if (showProgress) {
      pb <- progress::progress_bar$new(width = 75,
                                       force = T,
                                       format = "[:bar] :percent (:timepoint), ETA: :eta")}

    if (verbose) {message("Ready to run.")}
    start_time <- Sys.time()
    if (verbose) {message("BaCoN started (", format(Sys.time(), "%X"), ").")}

    vectorized_BaCoN <- \(vector, cf = corr_f) {
      out_vec <- rep(NA, base::length(vector))
      i <- !is.na(vector)
      out_vec[i] <- base::sapply(vector[i], \(.) {sum(vector > . - cf, na.rm = T)})
      out_vec}

    .in <- list(pos = input_matrix, neg = input_matrix)
    .in$pos[input_matrix < max(threshold - corr_f, -corr_f, na.rm = T)] <- NA
    .in$neg[input_matrix >= min(-threshold + corr_f, corr_f, na.rm = T)] <- NA

    .out <- list(main = array(NA, dim = dim(input_matrix),
                              dimnames = dimnames(input_matrix)))
    .out$pos$h <- .out$main; .out$pos$v <- .out$main
    .out$neg$h <- .out$main; .out$neg$v <- .out$main

    if (showProgress) {
      pb$update(0.05,
                tokens = list(timepoint = format(Sys.time(), "%X")))
    }

    roi <- which(rowSums(!is.na(.in$pos)) > 1)
    .out$pos$h[roi,] <- t(appl_func(.in$pos[roi,], 1, vectorized_BaCoN))

    if (showProgress) {
      pb$update(0.25,
                tokens = list(timepoint = format(Sys.time(), "%X")))
    }

    coi <- which(base::colSums(!is.na(.in$pos)) > 1)
    .out$pos$v[,coi] <- appl_func(.in$pos[,coi], 2, vectorized_BaCoN)


    vectorized_BaCoN <- \(vector, cf = corr_f) {
      out_vec <- rep(NA, base::length(vector))
      i <- !is.na(vector)
      out_vec[i] <- sapply(vector[i], \(.) {sum(vector < . + cf, na.rm = T)})
      out_vec}

    .out$pos$merge <- (.out$pos$v + .out$pos$h)

    if (showProgress) {
      pb$update(0.5, tokens = list(timepoint = format(Sys.time(), "%X")))
      }

    roi <- which(rowSums(!is.na(.in$neg)) > 1)
    .out$neg$h[roi,] <- t(appl_func(.in$neg[roi,], 1, vectorized_BaCoN))

    if (showProgress) {
      pb$update(0.75, tokens = list(timepoint = format(Sys.time(), "%X")))
      }

    coi <- which(colSums(!is.na(.in$neg)) > 1)
    .out$neg$v[,coi] <- appl_func(.in$neg[,coi], 2, vectorized_BaCoN)

    .out$neg$merge <- (.out$neg$v + .out$neg$h)

    y <- sum(dim(input_matrix))

    i_pos <- which(!is.na(.in$pos) & input_matrix >= threshold)
    .out$main[i_pos] <- 1 - (.out$pos$merge[i_pos] / y)

    i_neg <- which(!is.na(.in$neg) & input_matrix < -threshold)

    if (!base::length(base::intersect(i_pos, i_neg)) == 0) {
      print(base::intersect(i_pos, i_neg))
    }

    .out$main[i_neg] <- - 1 + (.out$neg$merge[i_neg] / y)

    if (showProgress) {
      pb$finished
    }

    if (verbose) {
      message(paste0("\nCompleted after ",
                     round(difftime(base::Sys.time(),
                                    start_time,
                                    units = "min"), 2), " minutes."))}
    return(.out$main)
  }
}
