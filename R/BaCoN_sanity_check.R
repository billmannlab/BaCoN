#' @export BaCoN_sanity_check

BaCoN_sanity_check <- function(cormat, bacon_mat, cf = 0.05, n = 1000) {

  y <- sum(dim(cormat)[1:2])

  .x <- data.frame(expression_gene = sample(rownames(cormat), n),
                   effect_gene = sample(colnames(cormat), n))

  .x$cor <- purrr::map2_dbl(.x$expression_gene, .x$effect_gene, \(.g1, .g2) {cormat[.g1,.g2]})
  .x$bacon <- purrr::map2_dbl(.x$expression_gene, .x$effect_gene, \(.g1, .g2) {bacon_mat[.g1,.g2]})


  .i <- which(.x$cor >= 0)

  .x[.i,"bacon_check"] <- purrr::map2_dbl(.x$expression_gene[.i], .x$effect_gene[.i],
                                              \(.g1, .g2) {1 - (sum(cormat[.g1,] > cormat[.g1,.g2] - cf, na.rm = T) + sum(cormat[,.g2] > cormat[.g1,.g2] - cf, na.rm = T)) / y})


  .i <- which(.x$cor < 0)

  .x[.i,"bacon_check"] <- purrr::map2_dbl(.x$expression_gene[.i], .x$effect_gene[.i],
                                             \(.g1, .g2) {-1 + (sum(cormat[.g1,] < cormat[.g1,.g2] + cf, na.rm = T) + sum(cormat[,.g2] < cormat[.g1,.g2] + cf, na.rm = T)) / y})


#  .d[cor >= 0, bacon_check := purrr::map2_dbl(expression_gene, effect_gene,
#                                              \(.g1, .g2) {1 - (sum(cormat[.g1,] > cormat[.g1,.g2] - cf, na.rm = T) + sum(cormat[,.g2] > cormat[.g1,.g2] - cf, na.rm = T)) / y})]
#  .d[cor < 0, bacon_check := purrr::map2_dbl(expression_gene, effect_gene,
#                                             \(.g1, .g2) {-1 + (sum(cormat[.g1,] < cormat[.g1,.g2] + cf, na.rm = T) + sum(cormat[,.g2] < cormat[.g1,.g2] + cf, na.rm = T)) / y})]

  return(all(.x$bacon == .x$bacon_check, na.rm = T))
}
