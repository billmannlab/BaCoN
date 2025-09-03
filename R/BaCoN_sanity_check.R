#' @export BaCoN_sanity_check

BaCoN_sanity_check <- function(cormat, bacon_mat, cf = 0.05, n = 1000) {

  y <- sum(dim(cormat)[1:2])

  .d <- data.table::data.table(expression_gene = sample(rownames(cormat), n),
                   effect_gene = sample(colnames(cormat), n))

  .d[, `:=`(cor = purrr::map2_dbl(expression_gene, effect_gene, \(.g1, .g2) {cormat[.g1,.g2]}),
            bacon = purrr::map2_dbl(expression_gene, effect_gene, \(.g1, .g2) {bacon_mat[.g1,.g2]}))]

  .d[cor >= 0, bacon_check := purrr::map2_dbl(expression_gene, effect_gene,
                                              \(.g1, .g2) {1 - (sum(cormat[.g1,] > cormat[.g1,.g2] - cf, na.rm = T) + sum(cormat[,.g2] > cormat[.g1,.g2] - cf, na.rm = T)) / y})]
  .d[cor < 0, bacon_check := purrr::map2_dbl(expression_gene, effect_gene,
                                             \(.g1, .g2) {-1 + (sum(cormat[.g1,] < cormat[.g1,.g2] + cf, na.rm = T) + sum(cormat[,.g2] < cormat[.g1,.g2] + cf, na.rm = T)) / y})]

  return(.d[, all(bacon == bacon_check)])
}
