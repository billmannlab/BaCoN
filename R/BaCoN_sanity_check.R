BaCoN_sanity_check <- \(cormat, bacon_mat, cf = 0.05, n = 1000) {

  y <- sum(dim(cormat)[1:2])
  .d <- data.table(expression_gene = base::sample(rownames(cormat), n),
                   effect_gene = base::sample(colnames(cormat), n))

  .d[, `:=`(cor = base::mapply(\(.g1, .g2) {cormat[.g1,.g2]}, expression_gene, effect_gene),
            bacon = base::mapply(\(.g1, .g2) {bacon_mat[.g1,.g2]}, expression_gene, effect_gene))]

  .d[cor >= 0, bacon_check := base::mapply(\(.g1, .g2) {
    1 - (sum(cormat[.g1,] > cormat[.g1,.g2]-cf, na.rm = T) + sum(cormat[,.g2] > cormat[.g1,.g2]-cf, na.rm = T)) / y}, expression_gene, effect_gene)]
  .d[cor < 0, bacon_check := base::mapply(\(.g1, .g2) {
    -1 + (sum(cormat[.g1,] < cormat[.g1,.g2]+cf, na.rm = T) + sum(cormat[,.g2] < cormat[.g1,.g2]+cf, na.rm = T)) / y}, expression_gene, effect_gene)]
  .d[, base::all(bacon == bacon_check)]
}
