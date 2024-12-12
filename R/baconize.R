baconize <- \(.vec, .cf) {
  .out <- base::rep(NA, base::length(.vec))
  i <- .vec >= 0 & !is.na(.vec)
  .out[i] <- base::sapply(.vec[i], \(x) base::sum(.vec > x - .cf, na.rm = T), simplify = T)
  if (base::sum(!i) > 0) {.out[!i] <- base::sapply(.vec[!i], \(x) base::sum(.vec < x + .cf, na.rm = T), simplify = T)}
  .out}
