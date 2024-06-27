
diversity_index <- \(set1, set2, .i = 0.5) {
  if (l(set1) != l(set2)) {message("Set size is not equal!")}
  x <- cumsum(sort(table(c(set1, set2)), decreasing = T))
  max_possible <- sum(l(set1), l(set2)) * .i
  l(x[which(x <= max_possible)]) / max_possible}
