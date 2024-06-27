
diversity_index <- \(set1, set2, .i = 0.5) {
  if (base::length(set1) != base::length(set2)) {message("Set size is not equal!")}
  x <- cumsum(base::sort(table(c(set1, set2)), decreasing = T))
  max_possible <- sum(base::length(set1), base::length(set2)) * .i
  return(base::length(x[which(x <= max_possible)]) / max_possible)}
