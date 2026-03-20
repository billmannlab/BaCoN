baconize <- function(.vec, .cf) {
  .out <- rep(NA, length(.vec))
  i <- .vec >= 0 & !is.na(.vec)
  .out[i] <- purrr::map_int(.vec[i], ~ sum(.vec > .x - .cf, na.rm = T))
  if (sum(!i) > 0) {
    .out[!i] <- purrr::map_int(.vec[!i], ~ sum(.vec < .x + .cf, na.rm = T))}
  .out}
