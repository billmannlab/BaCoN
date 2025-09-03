#' @export import_depmap

import_depmap <- function(filepath) {
  object <- as.matrix(data.table::fread(filepath), rownames = 1)
  colnames(object) <- stringr::str_split_i(colnames(object), " \\(", 1)
  return(object)
}
