#' @export sort_gene_pairs

sort_gene_pairs <- function(g1, g2, sep = ";", pairs, pair_sep = ";", invert = F) {
  
  stopifnot("invert needs to be logical." = is.logical(invert), 
            "No input pairs given!" = !(missing(pairs) & missing(g1) & missing(g2)), 
            "g1 and g2 have to be of equal length!" = (!missing(g1) & !missing(g2) & length(g1) == length(g2)), 
            "Either two gene vectors or a gene pair vector required" = !(missing(pairs) & (missing(g1) | missing(g2))))
  
  if (!missing(pairs) & missing(g1) & missing(g2)) {
    g1 <- stringr::str_split_i(pairs, pattern = pair_sep, 1)
    g2 <- stringr::str_split_i(pairs, pattern = pair_sep, 2)
  }
  
  return(switch(as.character(invert),
                "TRUE" = paste(pmax(g1, g2), pmin(g1, g2), sep = sep), 
                "FALSE" = paste(pmin(g1, g2), pmax(g1, g2), sep = sep)))
  
}
