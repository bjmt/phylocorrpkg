#' @export
cleanData <- function(x, tree = NULL, minFamilyCount = 2, maxFamilyCount = nrow(x) - 1) {
  if (!is.matrix(x)) stop("Input must be a matrix")
  if (!is.numeric(x) && !is.logical(x)) stop("Input must be a numeric or logical matrix")
  if (is.null(rownames(x)) || is.null(colnames(x))) stop("Column and row names are required.")
  x <- matrix(as.logical(x), nrow = nrow(x), dimnames = dimnames(x))
  if (!is.null(tree)) {
    if (nrow(x) != length(tree$tip.label)) warning("Input row count does not match tree tip count")
    if (any(!rownames(x) %in% tree$tip.label)) stop("Mismatching tree tip labels and input row names")
    x <- x[tree$tip.label, ]
  }
  famSums <- colSums(x)
  genomeSums <- rowSums(x)
  if (any(famSums < minFamilyCount)) {
    message("Found families with fewer calls than the allowed minimum, these will be removed.")
    x <- x[, famSums >= minFamilyCount]
  }
  if (any(famSums > maxFamilyCount)) {
    message("Found families with more calls than the allowed maximum, these will be removed.")
    x <- x[, famSums <= maxFamilyCount]
  }
  x
}

#' @export
cleanAnnotations <- function(annotations, annoIsSecondCol = TRUE) {
  if (annoIsSecondCol)
    lapply(split(annotations, annotations[[1]]), function(x) x[[2]])
  else
    lapply(split(annotations, annotations[[2]]), function(x) x[[1]])
}
