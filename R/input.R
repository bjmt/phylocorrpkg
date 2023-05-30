#' Input preparation
#'
#' These functions help in cleaning up and preparing the input presence-absence
#' matrix and the annotations.
#'
#' @param x A `matrix` or `data.frame` with family/gene/protein IDs as column
#'    names, and tree tip names as row names. The contents of the matrix should be
#'    logical presence-absence data, though integer/numeric data can be used and will
#'    be converted to logical.
#' @param tree A [ape::phylo-class] type tree object. The tip labels must be identical
#'    to the row names of the input matrix `x`.
#' @param minFamilyCount The minimum number of presence calls a family/gene/protein
#'    must have in order not to be filtered out. This helps greatly with increasing
#'    function prediction accuracy.
#' @param maxFamilyCount The maximum number of presence calls a family/gene/protein
#'    can have in order not to be filtered out. This helps remove those which are
#'    present in the entire tree, and thus do not contribute to function prediction
#'    accuracy.
#' @param annotations A `data.frame` with two columns. The first is a character vector
#'    of family/gene/protein IDs (matching the column names of the input `x`), and
#'    the second column is a character vector of annotations of any type. One
#'    family/gene/protein can have multiple annotations, present as additional rows.
#' @param annoIsSecondCol If the columns are switched then set this to `FALSE`.
#'
#' @return
#'    For [cleanData()], a presence-absence logical matrix.
#'
#'    For [cleanAnnotations()], a list of family/gene/protein IDs to annotations.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name InputFunctions

#' @rdname InputFunctions
#' @export
cleanData <- function(x, tree = NULL, minFamilyCount = 2, maxFamilyCount = nrow(x) - 1) {
  if (!is.matrix(x)) stop("Input must be a matrix")
  x <- as.matrix(x)
  if (!is.numeric(x) && !is.logical(x)) stop("Input must be a numeric or logical matrix")
  if (is.null(rownames(x)) || is.null(colnames(x))) stop("Column and row names are required.")
  x <- matrix(as.logical(x), nrow = nrow(x), dimnames = dimnames(x))
  if (!is.null(tree)) {
    if (nrow(x) != length(tree$tip.label)) warning("Input row count does not match tree tip count")
    if (any(!rownames(x) %in% tree$tip.label)) stop("Mismatching tree tip labels and input row names")
    x <- x[tree$tip.label, ]
  }
  famSums <- colSums(x)
  # genomeSums <- rowSums(x)
  if (any(famSums < minFamilyCount)) {
    message("Found families with fewer calls than the allowed minimum, these will be removed.")
    x <- x[, famSums >= minFamilyCount]
    famSums <- colSums(x)
  }
  if (any(famSums > maxFamilyCount)) {
    message("Found families with more calls than the allowed maximum, these will be removed.")
    x <- x[, famSums <= maxFamilyCount]
  }
  x
}

#' @rdname InputFunctions
#' @export
cleanAnnotations <- function(annotations, annoIsSecondCol = TRUE) {
  if (annoIsSecondCol)
    lapply(split(annotations, annotations[[1]]), function(x) x[[2]])
  else
    lapply(split(annotations, annotations[[2]]), function(x) x[[1]])
}
