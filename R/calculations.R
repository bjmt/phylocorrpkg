#' Metric calculation functions.
#'
#' Comparison functions to calculate the various metrics on the presence-absence
#' input matrix.
#'
#' @param x A single column from the presence-absence input matrix to compare
#'    to `y`.
#' @param y The second single column with which to compare to `x`.
#' @param z The presence-absence input matrix.
#' @param FUN The function perform the matrix-wide all-against-all
#'    calculations, optionally in parallel.
#'
#' @return
#'    For [getRuns()]: a list of `x` and `y`, runs-adjusted.
#'
#'    For all other functions: either a single value if comparing two columns,
#'    or a matrix of values if comparing all-by-all.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name CalculationFunctions

#' @rdname CalculationFunctions
#' @export
getRuns <- function(x, y) {
  z <- integer(length(x))
  z[x & y] <- 1
  z[x & !y] <- 2
  z[!x & y] <- 3
  z <- Rle(z)
  z <- runValue(z)
  x <- y <- logical(length(z))
  x[z %in% 1:2] <- TRUE
  y[z %in% c(1, 3)] <- TRUE
  list(x = x, y = y)
}

#' @rdname CalculationFunctions
#' @export
doParCalc <- function(z, FUN) {
  y <- foreach(i = seq_len(ncol(z))) %:%
    foreach(j = seq_len(ncol(z))[-seq_len(i)], .combine = "c") %dopar% {
      FUN(z[, i], z[, j])
    }
  for (i in seq_along(y)) {
    y[[i]] <- c(rep(NA_real_, ncol(z) - length(y[[i]])), y[[i]])
  }
  y <- matrix(unlist(y), ncol = ncol(z), dimnames = list(colnames(z), colnames(z)))
  y <- as.matrix(as.dist(y))
  diag(y) <- NA
  y
}

#' @rdname CalculationFunctions
#' @export
calcPCC <- function(z) {
  cor(z)
}

#' @rdname CalculationFunctions
#' @export
calcPCCPair <- function(x, y) {
  cor(x, y)
}

#' @rdname CalculationFunctions
#' @export
calcrPCCPair <- function(x, y) {
  z <- getRuns(x, y)
  cor(z$x, z$y)
}

#' @rdname CalculationFunctions
#' @export
calcHyperPPair <- function(x, y) {
  fisher.test(
    matrix(c(
      sum(x & y), sum(y & !x),
      sum(!y & x), sum(!y & !x)
    ), ncol = 2),
    alternative = "g"
  )$p.value
}

#' @rdname CalculationFunctions
#' @export
calcrHyperPPair <- function(x, y) {
  z <- getRuns(x, y)
  calcHyperPPair(z$x, z$y)
}

#' @rdname CalculationFunctions
#' @export
calcOvPair <- function(x, y) {
  sum(x & y)
}

#' @rdname CalculationFunctions
#' @export
calcJCPair <- function(x, y) {
  sum(x & y) / sum(x | y)
}

#' @rdname CalculationFunctions
#' @export
calcrJCPair <- function(x, y) {
  z <- getRuns(x, y)
  calcJCPair(z$x, z$y)
}

#' @rdname CalculationFunctions
#' @export
calcOccDiff <- function(z) {
  z <- colSums(z)
  z <- outer(z, z, FUN = function(x, y) abs(x - y))
  dimnames(z) <- list(colnames(z), colnames(z))
  z
}

#' @rdname CalculationFunctions
#' @export
calcOccDiffPair <- function(x, y) {
  abs(sum(x) - sum(y))
}
