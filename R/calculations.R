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

doParCalc <- function(x, FUN) {
  y <- foreach(i = seq_len(ncol(x))) %:%
    foreach(j = seq_len(ncol(x))[-seq_len(i)], .combine = "c") %dopar% {
      FUN(x[, i], x[, j])
    }
  for (i in seq_along(y)) {
    y[[i]] <- c(rep(NA_real_, ncol(x) - length(y[[i]])), y[[i]])
  }
  y <- matrix(unlist(y), ncol = ncol(x), dimnames = list(colnames(x), colnames(x)))
  y <- as.matrix(as.dist(y))
  diag(y) <- NA
  y
}

calcPCC <- function(x) {
  cor(x)
}

calcPCCPair <- function(x, y) {
  cor(x, y)
}

calcrPCCPair <- function(x, y) {
  z <- getRuns(x, y)
  cor(z$x, z$y)
}

calcHyperPPair <- function(x, y) {
  fisher.test(
    matrix(c(
      sum(x & y), sum(y & !x),
      sum(!y & x), sum(!y & !x)
    ), ncol = 2),
    alternative = "g"
  )$p.value
}

calcrHyperPPair <- function(x, y) {
  z <- getRuns(x, y)
  calcHyperP(z$x, z$y)
}

calcOvPair <- function(x, y) {
  sum(x & y)
}

calcJCPair <- function(x, y) {
  sum(x & y) / sum(x | y)
}

calcrJCPair <- function(x, y) {
  z <- getRuns(x, y)
  calcJC(z$x, z$y)
}

calcOccDiff <- function(x) {
  z <- colSums(x)
  z <- outer(z, z, FUN = function(x, y) abs(x - y))
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

calcOccDiffPair <- function(x, y) {
  abs(sum(x) - sum(y))
}
