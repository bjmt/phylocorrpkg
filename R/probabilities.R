#' Probability prediction functions.
#'
#' Functions to be used to generate the prediction functions/matrices.
#'
#' @param scores
#' @param annotations
#' @param highScoreIsBest
#' @param merged
#' @param maxProb
#' @param useMeanSmoothing
#' @param windowSize
#' @param df1
#' @param df2
#' @param len1
#' @param len2
#' @param lims1
#' @param lims2
#' @param bw1Match
#' @param bw1NoMatch
#' @param bw2Match
#' @param bw2NoMatch
#' @param CombinedDf
#' @param windowSize1
#' @param windowSize2
#' @param scores1
#' @param scores2
#' @param predFUN
#' @param predMatrix
#'
#' @return
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name ProbabilityFunctions

#' @rdname ProbabilityFunctions
#' @export
mergeScoresAndAnnotations <- function(scores, annotations, highScoreIsBest = TRUE) {
  scores <- as.matrix(scores)
  rownames(scores) <- colnames(scores)
  scores <- scores[rownames(scores) %in% names(annotations),
      colnames(scores) %in% names(annotations)]
  scores[upper.tri(scores)] <- NA
  diag(scores) <- NA
  scores2 <- melt(scores, na.rm = TRUE)
  colnames(scores2) <- c("Fam1", "Fam2", "Score")
  scores2$Anno1 <- annotations[scores2$Fam1]
  scores2$Anno2 <- annotations[scores2$Fam2]
  scores2$MatchingAnno <- mapply(function(x, y) any(x %in% y),
      scores2$Anno1, scores2$Anno2)
  scores2
}

#' @rdname ProbabilityFunctions
#' @export
getProbMatchFunSingle <- function(merged, highScoreIsBest = TRUE, maxProb = 1,
  useMeanSmoothing = FALSE, windowSize = max(1, round(0.0001 * nrow(merged)))) {
  if (anyNA(merged$Score))
    stop("Found NA values in Score column")
  if (any(is.infinite(merged$Score)))
    stop("Found non-finite values in score column")
  merged <- merged[order(merged$Score, decreasing = highScoreIsBest), ]
  merged$CumProbMatch <- cumsum(merged$MatchingAnno) / 1:nrow(merged)
  ProbMatch <- merged$CumProbMatch
  Scores <- merged$Score
  if (useMeanSmoothing) {
    ProbMatch <- as.vector(stats::filter(ProbMatch, rep(1 / windowSize, windowSize),
      sides = 2))
    Scores <- Scores[!is.na(ProbMatch)]
    ProbMatch <- ProbMatch[!is.na(ProbMatch)]
  }
  if (highScoreIsBest)
    approxfun(Scores, ProbMatch, yleft = 0, yright = maxProb)
  else
    approxfun(Scores, ProbMatch, yleft = maxProb, yright = 0)
}

#' @rdname ProbabilityFunctions
#' @export
getProbMatchDfDouble <- function(df1, df2, len1 = 21, len2 = 21,
  lims1 = range(df1$Score), lims2 = range(df2$Score),
  bw1Match = MASS::bandwidth.nrd(df1$Score[df1$MatchingAnno]) * 2,
  bw1NoMatch = MASS::bandwidth.nrd(df1$Score[!df1$MatchingAnno]) * 2,
  bw2Match = MASS::bandwidth.nrd(df2$Score[df2$MatchingAnno]) * 2,
  bw2NoMatch = MASS::bandwidth.nrd(df2$Score[!df2$MatchingAnno]) * 2) {

  df1$Score[df1$Score > lims1[2]] <- lims1[2]
  df1$Score[df1$Score < lims1[1]] <- lims1[1]
  df2$Score[df2$Score > lims2[2]] <- lims2[2]
  df2$Score[df2$Score < lims2[1]] <- lims2[1]

  dfTrue <- as.matrix(cbind(df1$Score[df1$MatchingAnno],
      df2$Score[df2$MatchingAnno]))
  dfFalse <- as.matrix(cbind(df1$Score[!df1$MatchingAnno],
      df2$Score[!df2$MatchingAnno]))

  seq1 <- seq(from = lims1[1], to = lims1[2], length.out = len1)
  seq2 <- seq(from = lims2[1], to = lims2[2], length.out = len2)

  binsTrue1 <- cut(dfTrue[, 1], seq1, include.lowest = TRUE)
  binsTrue2 <- cut(dfTrue[, 2], seq2, include.lowest = TRUE)
  binsFalse1 <- cut(dfFalse[, 1], seq1, include.lowest = TRUE)
  binsFalse2 <- cut(dfFalse[, 2], seq2, include.lowest = TRUE)

  kernTrue <- bkde2D(dfTrue, c(bw1Match, bw2Match), c(len1, len2),
    list(lims1, lims2))
  kernFalse <- bkde2D(dfTrue, c(bw1NoMatch, bw2NoMatch), c(len1, len2),
    list(lims1, lims2))

  kernTrueNorm <- (kernTrue$fhat / sum(kernTrue$fhat)) * nrow(dfTrue)
  kernFalseNorm <- (kernFalse$fhat / sum(kernFalse$fhat)) * nrow(dfFalse)

  kernRatio <- matrix(kernTrueNorm / (kernFalseNorm + kernTrueNorm),
    nrow(kernTrue$fhat), dimnames = list(seq1, seq2))

  dfPMF <- expand.grid(x = seq1, y = seq2)

  xInd <- match(as.character(dfPMF$x), rownames(kernRatio)) - 1
  yInd <- match(as.character(dfPMF$y), colnames(kernRatio)) - 1

  kernRatio[is.na(kernRatio)] <- 0
  kernRatioVec <- as.vector(kernRatio)
  kernRatioVec[is.na(kernRatioVec)] <- 0
  kernTrueNormVec <- as.vector(kernTrueNorm)
  kernFalseNormVec <- as.vector(kernFalseNorm)

  TFcounts <- calcCounts(xInd, yInd, kernRatio, kernRatioVec, kernTrueNormVec,
    kernFalseNormVec)

  dfPMF$CountTrue <- TFcounts$Tcount
  dfPMF$CountFalse <- TFcounts$Fcount
  dfPMF$CountTotal <- dfPMF$CountTrue + dfPMF$CountFalse
  dfPMF$MatchProb <- dfPMF$CountTrue / dfPMF$CountTotal

  dfPMF <- dfPMF[order(dfPMF$CountTotal, dfPMF$MatchProb), ]

}

#' @rdname ProbabilityFunctions
#' @export
getProbMatchMatrixDouble <- function(CombinedDf, useMeanSmoothing = FALSE,
  windowSize1 = 3, windowSize2 = 3) {
  CombinedMatrix <- acast(CombinedDf, x ~ y, value.var = "MatchProb")
  if (useMeanSmoothing) {
    CombinedRaster <- raster(CombinedMatrix)
    CombinedSmoothed <- as.matrix(focal(CombinedRaster,
        matrix(1, windowSize1, windowSize2), mean, pad = FALSE))
    dimnames(CombinedSmoothed) <- dimnames(CombinedMatrix)
    CombinedMatrix <- CombinedSmoothed
  }
  CombinedMatrix
}

#' @rdname ProbabilityFunctions
#' @export
calcMatchingProbsSingle <- function(scores, predFUN) {
  scores <- as.matrix(scores)
  rownames(scores) <- colnames(scores)
  diag(scores) <- NA
  matrix(predFUN(scores), nrow = nrow(scores), dimnames = dimnames(scores))
}

#' @rdname ProbabilityFunctions
#' @export
calcMatchingProbsDouble <- function(scores1, scores2, predMatrix) {
  cuts1 <- as.numeric(rownames(predMatrix))
  cuts1 <- c(-Inf, cuts1[-1] - diff(cuts1) / 2, Inf)
  cuts2 <- as.numeric(colnames(predMatrix))
  cuts2 <- c(-Inf, cuts2[-1] - diff(cuts2) / 2, Inf)
  scores1 <- as.matrix(scores1)
  scores2 <- as.matrix(scores2)
  rownames(scores1) <- colnames(scores1)
  rownames(scores2) <- colnames(scores2)
  scores1[upper.tri(scores1)] <- NA
  diag(scores1) <- NA
  scores2[upper.tri(scores2)] <- NA
  diag(scores2) <- NA
  scores1 <- melt(scores1, na.rm = TRUE)
  scores2 <- melt(scores2, na.rm = TRUE)
  ind1 <- findInterval(scores1[[3]], cuts1) - 1
  ind2 <- findInterval(scores2[[3]], cuts2) - 1
  data.frame(Fam1 = scores1[[1]], Fam2 = scores1[[2]],
    PMF = calcPMF2d(ind1, ind2, predMatrix))
}
