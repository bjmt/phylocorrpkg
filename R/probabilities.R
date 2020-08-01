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

#' @export
getProbMatchFunDouble <- function(df1, df2, highScoreIsBest1 = TRUE,
  highScoreIsBest2 = TRUE, maxProb1 = 1, maxProb1 = 1) {

}

#' @export
calcMatchingProbsSingle <- function(scores, predFUN) {
  scores <- as.matrix(scores)
  rownames(scores) <- colnames(scores)
  diag(scores) <- NA
  matrix(predFUN(scores), nrow = nrow(scores), dimnames = dimnames(scores))
}
