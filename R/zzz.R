#' phylocorrpkg: Phylogenetic correlations across large trees
#'
#' @description
#' Functional predictions from correlations across large phylogenetic
#' trees using presence-absence data.
#'
#' @docType package
#' @name phylocorrpkg-pkg
#'
#' @importFrom stats fisher.test cor as.dist approxfun
#' @importFrom S4Vectors Rle runValue
#' @importFrom reshape2 melt acast
#' @importFrom foreach foreach %:% %dopar%
#' @importFrom KernSmooth bkde2D
#' @importFrom raster focal raster
#' @useDynLib phylocorrpkg
NULL
