# Phylocorrelate -- the R package

The functions within this package have been provided so that one may make use of the methods which were used to generate the Phylocorrelate database.

Phylocorrelate is a tool that detects pairs of gene/protein families with similar phylogenetic distributions. Gene co-occurrence is important to detect as it may suggest an interaction between the genes, membership in the same pathway or complex, or more general functional associations or dependencies.

The online web-server ([PhyloCorrelate](https://phylocorrelate.uwaterloo.ca)) contains precomputed correlations for 27,372 gene families across 28,315 species, and uses the GTDB bacterial tree and gene occurrences from AnnoTree including entries from KEGG, PFAM, and TIGRFAM.

This R package extends the functionality of PhyloCorrelate, so that it can be used with any custom tree and dataset of genes (binary presence/absence matrix).

## Citation

If this package proves useful in your research, please cite:

Tremblay, B.J.M., Lobb, B., and Doxey, A.C. (2021). PhyloCorrelate: inferring bacterial gene-gene functional associations through large-scale phylogenetic profiling. _Bioinformatics_, 37(1), 17-22. DOI:[10.1093/bioinformatics/btaa1105](https://doi.org/10.1093/bioinformatics/btaa1105).

## Installation

```r
install.packages(c("remotes", "BiocManager"))
BiocManager::install("bjmt/phylocorrpkg")
```

Install the extra packages used in this guide:

```r
install.packages(c("ape", "doSNOW", "fst", "plot.matrix"))
```

The functions within the package are divided into three categories. To see the manual pages for each:

```r
library(phylocorrpkg)

?InputFunctions
?CalculationFunctions
?ProbabilityFunctions
```

## Data preparation

The `fst` format allows individual columns to be accessed very quickly, making it the optimal storage solution. However this means saving data as a `data.frames`, as well as saving the row names separately (if needed; the score matrices have identical column and row names, making this unnecessary).

To start, two things are needed: a phylogenetic tree, and a presence-absence matrix with tree tips as rows and families/genes/etc as column names. The presence-absence matrix can also be a matrix of counts, but all values greater than zero will be reduced to `TRUE`.

```r
library(phylocorrpkg)
library(fst)
library(ape)

tree <- read.tree("tree")
input <- read.table("input", header = TRUE, quote = "")

input <- cleanData(input, tree)

writeLines(rownames(input), "genomes.txt")
fst::write_fst(as.data.frame(input), "Table.fst")
```

If you wish, you can also use the input tables used to generate the PhyloCorrelate data. They can be downloaded [here](https://zenodo.org/record/3993422). The data preparation stage has already been performed on these files, and you can use them directly in the next section. The GO:BP/Pathway annotation data can also be downloaded and used in the annotation cleaning and probability matching section. It is recommended to use the TIGRFAM dataset, as the PFAM and KO datasets are much larger and will require much more intensive computation.

## Score calculation

Depending on the number of families, these calculations can consume a lot of RAM and CPU time. Parallelisation can be achieved using the `doParCalc()` function, though this is optional; the function will still work even if no cluster is registered. In order not to run out of RAM, it is recommended that objects are immediately written to disk and deleted from memory. There are several calculations available:

- Pearson correlation coefficient (PCC): `calcPCC()`
- Runs-adjusted PCC (rPCC): `calcrPCCPair()`
- Jaccard Coefficient (JC): `calcJCPair()`
- Runs-adjusted JC (rJC): `calcrJCPair()`
- Hypergeometric P-value (HyperP): `calcHyperPPair()`
- Runs-adjusted HyperP (rHyperP): `calcrHyperPPair()`
- Overlap between pairs (Ov): `calcOvPair()`
- Difference of occurrences between pairs (OccDiff): `calcOccDiff()`

You can also create your own comparison function: all it needs to do is be able to take two logical vectors and generate a single numerical score. You can use the `doParCalc()` function to apply this function to all pairs in the input matrix.

```r
library(phylocorrpkg)
library(fst)
library(foreach)
library(doSNOW)

cl <- makeCluster(8)
registerDoSNOW(cl)

input <- read_fst("Table.fst")

# calcPCC() is merely a wrapper for the stats::cor function which can be
# used directly without any parallelisation, as it is quite efficient
PCC <- calcPCC(input)
write_fst(as.data.frame(PCC), "PCC.fst")
rm(PCC)

rPCC <- doParCalc(input, calcrPCCPair)
write_fst(as.data.frame(rPCC), "rPCC.fst")
rm(rPCC)

JC <- doParCalc(input, calcJCPair)
write_fst(as.data.frame(JC), "JC.fst")
rm(JC)

rJC <- doParCalc(input, calcrJCPair)
write_fst(as.data.frame(rJC), "rJC.fst")
rm(rJC)

HyperP <- doParCalc(input, calcHyperPPair)
write_fst(as.data.frame(HyperP), "HyperP.fst")
rm(HyperP)

rHyperP <- doParCalc(input, calcrHyperPPair)
write_fst(as.data.frame(rHyperP), "rHyperP.fst")
rm(rHyperP)

Ov <- doParCalc(input, calcOvPair)
write_fst(as.data.frame(Ov), "Ov.fst")
rm(Ov)

# For calcOccDiff() parallelisation is not required, as this is a simple
# calculation
OccDiff <- calcOccDiff(input)
write_fst(as.data.frame(OccDiff), "OccDiff.fst")
rm(OccDiff)
```

## Matching annotation probability functions from single scores

This can be done with any annotation type. All that's needed is a `data.frame` with two columns: the first one with family names (matching those of the input column names), and the second with the annotation of interest. Families can have multiple annotations, appearing as duplicate rows. The probability functions can be obtained from any of the scores previously calculated. This probability prediction function can then be used to generate probabilities for the entire input score set.

In this example workflow, function match probabilities will be calculated using PCC scores.

```r
library(phylocorrpkg)
library(fst)

annotations <- read.table("annotations", header = TRUE, quote = "")
annotations <- cleanAnnotations(annotations)

PCC <- read_fst("PCC.fst")

# Only scores which are associated with annotations are kept to generate the
# prediction function
PCCdf <- mergeScoresAndAnnotations(PCC, annotations)
PCCfun <- getProbMatchFunSingle(PCCdf)

x <- seq(from = -1, to = 1, by = 0.01)
plot(x, PCCfun(x), type = "b", xlab = "PCC",
    ylab = "Probability of a Matching Function")

PCCPredictedProbs <- calcMatchingProbsSingle(PCC, PCCfun)
```

If there too few data points near the upper range of scores, the probability prediction function can behave unexpectedly for high scores. Using smoothing can help greatly with this (by setting `useMeanSmoothing = TRUE` in `getProbMatchFunSingle()`, and adjusting the `windowSize` parameter if needed). Additionally, if none of the probabilities from actual scores approach 1 then setting a lower `maxProb` can help generate more realistic probabilities for higher scores. Visually inspecting the Score VS Probability plot as above can greatly aid in determining the best parameters to use in the `getProbMatchFunSingle()` function.

## Matching annotation probability functions from two sets of scores

To increase the prediction accuracy, the probability prediction function can be generated from two different scores. Instead of a probability prediction function, this returns a 2D matrix of probabilities.

Note: for P-value type scores, it is recommend to transform them with `-log10()` before generating the probability prediction function. Also make sure to change any infinite values (resulting from P-values equalling zero) to a suitable high number.

In this example workflow, function match probabilities will be calculated from rHyperP and OccDiff scores.

```r
library(phylocorrpkg)
library(fst)
library(plot.matrix)

annotations <- read.table("annotations", header = TRUE, quote = "")
annotations <- cleanAnnotations(annotations)

rHyperP <- as.matrix(read_fst("rHyperP.fst"))
rHyperP <- matrix(-log10(rHyperP), nrow = nrow(rHyperP),
    dimnames = list(colnames(rHyperP), colnames(rHyperP)))
rHyperP[is.infinite(rHyperP)] <- 321
OccDiff <- read_fst("OccDiff.fst")

rHyperPdf <- mergeScoresAndAnnotations(rHyperP, annotations)
OccDiffdf <- mergeScoresAndAnnotations(OccDiff, annotations)

CombinedDf <- getProbMatchDfDouble(rHyperPdf, OccDiffdf)
CombinedMatrix <- getProbMatchMatrixDouble(CombinedDf)

plot(CombinedMatrix[nrow(CombinedMatrix):1, ], breaks = 10,
    xlab = "OccDiff", ylab = "rHyperP")

CombinedPredicted <- calcMatchingProbsDouble(rHyperP, OccDiff, CombinedMatrix)
```

There are two opportunities for smoothing when generating the final probability prediction matrix: in `getProbMatchDfDouble()` and `getProbMatchMatrixDouble()`. The former function is a 2D kernel function; here you can decide the size of the final matrix (`len1`, `len2`) and the bandwidth values. In the latter function, an optional mean smoothing step can be used (`useMeanSmoothing = TRUE`) controlled by `windowSize1` and `windowSize2`. Visually inspecting the resulting matrix with the above code can greatly help here.
