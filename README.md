# Phylocorrelate -- the R package

The functions within this package have been provided so that one may make use of the methods used to generate the Phylocorrelate database.

## Data preparation

The `fst` format allows individual columns to be accessed very quickly, making it the optimal storage solution. However this means saving data as a `data.frames`, as well as saving the row names separately.


```r
library(phylocorrpkg)
library(fst)
library(ape)
library(readr)

tree <- read.tree("tree")
input <- read.table("input", header = TRUE, quote = "")

input <- cleanData(input, tree)

write_lines(rownames(input), "genomes.txt")
fst::write_fst(as.data.frame(input), "Table.fst")
```

## Score calculation

Depending on the number of families, these calculations can consume a lot of RAM and CPU time. Parallelisation can be achieved using the `doParCalc()` function. In order not to run out of RAM, objects are immediately written to disk and deleted from memory.

```r
library(phylocorrpkg)
library(fst)
library(foreach)
library(doSNOW)

cl <- makeCluster(8)
registerDoSNOW(cl)

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

OccDiff <- calcOccDiff(input)
write_fst(as.data.frame(OccDiff), "OccDiff.fst")
rm(OccDiff)
```

## Matching annotation probability functions from single scores

This can be done with any annotation type. All that's needed is a `data.frame` with two columns: the first one with family names (matching those of the input column names), and the second with the annotation of interest. Families can have multiple annotations, appearing as duplicate rows. The probability functions can be obtained from any of the scores previously calculated.

```r
library(phylocorrpkg)
library(fst)

annotations <- read.table("annotations", header = TRUE, quote = "")
annotations <- cleanAnnotations(annotations)

PCC <- read_fst("PCC.fst")

PCCdf <- mergeScoresAndAnnotations(PCC, annotations)
PCCfun <- getProbMatchFunSingle(PCCdf)

x <- seq(from = -1, to = 1, by = 0.01)
plot(x, PCCfun(x), type = "b", xlab = "PCC",
    ylab = "Probability of a Matching Function")

PCCPredictedProbs <- calcMatchingProbsSingle(PCC, PCCfun)
```

If there too few data points near the upper range of scores, the probability prediction function can behave unexpectedly for high scores. Using smoothing can help greatly with this (by setting `useMeanSmoothing = TRUE` in `getProbMatchFunSingle()`, and adjusting the `windowSize` parameter if needed). Additionally, if none of the probabilities from actual scores approach 1 then setting a lower `maxProb` can help generate more realistic probabilities for higher scores. Visually inspecting the Score VS Probability plot as above can greatly aid in determining the best parameters to use in the `getProbMatchFunSingle()` function.

Note: for P-value type scores, it is recommend to transform them with `-log10()` before generating the probability prediction function.

