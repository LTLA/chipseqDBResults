---
title: Detecting differential enrichment of H3K9ac in murine B cells
author:
- name: Aaron T. L. Lun
  affiliation: 
  - &WEHI The Walter and Eliza Hall Institute of Medical Research, 1G Royal Parade, Parkville, VIC 3052, Melbourne, Australia
  - Department of Medical Biology, The University of Melbourne, Parkville, VIC 3010, Melbourne, Australia
- name: Gordon K. Smyth
  affiliation: 
  - *WEHI
  - Department of Mathematics and Statistics, The University of Melbourne, Parkville, VIC 3010, Melbourne, Australia
date: "2019-05-27"
vignette: >
  %\VignetteIndexEntry{2. Differential enrichment of H3K9ac in B cells}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
output:
  BiocStyle::html_document:
    toc_float: yes
    titlecaps: false
bibliography: ref.bib
---



# Overview

Here, we perform a window-based differential binding (DB) analysis to identify regions of differential H3K9ac enrichment between pro-B and mature B cells [@domingo2012bcell].
H3K9ac is associated with active promoters and tends to exhibit relatively narrow regions of enrichment relative to other marks such as H3K27me3.
We download the BAM files using the relevant function from the *[chipseqDBData](https://bioconductor.org/packages/3.10/chipseqDBData)* package^[BAM files are cached upon the first call to this function, so subsequent calls do not need to re-download the files.].
The experimental design contains two biological replicates for each of the two cell types.


```r
library(chipseqDBData)
acdata <- H3K9acData()
acdata
```

```
## DataFrame with 4 rows and 3 columns
##                  Name            Description
##           <character>            <character>
## 1    h3k9ac-proB-8113    pro-B H3K9ac (8113)
## 2    h3k9ac-proB-8108    pro-B H3K9ac (8108)
## 3 h3k9ac-matureB-8059 mature B H3K9ac (8059)
## 4 h3k9ac-matureB-8086 mature B H3K9ac (8086)
##                                                      Path
##                                               <character>
## 1    /tmp/RtmpnLp3hN/file36c3d36d627/h3k9ac-proB-8113.bam
## 2    /tmp/RtmpnLp3hN/file36c3d36d627/h3k9ac-proB-8108.bam
## 3 /tmp/RtmpnLp3hN/file36c3d36d627/h3k9ac-matureB-8059.bam
## 4 /tmp/RtmpnLp3hN/file36c3d36d627/h3k9ac-matureB-8086.bam
```

# Pre-processing checks 

## Examining mapping statistics

We use methods from the *[Rsamtools](https://bioconductor.org/packages/3.10/Rsamtools)* package to compute some mapping statistics for each BAM file.
Ideally, the proportion of mapped reads should be high (70-80% or higher), 
while the proportion of marked reads should be low (generally below 20%).


```r
library(Rsamtools)
diagnostics <- list()
for (bam in acdata$Path) {
    total <- countBam(bam)$records
    mapped <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE)))$records
    marked <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE, isDuplicate=TRUE)))$records
    diagnostics[[basename(bam)]] <- c(Total=total, Mapped=mapped, Marked=marked)
}
diag.stats <- data.frame(do.call(rbind, diagnostics))
diag.stats$Prop.mapped <- diag.stats$Mapped/diag.stats$Total*100
diag.stats$Prop.marked <- diag.stats$Marked/diag.stats$Mapped*100
diag.stats
```

```
##                            Total  Mapped Marked Prop.mapped Prop.marked
## h3k9ac-proB-8113.bam    10724526 8832006 434884    82.35335    4.923955
## h3k9ac-proB-8108.bam    10413135 7793913 252271    74.84694    3.236770
## h3k9ac-matureB-8059.bam 16675372 4670364 396785    28.00756    8.495805
## h3k9ac-matureB-8086.bam  6347683 4551692 141583    71.70635    3.110558
```

Note that all *[csaw](https://bioconductor.org/packages/3.10/csaw)* functions that read from a BAM file require BAM indices with `.bai` suffixes.
In this case, index files have already been downloaded by `H3K9acData()`, 
but users supplying their own files should take care to ensure that BAM indices are available with appropriate names.

## Obtaining the ENCODE blacklist for mm10

A number of genomic regions contain high artifactual signal in ChIP-seq experiments.
These often correspond to genomic features like telomeres or microsatellite repeats.
For example, multiple tandem repeats in the real genome are reported as a single unit in the genome build.
Alignment of all (non-specifically immunoprecipitated) reads from the former will result in artificially high coverage of the latter.
Moreover, differences in repeat copy numbers between conditions can lead to detection of spurious DB.

As such, these problematic regions must be removed prior to further analysis.
This is done with an annotated blacklist for the mm10 build of the mouse genome, 
constructed by identifying consistently problematic regions from ENCODE datasets [@encode2012encode].
We download this BED file and save it into a local cache with the *[BiocFileCache](https://bioconductor.org/packages/3.10/BiocFileCache)* package.
This allows it to be used again in later workflows without being re-downloaded.


```r
library(BiocFileCache)
bfc <- BiocFileCache("local", ask=FALSE)
black.path <- bfcrpath(bfc, file.path("https://www.encodeproject.org",
    "files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"))
```

Genomic intervals in the blacklist are loaded using the `import()` method from the *[rtracklayer](https://bioconductor.org/packages/3.10/rtracklayer)* package.
All reads mapped within the blacklisted intervals will be ignored during processing in *[csaw](https://bioconductor.org/packages/3.10/csaw)* by specifying the `discard=` parameter (see below).


```r
library(rtracklayer)
blacklist <- import(black.path)
blacklist
```

```
## GRanges object with 164 ranges and 0 metadata columns:
##         seqnames              ranges strand
##            <Rle>           <IRanges>  <Rle>
##     [1]    chr10     3110061-3110270      *
##     [2]    chr10   22142531-22142880      *
##     [3]    chr10   22142831-22143070      *
##     [4]    chr10   58223871-58224100      *
##     [5]    chr10   58225261-58225500      *
##     ...      ...                 ...    ...
##   [160]     chr9     3038051-3038300      *
##   [161]     chr9   24541941-24542200      *
##   [162]     chr9   35305121-35305620      *
##   [163]     chr9 110281191-110281400      *
##   [164]     chr9 123872951-123873160      *
##   -------
##   seqinfo: 19 sequences from an unspecified genome; no seqlengths
```

Any user-defined set of regions can be used as a blacklist in this analysis.

- For example, one could use predicted repeat regions from the UCSC genome annotation [@rosenbloom2015ucsc].
This tends to remove a greater number of problematic regions (especially microsatellites) compared to the ENCODE blacklist.
However, the size of the UCSC list means that genuine DB sites may also be removed.
Thus, the ENCODE blacklist is preferred for most applications.
- Alternatively, if negative control samples are available, they can be used to empirically identify problematic regions with the *[GreyListChIP](https://bioconductor.org/packages/3.10/GreyListChIP)* package.
These regions should be ignored as they have high coverage in the controls and are unlikely to be genuine binding sites.

## Setting up the read extraction parameters

In the *[csaw](https://bioconductor.org/packages/3.10/csaw)* package, the `readParam` object determines which reads are extracted from the BAM files.
The intention is to set this up once and to re-use it in all relevant functions.
For this analysis, reads are ignored if they map to blacklist regions or do not map to the standard set of mouse nuclear chromosomes^[In this case, we are not interested in the mitochondrial genome, as these should not be bound by histones anyway.].


```r
library(csaw)
standard.chr <- paste0("chr", c(1:19, "X", "Y"))
param <- readParam(minq=20, discard=blacklist, restrict=standard.chr)
```

Reads are also ignored if they have a mapping quality (MAPQ) score below 20^[This is more stringent than usual, to account for the fact that the short reads ued here (32-36 bp) are more difficult to accurately align.].
This avoids spurious results due to weak or non-unique alignments that should be assigned low MAPQ scores by the aligner.
Note that the range of MAPQ scores will vary between aligners, so some inspection of the BAM files is necessary to choose an appropriate value.

# Computing the average fragment length

Strand bimodality is often observed in ChIP-seq experiments involving narrow binding events like H3K9ac marking.
This refers to the presence of distinct subpeaks on each strand and is quantified with cross-correlation plots [@kharchenko2008design].
A strong peak in the cross-correlations should be observed if immunoprecipitation was successful.
The delay distance at the peak corresponds to the distance between forward- and reverse-strand subpeaks.
This is identified from Figure \@ref(fig:ccfplot) and is used as the average fragment length for this analysis.


```r
x <- correlateReads(acdata$Path, param=reform(param, dedup=TRUE))
frag.len <- maximizeCcf(x)
frag.len
```

```
## [1] 154
```


```r
plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l")
abline(v=frag.len, col="red")
text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")
```

<div class="figure">
<img src="h3k9ac_files/figure-html/ccfplot-1.png" alt="Cross-correlation function (CCF) against delay distance for the H3K9ac data set. The delay with the maximum correlation is shown as the red line." width="100%" />
<p class="caption">(\#fig:ccfplot)Cross-correlation function (CCF) against delay distance for the H3K9ac data set. The delay with the maximum correlation is shown as the red line.</p>
</div>

Only unmarked reads (i.e., not potential PCR duplicates) are used to calculate the cross-correlations.
This reduces noise from variable PCR amplification and decreases the size of the "phantom" peak at the read length [@landt2012chipseq].
However, general removal of marked reads is risky as it caps the signal in high-coverage regions of the genome.
This can result in loss of power to detect DB, or introduction of spurious DB when the same cap is applied to libraries of different sizes.
Thus, the marking status of each read will be ignored in the rest of the analysis, i.e., no duplicates will be removed in downstream steps.

# Counting reads into windows

*[csaw](https://bioconductor.org/packages/3.10/csaw)* uses a sliding window strategy to quantify protein binding intensity across the genome.
Each read is directionally extended to the average fragment length (Figure \@ref(fig:directional)) to represent the DNA fragment from which that read was sequenced.
Any position within the inferred fragment is a potential contact site for the protein of interest.
To quantify binding in a genomic window, the number of these fragments overlapping the window is counted.
The window is then moved to its next position on the genome and counting is repeated^[Each read is usually counted into multiple windows, which will introduce correlations between adjacent windows but will not otherwise affect the analysis.].
This is done for all samples such that a count is obtained for each window in each sample. 

<div class="figure">
<img src="h3k9ac_files/figure-html/directional-1.png" alt="Directional extension of reads by the average fragment length `ext` in single-end ChIP-seq data. Each extended read represents an imputed fragment, and the number of fragments overlapping a window of a given `width` is counted." width="100%"  class="widefigure" />
<p class="caption">(\#fig:directional)Directional extension of reads by the average fragment length `ext` in single-end ChIP-seq data. Each extended read represents an imputed fragment, and the number of fragments overlapping a window of a given `width` is counted.</p>
</div>

The `windowCounts()` function produces a `RangedSummarizedExperiment` object containing a matrix of such counts.
Each row corresponds to a window; each column represents a BAM file corresponding to a single sample^[Counting can be parallelized across files using the `BPPARAM=` argument.];
and each entry of the matrix represents the number of fragments overlapping a particular window in a particular sample. 


```r
win.data <- windowCounts(acdata$Path, param=param, width=150, ext=frag.len)
win.data
```

```
## class: RangedSummarizedExperiment 
## dim: 1671254 4 
## metadata(6): spacing width ... param final.ext
## assays(1): counts
## rownames: NULL
## rowData names(0):
## colnames: NULL
## colData names(4): bam.files totals ext rlen
```

To analyze H3K9ac data, a window size of 150 bp is used here.
This corresponds roughly to the length of the DNA in a nucleosome [@humburg2011chipseqr], which is the smallest relevant unit for studying histone mark enrichment.
The spacing between windows is set to the default of 50 bp, i.e., the start positions for adjacent windows are 50 bp apart.
Smaller spacings can be used to improve spatial resolution, but will increase memory usage and runtime by increasing the number of windows required to cover the genome.
This is unnecessary as increased resolution confers little practical benefit for this data set -- counts for very closely spaced windows will be practically identical.
Finally, windows with very low counts (by default, less than a sum of 10 across all samples) are removed to reduce memory usage.
This represents a preliminary filter to remove uninteresting windows corresponding to likely background regions.

# Filtering windows by abundance

As previously mentioned, low-abundance windows contain no binding sites and need to be filtered out.
This improves power by removing irrelevant tests prior to the multiple testing correction;
avoids problems with discreteness in downstream statistical methods;
and reduces computational work for further analyses.
Here, filtering is performed using the average abundance of each window [@mccarthy2012differential], which is defined as the average log-count per million for that window.
This performs well as an independent filter statistic for NB-distributed count data [@lun2014denovo].

The filter threshold is defined based on the assumption that most regions in the genome are not marked by H3K9ac.
Reads are counted into large bins and the median coverage across those bins is used as an estimate of the background abundance^[Large bins are necessary to obtain a precise estimate of background coverage, which would otherwise be too low in individual windows.].
This estimate is then compared to the average abundances of the windows, after rescaling to account for differences in the window and bin sizes.
A window is only retained if its coverage is 3-fold higher than that of the background regions,
i.e., the abundance of the window is greater than the background abundance estimate by log~2~(3) or more.
This removes a large number of windows that are weakly or not marked and are likely to be irrelevant.


```r
bins <- windowCounts(acdata$Path, bin=TRUE, width=2000, param=param)
filter.stat <- filterWindows(win.data, bins, type="global")
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical  982167  689087
```

The effect of the fold-change threshold is examined visually in Figure \@ref(fig:bghistplot).
The chosen threshold is greater than the abundances of most bins in the genome -- presumably, those that contain background regions.
This suggests that the filter will remove most windows lying within background regions.


```r
hist(filter.stat$back.abundances, main="", breaks=50,
    xlab="Background abundance (log2-CPM)")
threshold <- filter.stat$abundances[1] - filter.stat$filter[1] + log2(min.fc)
abline(v=threshold, col="red")
```

<div class="figure">
<img src="h3k9ac_files/figure-html/bghistplot-1.png" alt="Histogram of average abundances across all 2 kbp genomic bins. The filter threshold is shown as the red line." width="100%" />
<p class="caption">(\#fig:bghistplot)Histogram of average abundances across all 2 kbp genomic bins. The filter threshold is shown as the red line.</p>
</div>

The filtering itself is done by simply subsetting the `RangedSummarizedExperiment` object.


```r
filtered.data <- win.data[keep,]
```

# Normalizing for sample-specific trended biases

Normalization is required to eliminate confounding sample-specific biases prior to any comparisons between samples.
Here, a trended bias is present between samples in Figure \@ref(fig:trendplot).
This refers to a systematic fold-difference in per-window coverage between samples that changes according to the average abundance of the window.


```r
win.ab <- scaledAverage(filtered.data)
adjc <- calculateCPM(filtered.data, use.offsets=FALSE)
logfc <- adjc[,4] - adjc[,1]
smoothScatter(win.ab, logfc, ylim=c(-6, 6), xlim=c(0, 5),
    xlab="Average abundance", ylab="Log-fold change")

lfit <- smooth.spline(logfc~win.ab, df=5)
o <- order(win.ab)
lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)
```

<div class="figure">
<img src="h3k9ac_files/figure-html/trendplot-1.png" alt="Abundance-dependent trend in the log-fold change between two H3K9ac samples (mature B over pro-B), across all windows retained after filtering. A smoothed spline fitted to the log-fold change against the average abundance is also shown in red." width="100%" />
<p class="caption">(\#fig:trendplot)Abundance-dependent trend in the log-fold change between two H3K9ac samples (mature B over pro-B), across all windows retained after filtering. A smoothed spline fitted to the log-fold change against the average abundance is also shown in red.</p>
</div>

Trended biases cannot be removed by scaling methods like TMM normalization [@robinson2010scaling], as the amount of scaling required varies with the abundance of the window.
Rather, non-linear normalization methods must be used.
*[csaw](https://bioconductor.org/packages/3.10/csaw)* implements a version of the fast loess method [@ballman2004fast] that has been modified to handle count data [@lun2015csaw].
This produces a matrix of offsets that can be used during model fitting.


```r
filtered.data <- normOffsets(filtered.data)
offsets <- assay(filtered.data, "offset")
head(offsets)
```

```
##           [,1]      [,2]       [,3]       [,4]
## [1,] 0.5372350 0.3457091 -0.4860457 -0.3968984
## [2,] 0.5082160 0.3238545 -0.4601253 -0.3719453
## [3,] 0.5021301 0.3192584 -0.4546073 -0.3667813
## [4,] 0.6274643 0.4119203 -0.5571155 -0.4822690
## [5,] 0.7049858 0.4638779 -0.6039740 -0.5648898
## [6,] 0.7260917 0.4778642 -0.6174239 -0.5865321
```

The effect of non-linear normalization is visualized with another mean-difference plot.
Once the offsets are applied to adjust the log-fold changes, the trend is eliminated from the plot (Figure \@ref(fig:normplot)).
The cloud of points is also centred at a log-fold change of zero.
This indicates that normalization was successful in removing the differences between samples. 


```r
norm.adjc <- calculateCPM(filtered.data, use.offsets=TRUE)
norm.fc <- norm.adjc[,4]-norm.adjc[,1]
smoothScatter(win.ab, norm.fc, ylim=c(-6, 6), xlim=c(0, 5),
    xlab="Average abundance", ylab="Log-fold change")

lfit <- smooth.spline(norm.fc~win.ab, df=5)
lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)
```

<div class="figure">
<img src="h3k9ac_files/figure-html/normplot-1.png" alt="Effect of non-linear normalization on the trended bias between two H3K9ac samples. Normalized log-fold changes are shown for all windows retained after filtering. A smoothed spline fitted to the log-fold change against the average abundance is also shown in red." width="100%" />
<p class="caption">(\#fig:normplot)Effect of non-linear normalization on the trended bias between two H3K9ac samples. Normalized log-fold changes are shown for all windows retained after filtering. A smoothed spline fitted to the log-fold change against the average abundance is also shown in red.</p>
</div>

The implicit assumption of non-linear methods is that most windows at each abundance are not DB.
Any systematic difference between samples is attributed to bias and is removed.
The assumption of a non-DB majority is reasonable for this data set, given that the cell types being compared are quite closely related.
However, it is not appropriate in cases where large-scale DB is expected, as removal of the difference would result in loss of genuine DB.
An alternative normalization strategy for these situations will be described later in the [CBP analysis](https://bioconductor.org/packages/3.10/chipseqDB/vignettes/cbp.html#normalization-for-composition-biases).

# Statistical modelling of biological variability

## Setting up the design matrix 

Counts are modelled using negative binomial generalized linear models (NB GLMs) in the *[edgeR](https://bioconductor.org/packages/3.10/edgeR)* package [@mccarthy2012differential; @robinson2010edger].
The NB distribution is useful as it can handle low, discrete counts for each window.
The NB dispersion parameter allows modelling of biological variability between replicate samples. 
GLMs can also accommodate complex experimental designs, though a simple design is sufficient for this study.


```r
celltype <- acdata$Description
celltype[grep("pro", celltype)] <- "proB"
celltype[grep("mature", celltype)] <- "matureB"

celltype <- factor(celltype)
design <- model.matrix(~0+celltype)
colnames(design) <- levels(celltype)
design
```

```
##   matureB proB
## 1       0    1
## 2       0    1
## 3       1    0
## 4       1    0
## attr(,"assign")
## [1] 1 1
## attr(,"contrasts")
## attr(,"contrasts")$celltype
## [1] "contr.treatment"
```

As a general rule, the experimental design should contain at least two replicates in each of the biological conditions.
This ensures that the results for each condition are replicable and are not the result of technical artifacts such as PCR duplicates.
Obviously, more replicates will provide more power to detect DB accurately and reliability, albeit at the cost of time and experimental resources.

## Estimating the NB dispersion

The `RangedSummarizedExperiment` object is coerced into a `DGEList` object (plus offsets) for use in *[edgeR](https://bioconductor.org/packages/3.10/edgeR)*.
Estimation of the NB dispersion is performed using the `estimateDisp` function.
Specifically, a NB dispersion trend is fitted to all windows against the average abundance.
This means that empirical mean-dispersion trends can be flexibly modelled.


```r
library(edgeR)
y <- asDGEList(filtered.data)
str(y)
```

```
## Formal class 'DGEList' [package "edgeR"] with 1 slot
##   ..@ .Data:List of 3
##   .. ..$ : int [1:689087, 1:4] 6 6 7 12 15 17 24 22 25 24 ...
##   .. .. ..- attr(*, "dimnames")=List of 2
##   .. .. .. ..$ : chr [1:689087] "1" "2" "3" "4" ...
##   .. .. .. ..$ : chr [1:4] "Sample1" "Sample2" "Sample3" "Sample4"
##   .. ..$ :'data.frame':	4 obs. of  3 variables:
##   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1
##   .. .. ..$ lib.size    : int [1:4] 8392971 7269175 3792141 4241789
##   .. .. ..$ norm.factors: num [1:4] 1 1 1 1
##   .. ..$ : num [1:689087, 1:4] 16.1 16 16 16.2 16.2 ...
```

```r
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.04096 0.05252 0.06165 0.06075 0.07209 0.07395
```

The NB dispersion trend is visualized in Figure \@ref(fig:bcvplot) as the biological coefficient of variation (BCV), i.e., the square root of the NB dispersion.
Note that only the trended dispersion will be used in the downstream steps -- the common and tagwise values are only shown for diagnostic purposes.
Specifically, the common BCV provides an overall measure of the variability in the data set, averaged across all windows.
Data sets with common BCVs ranging from 10 to 20% are considered to have low variability, i.e., counts are highly reproducible.
The tagwise BCVs should also be dispersed above and below the fitted trend, indicating that the fit was successful.


```r
plotBCV(y)
```

<div class="figure">
<img src="h3k9ac_files/figure-html/bcvplot-1.png" alt="Abundance-dependent trend in the BCV for each window, represented by the blue line. Common (red) and tagwise estimates (black) are also shown." width="100%" />
<p class="caption">(\#fig:bcvplot)Abundance-dependent trend in the BCV for each window, represented by the blue line. Common (red) and tagwise estimates (black) are also shown.</p>
</div>

For most sequencing count data, we expect to see a decreasing trend that plateaus with increasing average abundance.
This reflects the greater reliability of large counts, where the effects of stochasticity and technical artifacts (e.g., mapping errors, PCR duplicates) are averaged out.
We observe no clear trend in Figure \@ref(fig:bcvplot) as the windows have already been filtered to the plateau.
This is still a satisfactory result as it indicates that the retained windows have low variability and more power to detect DB.

## Estimating the QL dispersion

Additional modelling is provided with the QL methods in *[edgeR](https://bioconductor.org/packages/3.10/edgeR)* [@lund2012ql].
This introduces a QL dispersion parameter for each window, which captures variability in the NB dispersion around the fitted trend for each window.
Thus, the QL dispersion can model window-specific variability, whereas the NB dispersion trend is averaged across many windows.
However, with limited replicates, there is not enough information for each window to stably estimate the QL dispersion.
This is overcome by sharing information between windows with empirical Bayes (EB) shrinkage.
The instability of the QL dispersion estimates is reduced by squeezing the estimates towards an abundance-dependent trend (Figure \@ref(fig:qlplot)).


```r
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
```

<div class="figure">
<img src="h3k9ac_files/figure-html/qlplot-1.png" alt="Effect of EB shrinkage on the raw QL dispersion estimate for each window (black) towards the abundance-dependent trend (blue) to obtain squeezed estimates (red)." width="100%" />
<p class="caption">(\#fig:qlplot)Effect of EB shrinkage on the raw QL dispersion estimate for each window (black) towards the abundance-dependent trend (blue) to obtain squeezed estimates (red).</p>
</div>

The extent of shrinkage is determined by the prior degrees of freedom (d.f.).
Large prior d.f. indicates that the dispersions were similar across windows, such that strong shrinkage to the trend could be performed to increase stability and power.
Small prior d.f. indicates that the dispersions were more variable.
In such cases, less squeezing is performed as strong shrinkage would be inappropriate.


```r
summary(fit$df.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2238 15.4949 15.4949 15.2544 15.4949 15.4949
```

Also note the use of `robust=TRUE` in the `glmQLFit()` call, which reduces the sensitivity of the EB procedures to outlier variances.
This is particularly noticeable in Figure \@ref(fig:qlplot) with highly variable windows that (correctly) do not get squeezed towards the trend. 

## Examining the data with MDS plots

Multi-dimensional scaling (MDS) plots are used to examine the similarities between samples. 
The distance between a pair of samples on this plot represents the overall log-fold change between those samples. 
Ideally, replicates should cluster together while samples from different conditions should be separate.
While the mature B replicates are less tightly grouped, samples still separate by cell type in Figure \@ref(fig:mdsplot).
This suggests that our downstream analysis will be able to detect significant differences in enrichment between cell types.


```r
plotMDS(norm.adjc, labels=celltype,
    col=c("red", "blue")[as.integer(celltype)])
```

<div class="figure">
<img src="h3k9ac_files/figure-html/mdsplot-1.png" alt="MDS plot with two dimensions for all samples in the H3K9ac data set. Samples are labelled and coloured according to the cell type." width="100%" />
<p class="caption">(\#fig:mdsplot)MDS plot with two dimensions for all samples in the H3K9ac data set. Samples are labelled and coloured according to the cell type.</p>
</div>

# Testing for DB and controlling the FDR

## Testing for DB with QL F-tests

Each window is tested for significant differences between cell types using the QL F-test [@lund2012ql].
This is superior to the likelihood ratio test that is typically used for GLMs, as the QL F-test accounts for the uncertainty in dispersion estimation.
One $p$-value is produced for each window, representing the evidence against the null hypothesis (i.e., that no DB is present in the window).
For this analysis, the comparison is parametrized such that the reported log-fold change for each window represents that of the coverage in pro-B cells over their mature B counterparts.


```r
contrast <- makeContrasts(proB-matureB, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
head(res$table)
```

```
##       logFC    logCPM        F     PValue
## 1 1.3657940 0.3096885 2.305514 0.14677966
## 2 1.3564260 0.2624458 2.304629 0.14685257
## 3 2.0015383 0.2502879 3.940307 0.06304923
## 4 2.0779767 0.5033489 5.576372 0.03003408
## 5 0.8842093 0.8051419 1.651922 0.21545238
## 6 0.9678112 0.8948537 2.054733 0.16936408
```

## Controlling the FDR across regions

One might attempt to control the FDR by applying the Benjamini-Hochberg (BH) method to the window-level $p$-values [@benjamini1995controlling].
However, the features of interest are not windows, but the genomic regions that they represent.
Control of the FDR across windows does not guarantee control of the FDR across regions [@lun2014denovo].
The latter is arguably more relevant for the final interpretation of the results.

We instead control the region-level FDR by aggregating windows into regions and combining the $p$-values.
Here, adjacent windows less than 100 bp apart are aggregated into clusters.
Each cluster represents a genomic region.
Smaller values of `tol` allow distinct marking events to kept separate,
while larger values provide a broader perspective, e.g., by considering adjacent co-regulated sites as a single entity.
Chaining effects are mitigated by setting a maximum cluster width of 5 kbp.


```r
merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)
merged$region
```

```
## GRanges object with 41616 ranges and 0 metadata columns:
##           seqnames            ranges strand
##              <Rle>         <IRanges>  <Rle>
##       [1]     chr1   4775451-4775750      *
##       [2]     chr1   4785001-4786300      *
##       [3]     chr1   4807251-4807750      *
##       [4]     chr1   4808001-4808600      *
##       [5]     chr1   4857051-4858950      *
##       ...      ...               ...    ...
##   [41612]     chrY 73038001-73038400      *
##   [41613]     chrY 75445801-75446200      *
##   [41614]     chrY 88935951-88936350      *
##   [41615]     chrY 90554201-90554400      *
##   [41616]     chrY 90812801-90813100      *
##   -------
##   seqinfo: 21 sequences from an unspecified genome
```

A combined $p$-value is computed for each cluster using the method of @simes1986, 
based on the $p$-values of the constituent windows.
This represents the evidence against the global null hypothesis for each cluster, i.e., that no DB exists in any of its windows.
Rejection of this global null indicates that the cluster (and the region that it represents) contains DB.
Applying the BH method to the combined $p$-values allows the region-level FDR to be controlled.


```r
tabcom <- combineTests(merged$id, res$table)
head(tabcom)
```

```
## DataFrame with 6 rows and 6 columns
##    nWindows  logFC.up logFC.down              PValue                FDR
##   <integer> <integer>  <integer>           <numeric>          <numeric>
## 1         3         3          0   0.146852572988762   0.24642809183469
## 2        24         9          0  0.0882966655233521  0.168735548166406
## 3         8         1          3   0.526424531180179  0.648041273430584
## 4        10         1          2   0.729650939287535  0.829875744449031
## 5        36        14          6  0.0208882037608781 0.0605414172269767
## 6         3         0          3 0.00884432996167321  0.034553893018911
##     direction
##   <character>
## 1          up
## 2          up
## 3       mixed
## 4       mixed
## 5          up
## 6        down
```

Each row of the output table contains the statistics for a single cluster, 
including the combined *p*-value before and after the BH correction.
Additional fields include `nWindows`, the total number of windows in the cluster; 
`logFC.up`, the number of windows with a DB log-fold change above 0.5; 
and `log.FC.down`, the number fof windows with a log-fold change below -0.5.

## Examining the scope and direction of DB

We can easily calculate the total number of DB regions at a FDR of 5%.


```r
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)
```

```
##    Mode   FALSE    TRUE 
## logical   28515   13101
```

Determining the direction of DB is more complicated, as clusters may contain windows that are changing in opposite directions.
One approach is to use the direction of DB from the windows that contribute most to the combined $p$-value, as reported in the `direction` field for each cluster.
If significance is driven by windows changing in both directions, the direction for the cluster is defined as `"mixed"`.
Otherwise, the reported direction is the same as that of the windows, i.e., `"up"` or `"down"`.


```r
table(tabcom$direction[is.sig])
```

```
## 
##  down mixed    up 
##  8580   154  4367
```

Another approach is to use the log-fold change of the most significant window (identified with the `getBestTest()` function) as a proxy for the log-fold change of the cluster.


```r
tabbest <- getBestTest(merged$id, res$table)
head(tabbest)
```

```
## DataFrame with 6 rows and 6 columns
##        best              logFC            logCPM                F
##   <integer>          <numeric>         <numeric>        <numeric>
## 1         3   2.00153829256075 0.250287926080891  3.9403070089186
## 2        15   6.45488225628961 0.712521465636628 11.9826122454754
## 3        35    1.1783996686851 0.727376262089356 2.51421074291099
## 4        43 -0.908825402006814   1.0234078969562 2.74637354889958
## 5        60   6.57273805081489 0.809667875887879 14.9826446406268
## 6        82  -5.64002491507422 0.584285981417552 14.0521867903502
##               PValue                FDR
##            <numeric>          <numeric>
## 1  0.189147676521167  0.335560137526852
## 2 0.0882966655233521  0.190582668484939
## 3                  1                  1
## 4                  1                  1
## 5 0.0464345919412199  0.121535973473321
## 6 0.0176886599233464 0.0636112904186881
```

In the above table, each row contains the statistics for each cluster.
Of interest are the `best` and `logFC` fields.
The former is the index of the window that is the most significant in each cluster, 
while the latter is the log-fold change of that window.
This is used to obtain a summary of the direction of DB across all clusters.


```r
is.sig.pos <- (tabbest$logFC > 0)[is.sig]
summary(is.sig.pos)
```

```
##    Mode   FALSE    TRUE 
## logical    8664    4437
```

This approach is generally satisfactory, though it will not capture multiple changes in opposite directions^[Try `mixedClusters()` to formally detect clusters that contain significant changes in both directions.].
It also tends to overstate the magnitude of the log-fold change in each cluster.

# Saving results to file

One approach to saving results is to store all statistics in the metadata of a `GRanges` object.
This is useful as it keeps the statistics and coordinates together for each cluster, 
avoiding problems with synchronization in downstream steps.
We also store the midpoint and log-fold change of the most significant window in each cluster.
The updated `GRanges` object is then saved to file as a serialized R object with the `saveRDS` function.


```r
out.ranges <- merged$region
mcols(out.ranges) <- data.frame(tabcom,
    best.pos=mid(ranges(rowRanges(filtered.data[tabbest$best]))),
    best.logFC=tabbest$logFC)
saveRDS(file="h3k9ac_results.rds", out.ranges)
```

For input into other programs like genome browsers, results need to be saved in a more conventional format.
Here, coordinates of DB regions are saved in BED format via *[rtracklayer](https://bioconductor.org/packages/3.10/rtracklayer)*, using the log-transformed FDR as the score.


```r
simplified <- out.ranges[is.sig]
simplified$score <- -10*log10(simplified$FDR)
export(con="h3k9ac_results.bed", object=simplified)
```

# Interpreting the DB results

## Adding gene-centric annotation

### Using the `detailRanges` function

*[csaw](https://bioconductor.org/packages/3.10/csaw)* provides its own annotation function, `detailRanges()`.
This identifies all genic features overlapping each region and reports them in a compact string form.
Briefly, features are reported as `SYMBOL:STRAND:TYPE` where `SYMBOL` represents the gene symbol;
`STRAND` reports the strand of the gene; and `TYPE` reports the type(s) of overlapped feature, 
e.g., `E` for exons, `P` for promoters, `I` for introns^[Introns are only reported if an exon is not overlapped.].


```r
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
anno <- detailRanges(out.ranges, orgdb=org.Mm.eg.db,
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene)
head(anno$overlap)
```

```
## [1] "Mrpl15:-:E"            "Mrpl15:-:PE"           "Lypla1:+:P"           
## [4] "Lypla1:+:PE"           "Lypla1:+:I,Tcea1:+:PE" "Rgs20:-:I"
```

Annotated features that flank the region of interest are also reported.
The description for each feature is formatted as described above but the `TYPE` instead represents the distance (in base pairs) between the closest exon of the gene and the region.
By default, only flanking features within 5 kbp of each region are considered.


```r
head(anno$left)
```

```
## [1] "Mrpl15:-:935" "Mrpl15:-:896" ""             "Lypla1:+:19" 
## [5] ""             ""
```

```r
head(anno$right)
```

```
## [1] "Mrpl15:-:627" ""             "Lypla1:+:38"  ""            
## [5] ""             ""
```

The annotation for each region is stored in the metadata of the `GRanges` object.
The compact string form is useful for human interpretation, as it allows rapid examination of all genic features neighbouring each region.


```r
meta <- mcols(out.ranges)
mcols(out.ranges) <- data.frame(meta, anno)
```

### Using the *[ChIPpeakAnno](https://bioconductor.org/packages/3.10/ChIPpeakAnno)* package

As its name suggests, the *[ChIPpeakAnno](https://bioconductor.org/packages/3.10/ChIPpeakAnno)* package is designed to annotate peaks from ChIP-seq experiments [@zhu2010chippeakanno].
A `GRanges` object containing all regions of interest is supplied to the relevant function after removing all previous metadata fields to reduce clutter.
The gene closest to each region is then reported.
Gene coordinates are taken from the NCBI mouse 38 annotation, which is roughly equivalent to the annotation in the mm10 genome build.


```r
library(ChIPpeakAnno)
data(TSS.mouse.GRCm38)
minimal <- out.ranges
elementMetadata(minimal) <- NULL
anno.regions <- annotatePeakInBatch(minimal, AnnotationData=TSS.mouse.GRCm38)
colnames(elementMetadata(anno.regions))
```

```
## [1] "peak"                     "feature"                 
## [3] "start_position"           "end_position"            
## [5] "feature_strand"           "insideFeature"           
## [7] "distancetoFeature"        "shortestDistance"        
## [9] "fromOverlappingOrNearest"
```

Alternatively, identification of all overlapping features within, say, 5 kbp can be achieved by setting `maxgap=5000` and `output="overlapping"` in `annotatePeakInBatch`.
This will report each overlapping feature in a separate entry of the returned `GRanges` object, i.e., each input region may have multiple output values.
In contrast, `detailRanges()` will report all overlapping features for a region as a single string, i.e., each input region has one output value.
Which is preferable depends on the purpose of the annotation -- the `detailRanges()` output is more convenient for direct annotation of a DB list, while the `annotatePeakInBatch()` output contains more information and is more convenient for further manipulation.

## Reporting gene-based results

Another approach to annotation is to flip the problem around such that DB statistics are reported directly for features of interest like genes.
This is more convenient when the DB analysis needs to be integrated with, e.g., differential expression analyses of matched RNA-seq data.
In the code below, promoter coordinates and gene symbols are obtained from various annotation objects.


```r
prom <- suppressWarnings(promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,
    upstream=3000, downstream=1000, columns=c("tx_name", "gene_id")))
entrez.ids <- sapply(prom$gene_id, FUN=function(x) x[1]) # Using the first Entrez ID.
gene.name <- select(org.Mm.eg.db, keys=entrez.ids, keytype="ENTREZID", column="SYMBOL")
prom$gene_name <- gene.name$SYMBOL[match(entrez.ids, gene.name$ENTREZID)]
head(prom)
```

```
## GRanges object with 6 ranges and 3 metadata columns:
##                        seqnames          ranges strand |
##                           <Rle>       <IRanges>  <Rle> |
##   ENSMUST00000193812.1     chr1 3070253-3074252      + |
##   ENSMUST00000082908.1     chr1 3099016-3103015      + |
##   ENSMUST00000192857.1     chr1 3249757-3253756      + |
##   ENSMUST00000161581.1     chr1 3463587-3467586      + |
##   ENSMUST00000192183.1     chr1 3528795-3532794      + |
##   ENSMUST00000193244.1     chr1 3677155-3681154      + |
##                                     tx_name         gene_id   gene_name
##                                 <character> <CharacterList> <character>
##   ENSMUST00000193812.1 ENSMUST00000193812.1            <NA>        <NA>
##   ENSMUST00000082908.1 ENSMUST00000082908.1            <NA>        <NA>
##   ENSMUST00000192857.1 ENSMUST00000192857.1            <NA>        <NA>
##   ENSMUST00000161581.1 ENSMUST00000161581.1            <NA>        <NA>
##   ENSMUST00000192183.1 ENSMUST00000192183.1            <NA>        <NA>
##   ENSMUST00000193244.1 ENSMUST00000193244.1            <NA>        <NA>
##   -------
##   seqinfo: 66 sequences (1 circular) from mm10 genome
```

All windows overlapping each promoter are defined as a cluster.
DB statistics are computed for each cluster/promoter with the `combineOverlaps()` function, which behaves like the `combineTests()` function for `Hits` input objects.
This directly yields DB results for the annotated features.
Promoters with no overlapping windows are assigned `NA` values for the various fields and are filtered out below for demonstration purposes.


```r
olap <- findOverlaps(prom, rowRanges(filtered.data))
tabprom <- combineOverlaps(olap, res$table)
with.anno <- data.frame(ID=prom$tx_name, Gene=prom$gene_name, tabprom)
head(with.anno[!is.na(with.anno$PValue),])
```

```
##                       ID   Gene nWindows logFC.up logFC.down    PValue
## 14  ENSMUST00000134384.7 Lypla1       18        2          5 0.7004649
## 15 ENSMUST00000027036.10 Lypla1       18        2          5 0.7004649
## 16  ENSMUST00000150971.7 Lypla1       18        2          5 0.7004649
## 17  ENSMUST00000155020.1 Lypla1       18        2          5 0.7004649
## 18  ENSMUST00000119612.8 Lypla1       18        2          5 0.7004649
## 19  ENSMUST00000137887.7 Lypla1       18        2          5 0.7004649
##          FDR direction
## 14 0.7393987     mixed
## 15 0.7393987     mixed
## 16 0.7393987     mixed
## 17 0.7393987     mixed
## 18 0.7393987     mixed
## 19 0.7393987     mixed
```

Note that this strategy is distinct from counting reads across promoters.
Using promoter-level counts would not provide enough spatial resolution to detect sharp binding events that only occur in a subinterval of the promoter.
In particular, detection may be compromised by non-specific background or the presence of multiple opposing DB events in the same promoter.
Combining window-level statistics is preferable as resolution is maintained for optimal performance.

# Visualizing DB results

## Overview

Here, the *[Gviz](https://bioconductor.org/packages/3.10/Gviz)* package is used to visualize read coverage across the data set at regions of interest [@hahne2016visualizing].
Coverage in each BAM file will be represented by a single track.
Several additional tracks will also be included in each plot.
One is the genome axis track, to display the genomic coordinates across the plotted region.
The other is the annotation track containing gene models, with gene IDs replaced by symbols (where possible) for easier reading.


```r
library(Gviz)
gax <- GenomeAxisTrack(col="black", fontsize=15, size=2)
greg <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, showId=TRUE,
    geneSymbol=TRUE, name="", background.title="transparent")
symbols <- unlist(mapIds(org.Mm.eg.db, gene(greg), "SYMBOL",
    "ENTREZID", multiVals = "first"))
symbol(greg) <- symbols[gene(greg)]
```

We will also sort the DB regions by p-value for easier identification of regions of interest.


```r
o <- order(out.ranges$PValue)
sorted.ranges <- out.ranges[o]
sorted.ranges
```

```
## GRanges object with 41616 ranges and 11 metadata columns:
##           seqnames              ranges strand |  nWindows  logFC.up
##              <Rle>           <IRanges>  <Rle> | <integer> <integer>
##       [1]    chr17   34285101-34290050      * |        97         0
##       [2]     chr9 109050201-109053150      * |        57         0
##       [3]    chr17   34261151-34265850      * |        92        11
##       [4]    chr17   34306001-34308650      * |        51         0
##       [5]    chr18   60802751-60805750      * |        55         0
##       ...      ...                 ...    ... .       ...       ...
##   [41612]    chr18   23751901-23753200      * |        22         0
##   [41613]    chr12   83922051-83922650      * |        10         2
##   [41614]    chr15   99395101-99395650      * |         8         0
##   [41615]     chr3   67504201-67504500      * |         4         0
##   [41616]     chr4   43043401-43043700      * |         4         0
##           logFC.down               PValue                  FDR   direction
##            <integer>            <numeric>            <numeric> <character>
##       [1]         97 4.04797773668499e-11 1.22569557219491e-06        down
##       [2]         57 7.13783570863201e-11 1.22569557219491e-06        down
##       [3]         78 8.83575239471534e-11 1.22569557219491e-06        down
##       [4]         51 1.23282008545495e-10 1.28262601690733e-06        down
##       [5]         55 2.06286478303051e-10 1.54430072588262e-06        down
##       ...        ...                  ...                  ...         ...
##   [41612]          2    0.999832725328934    0.999908062838253       mixed
##   [41613]          0      0.9998854212888    0.999908062838253       mixed
##   [41614]          0    0.999908062838253    0.999908062838253       mixed
##   [41615]          0    0.999908062838253    0.999908062838253          up
##   [41616]          0    0.999908062838253    0.999908062838253       mixed
##            best.pos         best.logFC                  overlap
##           <integer>          <numeric>                 <factor>
##       [1]  34287575  -7.18686332236642               H2-Aa:-:PE
##       [2] 109051575  -6.19603369122881              Shisa5:+:PE
##       [3]  34262025  -7.70114852451639              H2-Ab1:+:PE
##       [4]  34306075  -5.80798257689994              H2-Eb1:+:PE
##       [5]  60804525  -5.98346376178937                Cd74:+:PE
##       ...       ...                ...                      ...
##   [41612]  23752525  -0.77704982053454 Gm15972:-:PE,Mapre2:+:PE
##   [41613]  83922125  0.880874875114293                 Numb:-:P
##   [41614]  99395425 -0.411300240851034               Tmbim6:+:I
##   [41615]  67504275  0.491618329274888              Rarres1:-:I
##   [41616]  43043575   0.17425353617959              Fam214b:-:I
##                     left                                         right
##                 <factor>                                      <factor>
##       [1]    H2-Aa:-:565                                              
##       [2]                Atrip-trex1:-:4783,Trex1:-:4788,Shisa5:+:1713
##       [3]                                                H2-Ab1:+:1252
##       [4]                                                 H2-Eb1:+:941
##       [5]                                                  Cd74:+:2158
##       ...            ...                                           ...
##   [41612]   Gm15972:-:78                                  Mapre2:+:525
##   [41613]     Numb:-:117                                              
##   [41614]  Tmbim6:+:1371                                 Tmbim6:+:4007
##   [41615]                                                             
##   [41616] Fam214b:-:3106                                Fam214b:-:1948
##   -------
##   seqinfo: 21 sequences from an unspecified genome
```

## Simple DB across a broad region

We start by visualizing one of the top-ranking DB regions.
This represents a simple DB event where the entire region changes in one direction (Figure \@ref(fig:simplebroadplot)).
Specifically, it represents an increase in H3K9ac marking at the *H2-Aa* locus in mature B cells.
This is consistent with the expected biology --
    H3K9ac is a mark of active gene expression [@karmodiya2012h3k9]
    and MHCII components are upregulated in mature B cells [@hoffman2002changes].


```r
cur.region <- sorted.ranges[1]
cur.region
```

```
## GRanges object with 1 range and 11 metadata columns:
##       seqnames            ranges strand |  nWindows  logFC.up logFC.down
##          <Rle>         <IRanges>  <Rle> | <integer> <integer>  <integer>
##   [1]    chr17 34285101-34290050      * |        97         0         97
##                     PValue                  FDR   direction  best.pos
##                  <numeric>            <numeric> <character> <integer>
##   [1] 4.04797773668499e-11 1.22569557219491e-06        down  34287575
##              best.logFC    overlap        left    right
##               <numeric>   <factor>    <factor> <factor>
##   [1] -7.18686332236642 H2-Aa:-:PE H2-Aa:-:565         
##   -------
##   seqinfo: 21 sequences from an unspecified genome
```



One track is plotted for each sample, in addition to the coordinate and annotation tracks.
Coverage is plotted in terms of sequencing depth-per-million at each base.
This corrects for differences in library sizes between tracks.


```r
collected <- list()
lib.sizes <- filtered.data$totals/1e6
for (i in seq_along(acdata$Path)) {
    reads <- extractReads(bam.file=acdata$Path[i], cur.region, param=param)
    cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
    collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,10),
        name=acdata$Description[i], col.axis="black", col.title="black",
        fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))
```

<div class="figure">
<img src="h3k9ac_files/figure-html/simplebroadplot-1.png" alt="Coverage tracks for a simple DB event between pro-B and mature B cells, across a broad region in the H3K9ac data set. Read coverage for each sample is shown as a per-million value at each base." width="768" />
<p class="caption">(\#fig:simplebroadplot)Coverage tracks for a simple DB event between pro-B and mature B cells, across a broad region in the H3K9ac data set. Read coverage for each sample is shown as a per-million value at each base.</p>
</div>

## Complex DB across a broad region

Complex DB refers to situations where multiple DB events are occurring within the same enriched region.
These are identified as those clusters that contain windows changing in both directions^[Technically, we should use `mixedClusters()` for rigorous identification of regions with significant changes in both directions. However, for simplicity, we'll just use a more _ad hoc_ approach here.].
Here, one of the top-ranking complex clusters is selected for visualization.


```r
complex <- sorted.ranges$logFC.up > 0 & sorted.ranges$logFC.down > 0
cur.region <- sorted.ranges[complex][2]
cur.region
```

```
## GRanges object with 1 range and 11 metadata columns:
##       seqnames              ranges strand |  nWindows  logFC.up logFC.down
##          <Rle>           <IRanges>  <Rle> | <integer> <integer>  <integer>
##   [1]     chr5 122987201-122991450      * |        83        18         43
##                     PValue                  FDR   direction  best.pos
##                  <numeric>            <numeric> <character> <integer>
##   [1] 1.30976135102916e-08 1.33826057750574e-05        down 122990925
##              best.logFC                       overlap         left
##               <numeric>                      <factor>     <factor>
##   [1] -5.48534588563145 A930024E05Rik:+:PE,Kdm2b:-:PE Kdm2b:-:2230
##                      right
##                   <factor>
##   [1] A930024E05Rik:+:2913
##   -------
##   seqinfo: 21 sequences from an unspecified genome
```



This region contains a bidirectional promoter where different genes are marked in the different cell types (Figure \@ref(fig:complexplot)).
Upon differentiation to mature B cells, loss of marking in one part of the region is balanced by a gain in marking in another part of the region.
This represents a complex DB event that would not be detected if reads were counted across the entire region.


```r
collected <- list()
for (i in seq_along(acdata$Path)) {
    reads <- extractReads(bam.file=acdata$Path[i], cur.region, param=param)
    cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
    collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,3),
        name=acdata$Description[i], col.axis="black", col.title="black",
        fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))
```

<div class="figure">
<img src="h3k9ac_files/figure-html/complexplot-1.png" alt="Coverage tracks for a complex DB event in the H3K9ac data set, shown as per-million values." width="768" />
<p class="caption">(\#fig:complexplot)Coverage tracks for a complex DB event in the H3K9ac data set, shown as per-million values.</p>
</div>

## Simple DB across a small region

Both of the examples above involve differential marking within broad regions spanning several kilobases.
This is consistent with changes in the marking profile across a large number of nucleosomes.
However, H3K9ac marking can also be concentrated into small regions, involving only a few nucleosomes.
*[csaw](https://bioconductor.org/packages/3.10/csaw)* is equally capable of detecting sharp DB within these small regions.
This is demonstrated by examining those clusters that contain a smaller number of windows.


```r
sharp <- sorted.ranges$nWindows < 20
cur.region <- sorted.ranges[sharp][1]
cur.region
```

```
## GRanges object with 1 range and 11 metadata columns:
##       seqnames            ranges strand |  nWindows  logFC.up logFC.down
##          <Rle>         <IRanges>  <Rle> | <integer> <integer>  <integer>
##   [1]    chr16 36665551-36666200      * |        11         0         11
##                     PValue                  FDR   direction  best.pos
##                  <numeric>            <numeric> <character> <integer>
##   [1] 1.29839663897595e-08 1.33826057750574e-05        down  36665925
##              best.logFC   overlap        left    right
##               <numeric>  <factor>    <factor> <factor>
##   [1] -4.93341819257933 Cd86:-:PE Cd86:-:3937         
##   -------
##   seqinfo: 21 sequences from an unspecified genome
```



Marking is increased for mature B cells within a 500 bp region (Figure \@ref(fig:simplesharpplot)), which is sharper than the changes in the previous two examples.
This also coincides with the promoter of the *Cd86* gene.
Again, this makes biological sense as CD86 is involved in regulating immunoglobulin production in activated B-cells [@podojil2003selective].


```r
collected <- list()
for (i in seq_along(acdata$Path)) {
    reads <- extractReads(bam.file=acdata$Path[i], cur.region, param=param)
    cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
    collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,3),
        name=acdata$Description[i], col.axis="black", col.title="black",
        fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))
```

<div class="figure">
<img src="h3k9ac_files/figure-html/simplesharpplot-1.png" alt="Coverage tracks for a sharp and simple DB event in the H3K9ac data set, shown as per-million values." width="768" />
<p class="caption">(\#fig:simplesharpplot)Coverage tracks for a sharp and simple DB event in the H3K9ac data set, shown as per-million values.</p>
</div>

Note that the window size will determine whether sharp or broad events are preferentially detected.
Larger windows provide more power to detect broad events (as the counts are higher), while smaller windows provide more resolution to detect sharp events.
Optimal detection of all features can be obtained by performing analyses with multiple window sizes and consolidating the results^[See `?consolidateWindows` and `?consolidateTests` for further information.], though -- for brevity -- this will not be described here.
In general, smaller window sizes are preferred as strong DB events with sufficient coverage will always be detected.
For larger windows, detection may be confounded by other events within the window that distort the log-fold change in the counts between conditions.

# Session information


```r
sessionInfo()
```

```
## R version 3.6.0 Patched (2019-05-02 r76458)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.2 LTS
## 
## Matrix products: default
## BLAS:   /home/luna/Software/R/R-3-6-branch-dev/lib/libRblas.so
## LAPACK: /home/luna/Software/R/R-3-6-branch-dev/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] Gviz_1.29.0                             
##  [2] ChIPpeakAnno_3.19.2                     
##  [3] VennDiagram_1.6.20                      
##  [4] futile.logger_1.4.3                     
##  [5] TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.7
##  [6] GenomicFeatures_1.37.1                  
##  [7] org.Mm.eg.db_3.8.2                      
##  [8] AnnotationDbi_1.47.0                    
##  [9] edgeR_3.27.3                            
## [10] limma_3.41.2                            
## [11] csaw_1.19.1                             
## [12] SummarizedExperiment_1.15.1             
## [13] DelayedArray_0.11.0                     
## [14] BiocParallel_1.19.0                     
## [15] matrixStats_0.54.0                      
## [16] Biobase_2.45.0                          
## [17] rtracklayer_1.45.1                      
## [18] BiocFileCache_1.9.0                     
## [19] dbplyr_1.4.0                            
## [20] Rsamtools_2.1.2                         
## [21] Biostrings_2.53.0                       
## [22] XVector_0.25.0                          
## [23] GenomicRanges_1.37.7                    
## [24] GenomeInfoDb_1.21.1                     
## [25] IRanges_2.19.5                          
## [26] S4Vectors_0.23.3                        
## [27] BiocGenerics_0.31.2                     
## [28] chipseqDBData_1.1.0                     
## [29] knitr_1.23                              
## [30] BiocStyle_2.13.0                        
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.4-1              seqinr_3.4-5                 
##   [3] htmlTable_1.13.1              biovizBase_1.33.0            
##   [5] base64enc_0.1-3               dichromat_2.0-0              
##   [7] rstudioapi_0.10               bit64_0.9-7                  
##   [9] interactiveDisplayBase_1.23.0 codetools_0.2-16             
##  [11] splines_3.6.0                 ade4_1.7-13                  
##  [13] Formula_1.2-3                 cluster_2.0.9                
##  [15] GO.db_3.8.2                   graph_1.63.0                 
##  [17] shiny_1.3.2                   BiocManager_1.30.4           
##  [19] compiler_3.6.0                httr_1.4.0                   
##  [21] backports_1.1.4               assertthat_0.2.1             
##  [23] Matrix_1.2-17                 lazyeval_0.2.2               
##  [25] later_0.8.0                   formatR_1.6                  
##  [27] acepack_1.4.1                 htmltools_0.3.6              
##  [29] prettyunits_1.0.2             tools_3.6.0                  
##  [31] gtable_0.3.0                  glue_1.3.1                   
##  [33] GenomeInfoDbData_1.2.1        dplyr_0.8.1                  
##  [35] rappdirs_0.3.1                Rcpp_1.0.1                   
##  [37] multtest_2.41.0               ExperimentHub_1.11.1         
##  [39] xfun_0.7                      stringr_1.4.0                
##  [41] mime_0.6                      ensembldb_2.9.1              
##  [43] statmod_1.4.30                XML_3.98-1.19                
##  [45] idr_1.2                       AnnotationHub_2.17.3         
##  [47] zlibbioc_1.31.0               MASS_7.3-51.4                
##  [49] scales_1.0.0                  BSgenome_1.53.0              
##  [51] VariantAnnotation_1.31.3      hms_0.4.2                    
##  [53] promises_1.0.1                ProtGenerics_1.17.2          
##  [55] RBGL_1.61.0                   AnnotationFilter_1.9.0       
##  [57] lambda.r_1.2.3                RColorBrewer_1.1-2           
##  [59] yaml_2.2.0                    curl_3.3                     
##  [61] gridExtra_2.3                 memoise_1.1.0                
##  [63] ggplot2_3.1.1                 rpart_4.1-15                 
##  [65] biomaRt_2.41.0                latticeExtra_0.6-28          
##  [67] stringi_1.4.3                 RSQLite_2.1.1                
##  [69] highr_0.8                     checkmate_1.9.3              
##  [71] rlang_0.3.4                   pkgconfig_2.0.2              
##  [73] bitops_1.0-6                  evaluate_0.13                
##  [75] lattice_0.20-38               purrr_0.3.2                  
##  [77] htmlwidgets_1.3               GenomicAlignments_1.21.2     
##  [79] bit_1.1-14                    tidyselect_0.2.5             
##  [81] plyr_1.8.4                    magrittr_1.5                 
##  [83] bookdown_0.10                 R6_2.4.0                     
##  [85] Hmisc_4.2-0                   DBI_1.0.0                    
##  [87] foreign_0.8-71                pillar_1.4.0                 
##  [89] nnet_7.3-12                   survival_2.44-1.1            
##  [91] RCurl_1.95-4.12               tibble_2.1.1                 
##  [93] crayon_1.3.4                  futile.options_1.0.1         
##  [95] KernSmooth_2.23-15            rmarkdown_1.13               
##  [97] progress_1.2.2                locfit_1.5-9.1               
##  [99] data.table_1.12.2             blob_1.1.1                   
## [101] digest_0.6.19                 xtable_1.8-4                 
## [103] httpuv_1.5.1                  regioneR_1.17.2              
## [105] munsell_0.5.0
```

# References
