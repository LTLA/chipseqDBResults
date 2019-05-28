---
title: "Detecting differential binding of CBP in mouse fibroblasts"
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
  %\VignetteIndexEntry{3. Differential binding of CBP in fibroblasts}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
output:
  BiocStyle::html_document:
    toc_float: yes
    titlecaps: false
bibliography: ref.bib
---



# Overview

Here, we perform a window-based DB analysis to identify differentially bound (DB) regions for CREB-binding protein (CBP).
This provides an example of how to use *[csaw](https://bioconductor.org/packages/3.10/csaw)* for transcription factor (TF) data, 
complementing our [previous analysis](https://bioconductor.org/packages/3.10/chipseqDB/vignettes/h3k9ac.html) on a histone mark.
We use CBP ChIP-seq data from a study comparing wild-type (WT) and CBP knock-out (KO) animals [@kasper2014genomewide], with two biological replicates for each genotype.
BAM files and indices are downloaded using *[chipseqDBData](https://bioconductor.org/packages/3.10/chipseqDBData)* and cached for later use.


```r
library(chipseqDBData)
cbpdata <- CBPData()
```

```r
cbpdata
```

```
## DataFrame with 4 rows and 3 columns
##          Name       Description
##   <character>       <character>
## 1  SRR1145787 CBP wild-type (1)
## 2  SRR1145788 CBP wild-type (2)
## 3  SRR1145789 CBP knock-out (1)
## 4  SRR1145790 CBP knock-out (2)
##                                             Path
##                                      <character>
## 1 /tmp/Rtmpcjbzy4/file8882c32cc39/SRR1145787.bam
## 2 /tmp/Rtmpcjbzy4/file8882c32cc39/SRR1145788.bam
## 3 /tmp/Rtmpcjbzy4/file8882c32cc39/SRR1145789.bam
## 4 /tmp/Rtmpcjbzy4/file8882c32cc39/SRR1145790.bam
```

# Pre-processing checks

We check some mapping statistics for the CBP dataset with *[Rsamtools](https://bioconductor.org/packages/3.10/Rsamtools)*, as previously described.


```r
library(Rsamtools)
diagnostics <- list()
for (bam in cbpdata$Path) {
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
##                   Total   Mapped  Marked Prop.mapped Prop.marked
## SRR1145787.bam 28525952 24289396 2022868    85.14842    8.328194
## SRR1145788.bam 25514465 21604007 1939224    84.67356    8.976224
## SRR1145789.bam 34476967 29195883 2412650    84.68228    8.263665
## SRR1145790.bam 32624587 27348488 2617879    83.82784    9.572299
```

We construct a `readParam` object to standardize the parameter settings in this analysis.
The ENCODE blacklist is again used^[Assuming you ran the previous workflow, this will be retrieved from cache rather than being downloaded again.] to remove reads in problematic regions [@encode2012encode].


```r
library(BiocFileCache)
bfc <- BiocFileCache("local", ask=FALSE)
black.path <- bfcrpath(bfc, file.path("https://www.encodeproject.org",
    "files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"))

library(rtracklayer)
blacklist <- import(black.path)
```

We set the minimum mapping quality score to 10 to remove poorly or non-uniquely aligned reads.


```r
library(csaw)
param <- readParam(minq=10, discard=blacklist)
param
```

```
##     Extracting reads in single-end mode
##     Duplicate removal is turned off 
##     Minimum allowed mapping score is 10 
##     Reads are extracted from both strands
##     No restrictions are placed on read extraction
##     Reads in 164 regions will be discarded
```

# Computing the average fragment length

The average fragment length is estimated by maximizing the cross-correlation function (Figure \@ref(fig:ccfplot)), as previously described.
Generally, cross-correlations for TF datasets are sharper than for histone marks as the TFs typically contact a smaller genomic interval.
This results in more pronounced strand bimodality in the binding profile.


```r
x <- correlateReads(cbpdata$Path, param=reform(param, dedup=TRUE))
frag.len <- maximizeCcf(x)
frag.len
```

```
## [1] 161
```


```r
plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l")
abline(v=frag.len, col="red")
text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")
```

<div class="figure">
<img src="cbp_files/figure-html/ccfplot-1.png" alt="Cross-correlation function (CCF) against delay distance for the CBP dataset. The delay with the maximum correlation is shown as the red line." width="100%" />
<p class="caption">(\#fig:ccfplot)Cross-correlation function (CCF) against delay distance for the CBP dataset. The delay with the maximum correlation is shown as the red line.</p>
</div>

# Counting reads into windows

Reads are then counted into sliding windows using *[csaw](https://bioconductor.org/packages/3.10/csaw)* [@lun2015csaw].
For TF data analyses, smaller windows are necessary to capture sharp binding sites.
A large window size will be suboptimal as the count for a particular site will be "contaminated" by non-specific background in the neighbouring regions.
In this case, a window size of 10 bp is used.


```r
win.data <- windowCounts(cbpdata$Path, param=param, width=10, ext=frag.len)
win.data
```

```
## class: RangedSummarizedExperiment 
## dim: 9952827 4 
## metadata(6): spacing width ... param final.ext
## assays(1): counts
## rownames: NULL
## rowData names(0):
## colnames: NULL
## colData names(4): bam.files totals ext rlen
```

The default spacing of 50 bp is also used here.
This may seem inappropriate given that the windows are only 10 bp.
However, reads lying in the interval between adjacent windows will still be counted into several windows.
This is because reads are extended to the value of `frag.len`, which is substantially larger than the 50 bp spacing^[Smaller spacings can be used but will provide little benefit given that each extended read already overlaps multiple windows.].

# Normalization for composition biases

Composition biases are introduced when the amount of DB in each condition is unbalanced [@robinson2010scaling; @lun2014denovo].
More binding in one condition means that more reads are sequenced at the binding sites, leaving fewer reads for the rest of the genome.
This suppresses the genomic coverage at non-DB sites, resulting in spurious differences between samples. 
We expect unbalanced DB in this dataset as CBP function should be compromised in the KO cells,
such that most - if not all - of the DB sites should exhibit increased CBP binding in the WT condition.

To remove this bias, we assign reads to large genomic bins and assume that most bins represent non-DB background regions.
Any systematic differences in the coverage of those bins is attributed to composition bias and is normalized out.
Specifically, the trimmed mean of M-values (TMM) method [@robinson2010scaling] is applied to compute normalization factors from the bin counts.
These factors are stored in `win.data`^[See the `se.out=` argument.] so that they will be applied during the DB analysis with the window counts.


```r
bins <- windowCounts(cbpdata$Path, bin=TRUE, width=10000, param=param)
win.data <- normFactors(bins, se.out=win.data)
(normfacs <- win.data$norm.factors)
```

```
## [1] 1.0125617 0.9083253 1.0443668 1.0410799
```

We visualize the effect of normalization with mean-difference plots between pairs of samples (Figure \@ref(fig:compoplot)).
The dense cloud in each plot represents the majority of bins in the genome.
These are assumed to mostly contain background regions.
A non-zero log-fold change for these bins indicates that composition bias is present between samples. 
The red line represents the log-ratio of normalization factors and passes through the centre of the cloud in each plot,
indicating that the bias has been successfully identified and removed.


```r
bin.ab <- scaledAverage(bins)
adjc <- calculateCPM(bins, use.norm.factors=FALSE)

par(cex.lab=1.5, mfrow=c(1,3))
smoothScatter(bin.ab, adjc[,1]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (1 vs 4)")
abline(h=log2(normfacs[1]/normfacs[4]), col="red")

smoothScatter(bin.ab, adjc[,2]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (2 vs 4)")
abline(h=log2(normfacs[2]/normfacs[4]), col="red")

smoothScatter(bin.ab, adjc[,3]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (3 vs 4)")
abline(h=log2(normfacs[3]/normfacs[4]), col="red")
```

<div class="figure">
<img src="cbp_files/figure-html/compoplot-1.png" alt="Mean-difference plots for the bin counts, comparing sample 4 to all other samples. The red line represents the log-ratio of the normalization factors between samples." width="1152" />
<p class="caption">(\#fig:compoplot)Mean-difference plots for the bin counts, comparing sample 4 to all other samples. The red line represents the log-ratio of the normalization factors between samples.</p>
</div>

Note that this normalization strategy is quite different from that in the H3K9ac analysis.
Here, systematic DB in one direction is expected between conditions, given that CBP function is lost in the KO genotype.
This means that the assumption of a non-DB majority (required for non-linear normalization of the [H3K9ac data](https://bioconductor.org/packages/3.10/chipseqDB/vignettes/h3k9ac.html#normalizing-for-sample-specific-trended-biases)) is not valid.
No such assumption is made by the binned-TMM approach described above, which makes it more appropriate for use in the CBP analysis.

# Filtering of low-abundance windows

Removal of low-abundance windows is performed as [previously described](https://bioconductor.org/packages/3.10/chipseqDB/vignettes/h3k9ac.html#filtering-windows-by-abundance), by computing the coverage in each window relative to a global estimate of background enrichment.
The majority of windows in background regions are filtered out upon applying a modest fold-change threshold.
This leaves a small set of relevant windows for further analysis.


```r
filter.stat <- filterWindows(win.data, bins, type="global")
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical 9653836  298991
```

```r
filtered.data <- win.data[keep,]
```

Note that the 10 kbp bins are used here for filtering, while smaller 2 kbp bins were used in the corresponding step for the H3K9ac analysis.
This is purely for convenience -- the 10 kbp counts for this dataset were previously loaded for normalization, and can be re-used during filtering to save time.
Changes in bin size will have little impact on the results, so long as the bins (and their counts) are large enough for precise estimation of the background abundance.
While smaller bins provide greater spatial resolution, this is irrelevant for quantifying coverage in large background regions that span most of the genome.

# Statistical modelling of biological variability

Counts for each window are modelled using *[edgeR](https://bioconductor.org/packages/3.10/edgeR)* [@mccarthy2012differential; @robinson2010edger].
We convert our `RangedSummarizedExperiment` object into a `DGEList`.


```r
library(edgeR)
y <- asDGEList(filtered.data)
summary(y)
```

```
##         Length  Class      Mode   
## counts  1195964 -none-     numeric
## samples       3 data.frame list
```

We then construct a design matrix for our experimental design.
Again, we have a simple one-way layout with two groups of two replicates.


```r
genotype <- cbpdata$Description
genotype[grep("wild-type", genotype)] <- "wt"
genotype[grep("knock-out", genotype)] <- "ko"

genotype <- factor(genotype)
design <- model.matrix(~0+genotype)
colnames(design) <- levels(genotype)
design
```

```
##   ko wt
## 1  0  1
## 2  0  1
## 3  1  0
## 4  1  0
## attr(,"assign")
## [1] 1 1
## attr(,"contrasts")
## attr(,"contrasts")$genotype
## [1] "contr.treatment"
```

We estimate the negative binomial (NB) and quasi-likelihood (QL) dispersions for each window [@lund2012ql].
The estimated NB dispersions (Figure \@ref(fig:bcvplot)) are substantially larger than those observed in the H3K9ac dataset.
They also exhibit an unusual increasing trend with respect to abundance.


```r
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1196  0.1645  0.1843  0.1889  0.2159  0.2549
```

```r
plotBCV(y)
```

<div class="figure">
<img src="cbp_files/figure-html/bcvplot-1.png" alt="Abundance-dependent trend in the biological coefficient of variation (i.e., the root-NB dispersion) for each window, represented by the blue line. Common (red) and tagwise estimates (black) are also shown." width="100%" />
<p class="caption">(\#fig:bcvplot)Abundance-dependent trend in the biological coefficient of variation (i.e., the root-NB dispersion) for each window, represented by the blue line. Common (red) and tagwise estimates (black) are also shown.</p>
</div>

The estimated prior d.f. is also infinite, meaning that all the QL dispersions are equal to the trend (Figure \@ref(fig:qlplot)).


```r
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   18260     Inf     Inf     Inf     Inf     Inf
```

```r
plotQLDisp(fit)
```

<div class="figure">
<img src="cbp_files/figure-html/qlplot-1.png" alt="Effect of EB shrinkage on the raw QL dispersion estimate for each window (black) towards the abundance-dependent trend (blue) to obtain squeezed estimates (red). Quarter-root estimates are shown for greater dynamic range." width="100%" />
<p class="caption">(\#fig:qlplot)Effect of EB shrinkage on the raw QL dispersion estimate for each window (black) towards the abundance-dependent trend (blue) to obtain squeezed estimates (red). Quarter-root estimates are shown for greater dynamic range.</p>
</div>

These results are consistent with the presence of a systematic difference in CBP enrichment between the WT replicates.
An increasing trend in Figure \@ref(fig:bcvplot) is typical after normalization for composition biases,
where replicates exhibit some differences in efficiency that manifest as increased dispersions at high abundance.
The dispersions for all windows are inflated to a similarly large value by this difference, 
manifesting as low variability in the dispersions across windows.
This effect is illustrated in Figure \@ref(fig:mdsplot) where the WT samples are clearly separated in both dimensions.


```r
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("red", "blue")[as.integer(genotype)])
```

<div class="figure">
<img src="cbp_files/figure-html/mdsplot-1.png" alt="MDS plot with two dimensions for all samples in the CBP dataset. Samples are labelled and coloured according to the genotype. A larger top set of windows was used to improve the visualization of the genome-wide differences between the WT samples." width="100%" />
<p class="caption">(\#fig:mdsplot)MDS plot with two dimensions for all samples in the CBP dataset. Samples are labelled and coloured according to the genotype. A larger top set of windows was used to improve the visualization of the genome-wide differences between the WT samples.</p>
</div>

The presence of a large batch effect between replicates is not ideal.
Nonetheless, we can still proceed with the DB analysis - albeit with some loss of power due to the inflated NB dispersions - 
given that there are strong differences between genotypes in Figure \@ref(fig:mdsplot),

# Testing for DB

We test for a significant difference in binding between genotypes in each window using the QL F-test.


```r
contrast <- makeContrasts(wt-ko, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
```

Windows are clustered into regions and the region-level FDR is controlled using Simes' method [@simes1986; @lun2014denovo].


```r
merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)
tabcom <- combineTests(merged$id, res$table)
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)
```

```
##    Mode   FALSE    TRUE 
## logical   59583    1832
```

All significant regions have increased CBP binding in the WT genotype.
This is expected given that protein function should be lost in the KO genotype.


```r
table(tabcom$direction[is.sig])
```

```
## 
##   up 
## 1832
```

```r
# Direction according the best window in each cluster.
tabbest <- getBestTest(merged$id, res$table)
is.sig.pos <- (tabbest$logFC > 0)[is.sig]
summary(is.sig.pos)
```

```
##    Mode    TRUE 
## logical    1832
```

These results are saved to file in the form of a serialized R object for later inspection.


```r
out.ranges <- merged$region
mcols(out.ranges) <- data.frame(tabcom,
    best.pos=mid(ranges(rowRanges(filtered.data[tabbest$best]))),
    best.logFC=tabbest$logFC)
saveRDS(file="cbp_results.rds", out.ranges)
```

# Annotation and visualization

Annotation for each region is added using the `detailRanges()` function,
as [previously described](https://bioconductor.org/packages/3.10/chipseqDB/vignettes/h3k9ac.html#using-the-detailranges-function).


```r
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
anno <- detailRanges(out.ranges, orgdb=org.Mm.eg.db,
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene)
mcols(out.ranges) <- cbind(mcols(out.ranges), anno)
```

We will visualize one of the top-ranked DB regions here.
This corresponds to a simple DB event as all windows are changing in the same direction, i.e., up in the WT.
The binding region is also quite small relative to some of the H3K9ac examples, 
consistent with sharp TF binding to a specific recognition site.


```r
o <- order(out.ranges$PValue)    
cur.region <- out.ranges[o[2]]
cur.region
```

```
## GRanges object with 1 range and 11 metadata columns:
##       seqnames            ranges strand |  nWindows  logFC.up logFC.down
##          <Rle>         <IRanges>  <Rle> | <integer> <integer>  <integer>
##   [1]    chr16 70313851-70314860      * |        21        21          0
##                     PValue                  FDR   direction  best.pos
##                  <numeric>            <numeric> <character> <integer>
##   [1] 1.32085091386953e-13 2.80252861892365e-09          up  70314505
##             best.logFC     overlap        left       right
##              <numeric> <character> <character> <character>
##   [1] 4.39455538103856   Gbe1:+:PE                        
##   -------
##   seqinfo: 66 sequences from an unspecified genome
```



We use *[Gviz](https://bioconductor.org/packages/3.10/Gviz)* [@hahne2016visualizing] to plot the results.
As in the [H3K9ac analysis](https://bioconductor.org/packages/3.10/chipseqDB/vignettes/h3k9ac.html#visualizing-db-results), 
we set up some tracks to display genome coordinates and gene annotation.


```r
library(Gviz)
gax <- GenomeAxisTrack(col="black", fontsize=15, size=2)
greg <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, showId=TRUE,
    geneSymbol=TRUE, name="", background.title="transparent")
symbols <- unlist(mapIds(org.Mm.eg.db, gene(greg), "SYMBOL",
    "ENTREZID", multiVals = "first"))
symbol(greg) <- symbols[gene(greg)]
```

We visualize two tracks for each sample -- one for the forward-strand coverage, another for the reverse-strand coverage.
This allows visualization of the strand bimodality that is characteristic of genuine TF binding sites.
In Figure \@ref(fig:tfplot), two adjacent sites are present at the *Gbe1* promoter, both of which exhibit increased binding in the WT genotype.
Coverage is also substantially different between the WT replicates, consistent with the presence of a batch effect.


```r
collected <- list()
lib.sizes <- filtered.data$totals/1e6

for (i in seq_along(cbpdata$Path)) {
    reads <- extractReads(bam.file=cbpdata$Path[i], cur.region, param=param)
    pcov <- as(coverage(reads[strand(reads)=="+"])/lib.sizes[i], "GRanges")
    ncov <- as(coverage(reads[strand(reads)=="-"])/-lib.sizes[i], "GRanges")
    ptrack <- DataTrack(pcov, type="histogram", lwd=0, ylim=c(-5, 5),
        name=cbpdata$Description[i], col.axis="black", col.title="black",
        fill="blue", col.histogram=NA)
    ntrack <- DataTrack(ncov, type="histogram", lwd=0, ylim=c(-5, 5),
        fill="red", col.histogram=NA)
    collected[[i]] <- OverlayTrack(trackList=list(ptrack, ntrack))
}

plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))
```

<div class="figure">
<img src="cbp_files/figure-html/tfplot-1.png" alt="Coverage tracks for TF binding sites that are differentially bound in the WT (top two tracks) against the KO (last two tracks). Blue and red tracks represent forward- and reverse-strand coverage, respectively, on a per-million scale (capped at 5 in SRR1145788, for visibility)." width="768" />
<p class="caption">(\#fig:tfplot)Coverage tracks for TF binding sites that are differentially bound in the WT (top two tracks) against the KO (last two tracks). Blue and red tracks represent forward- and reverse-strand coverage, respectively, on a per-million scale (capped at 5 in SRR1145788, for visibility).</p>
</div>

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
##  [2] org.Mm.eg.db_3.8.2                      
##  [3] TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.7
##  [4] GenomicFeatures_1.37.1                  
##  [5] AnnotationDbi_1.47.0                    
##  [6] edgeR_3.27.3                            
##  [7] limma_3.41.2                            
##  [8] csaw_1.19.1                             
##  [9] SummarizedExperiment_1.15.1             
## [10] DelayedArray_0.11.0                     
## [11] BiocParallel_1.19.0                     
## [12] matrixStats_0.54.0                      
## [13] Biobase_2.45.0                          
## [14] rtracklayer_1.45.1                      
## [15] BiocFileCache_1.9.0                     
## [16] dbplyr_1.4.0                            
## [17] Rsamtools_2.1.2                         
## [18] Biostrings_2.53.0                       
## [19] XVector_0.25.0                          
## [20] GenomicRanges_1.37.7                    
## [21] GenomeInfoDb_1.21.1                     
## [22] IRanges_2.19.5                          
## [23] S4Vectors_0.23.3                        
## [24] BiocGenerics_0.31.2                     
## [25] chipseqDBData_1.1.0                     
## [26] knitr_1.23                              
## [27] BiocStyle_2.13.0                        
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.4-1              biovizBase_1.33.0            
##  [3] htmlTable_1.13.1              base64enc_0.1-3              
##  [5] dichromat_2.0-0               rstudioapi_0.10              
##  [7] bit64_0.9-7                   interactiveDisplayBase_1.23.0
##  [9] codetools_0.2-16              splines_3.6.0                
## [11] Formula_1.2-3                 cluster_2.0.9                
## [13] shiny_1.3.2                   BiocManager_1.30.4           
## [15] compiler_3.6.0                httr_1.4.0                   
## [17] backports_1.1.4               assertthat_0.2.1             
## [19] Matrix_1.2-17                 lazyeval_0.2.2               
## [21] later_0.8.0                   acepack_1.4.1                
## [23] htmltools_0.3.6               prettyunits_1.0.2            
## [25] tools_3.6.0                   gtable_0.3.0                 
## [27] glue_1.3.1                    GenomeInfoDbData_1.2.1       
## [29] dplyr_0.8.1                   rappdirs_0.3.1               
## [31] Rcpp_1.0.1                    ExperimentHub_1.11.1         
## [33] xfun_0.7                      stringr_1.4.0                
## [35] mime_0.6                      ensembldb_2.9.1              
## [37] statmod_1.4.30                XML_3.98-1.19                
## [39] AnnotationHub_2.17.3          zlibbioc_1.31.0              
## [41] scales_1.0.0                  BSgenome_1.53.0              
## [43] VariantAnnotation_1.31.3      ProtGenerics_1.17.2          
## [45] hms_0.4.2                     promises_1.0.1               
## [47] AnnotationFilter_1.9.0        RColorBrewer_1.1-2           
## [49] yaml_2.2.0                    curl_3.3                     
## [51] memoise_1.1.0                 gridExtra_2.3                
## [53] ggplot2_3.1.1                 biomaRt_2.41.0               
## [55] rpart_4.1-15                  latticeExtra_0.6-28          
## [57] stringi_1.4.3                 RSQLite_2.1.1                
## [59] highr_0.8                     checkmate_1.9.3              
## [61] rlang_0.3.4                   pkgconfig_2.0.2              
## [63] bitops_1.0-6                  evaluate_0.13                
## [65] lattice_0.20-38               purrr_0.3.2                  
## [67] GenomicAlignments_1.21.2      htmlwidgets_1.3              
## [69] bit_1.1-14                    tidyselect_0.2.5             
## [71] plyr_1.8.4                    magrittr_1.5                 
## [73] bookdown_0.10                 R6_2.4.0                     
## [75] Hmisc_4.2-0                   DBI_1.0.0                    
## [77] pillar_1.4.0                  foreign_0.8-71               
## [79] survival_2.44-1.1             RCurl_1.95-4.12              
## [81] nnet_7.3-12                   tibble_2.1.1                 
## [83] crayon_1.3.4                  KernSmooth_2.23-15           
## [85] rmarkdown_1.13                progress_1.2.2               
## [87] locfit_1.5-9.1                data.table_1.12.2            
## [89] blob_1.1.1                    digest_0.6.19                
## [91] xtable_1.8-4                  httpuv_1.5.1                 
## [93] munsell_0.5.0
```

# References

