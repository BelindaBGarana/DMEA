
# Drug Mechanism Enrichment Analysis

## Introduction

Though there are many algorithms for drug repurposing, most only
evaluate individual drugs which can be prone to off-target effects. In
contrast, our approach is to evaluate drug mechanisms-of-action by
aggregating results across many drugs. Though mechanisms-of-action have
also been evaluated in a similar manner by the Connectivity Map
(clue.io), Drugmonizome (<https://maayanlab.cloud/drugmonizome>), and
DrugEnrichr(<https://maayanlab.cloud/DrugEnrichr>), these tools only
allow analysis of select public datasets and do not consider the ranks
or directionality of input features. With our tool Drug Mechanism
Enrichment Analysis (DMEA), users can query any dataset of interest
using either an input gene signature or drug rank list to identify
enriched mechanisms-of-action for drug repurposing efforts. For more
information, please visit our website:
<https://belindabgarana.github.io/DMEA>

## Installation

To install this package, start R and enter:

``` r
if (!require("BiocManager", quietlly = TRUE))
  install.packages("BiocManager")
  
BiocManager::install("DMEA")
```

Alternatively, to install this package directly from our GitHub
repository, start R and enter:

``` r
if (!require("devtools", quietlly = TRUE))
  install.packages("devtools")
  
devtools::install_github("BelindaBGarana/DMEA")
```

If you are using Windows OS, you may need to change the above code to:

``` r
if (!require("devtools", quietlly = TRUE))
  install.packages("devtools")
  
devtools::install_github("BelindaBGarana/DMEA", build = FALSE)
```

## Input Option 1: Drug Rank List

With our drugSEA function, you can calculate the enrichment of drug sets
in an input drug rank list. This drug rank list can be derived from
anywhere (e.g., various drug repurposing algorithms). The only
requirements are that: \* Each drug has only one numeric, nonzero rank
metric \* At least 2 drug sets are represented in your list

Required inputs: \* data: a dataframe containing one column for drug
names and another column for drug ranks \* gmt: gmt object containing
drug set information. If not provided, a third column in your data input
parameter must include set annotations; if so, the name of this column
must be denoted by the set.type parameter. If you do not have your own
drug set annotations, you can use our PRISM moa annotations on our
GitHub repository at:
<https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/MOA_gmt_file_n6_no_special_chars.gmt>

Example: 1. Load the DMEA package into your R environment

``` r
library(DMEA)
```

2.  Import a data frame containing your drug rank list

``` r
# prepare dataframe containing rank metrics and set annotations for each drug
Drug <- paste0("Drug_", seq(from = 1, to = 24)) # Drug is our default for the drug parameter
data <- as.data.frame(Drug)
data$Pearson.est <- rnorm(length(Drug), sd = 0.25) # Pearson.est is our default for the rank.metric parameter
data$moa <- rep(paste("Set", LETTERS[seq(from = 1, to = 4)]), 6) # moa is our default for the set.type parameter

# perform drugSEA and store results
# if you use other column names for your data parameter, 
# then make sure to set drug, rank.metric, and set.type parameters
drugSEA.test <- DMEA::drugSEA(data)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## 
    ## Attaching package: 'snow'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, clusterSplit, makeCluster, parApply,
    ##     parCapply, parLapply, parRapply, parSapply, splitIndices,
    ##     stopCluster

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## 
    ## Attaching package: 'testthat'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     matches

    ## [1] "Generating gmt object for enrichment analysis..."
    ## [1] "Running enrichment analysis..."

    ## Warning in DMEA::drugSEA(data): No enrichments met the FDR cut-off to produce
    ## mountain plots

## Input Option 2: Gene Signature

With our DMEA function, you can identify selectively toxic drug sets
based on the correlations between drug sensitivity scores (e.g., AUC)
and weighted gene voting scores generated with an input gene signature.
This requires having drug screen and expression datasets which share the
same samples (e.g., cell lines). To run this analysis, it is required
that: \* Each drug has only one numeric, nonzero correlation metric \*
At least 2 drug sets are represented in your drug screen data \* Each
gene has only one numeric weight value \* At least one gene symbol
matches between the input gene signature and expression data frames \*
At least 3 samples (e.g., cell lines) are in both the expression and
drug screen sensitivity score data frames

Required inputs: \* drug.sensitivity: a dataframe containing drug
sensitivity scores (e.g., AUC) with drug names as column names, except
for one column which contains sample names \* gmt: gmt object containing
drug set information. If not provided, a drug.info parameter with set
membership information must be provided. If you do not have your own
drug set annotations, you can use our PRISM moa annotations on our
GitHub repository at:
<https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/MOA_gmt_file_n6_no_special_chars.gmt>
\* expression: a dataframe containing normalize gene expression with
gene symbols as column names, except for one column which contains
sample names \* weights: a dataframe containing gene symbols in one
column (default: column 1) and weight values in another column (default:
column 2)

Example: 1. Load the DMEA package into your R environment

``` r
library(DMEA)
```

2.  Import a data frame containing your drug rank list

``` r
# prepare drug sensitivity data frame
Sample_ID <- seq(from = 1, to = 21) # create list of sample names
drug.sensitivity <- as.data.frame(Sample_ID)
Drug <- paste0("Drug_", seq(from = 1, to = 24)) # create list of drug names
for(i in 1:length(Drug)){ # give each drug values representative of AUC sensitivity scores
  drug.sensitivity[,c(Drug[i])] <- rnorm(length(Sample_ID), mean = 0.83, sd = 0.166)
}

# prepare drug.info data frame
# alternatively, a gmt object could be provided
drug.info <- as.data.frame(Drug)
drug.info$moa <- rep(paste("Set", LETTERS[seq(from = 1, to = 4)]), 6) # moa is our default for the set.type parameter

# prepare gene weight data frame
Gene <- paste0("Gene_", seq(from = 1, to = 50)) # create list of gene symbols
weights <- as.data.frame(Gene) # by 
weights$Rank_metric <- rnorm(length(Gene)) # give each gene a weight

# prepare expression data frame
# by default, column 1 of your expression data frame is the column name from which sample names are gathered
# the column containing sample names in your drug sensitivity data frame should have the same name
expression <- as.data.frame(Sample_ID)
for(i in 1:length(Gene)){ # give each gene values representative of normalized RNA expression
  expression[,c(Gene[i])] <- rnorm(length(Sample_ID), sd = 0.5)
}

# perform DMEA and store results
DMEA.test <- DMEA::DMEA(drug.sensitivity, expression = expression, weights = weights, drug.info = drug.info)
```

    ## [1] "Calculating Weighted Voting scores..."

    ## 
    ## Attaching package: 'tidyselect'

    ## The following object is masked from 'package:testthat':
    ## 
    ##     matches

    ## [1] "Running correlations and regressions..."

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## Warning in rank.corr(data = WV.result.drug.sensitivity, variable = "Drug", : No
    ## correlations met the FDR cut-off to produce scatter plots

    ## [1] "Generating gmt object for enrichment analysis..."
    ## [1] "Running enrichment analysis..."

    ## Warning in drugSEA(data = corr.output, gmt, rank.metric = rank.metric,
    ## stat.type = stat.type, : No enrichments met the FDR cut-off to produce mountain
    ## plots

## Session Info

``` r
sessionInfo()
```

    ## R version 4.2.3 (2023-03-15)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3    reshape2_1.4.4   qvalue_2.30.0    tidyselect_1.2.0
    ##  [5] ggrepel_0.9.3    aplot_0.1.10     cowplot_1.1.1    ggplot2_3.4.1   
    ##  [9] testthat_3.1.7   sjmisc_2.8.9     doSNOW_1.0.20    iterators_1.0.14
    ## [13] foreach_1.5.2    snow_0.4-4       dplyr_1.1.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] sjlabelled_1.2.0   xfun_0.38          splines_4.2.3      ggfun_0.0.9       
    ##  [5] colorspace_2.1-0   vctrs_0.6.1        generics_0.1.3     htmltools_0.5.5   
    ##  [9] yaml_2.3.7         utf8_1.2.3         gridGraphics_0.5-1 rlang_1.1.0       
    ## [13] pillar_1.9.0       glue_1.6.2         withr_2.5.0        plyr_1.8.8        
    ## [17] lifecycle_1.0.3    stringr_1.5.0      munsell_0.5.0      gtable_0.3.3      
    ## [21] codetools_0.2-19   evaluate_0.20      knitr_1.42         fastmap_1.1.1     
    ## [25] fansi_1.0.4        Rcpp_1.0.10        scales_1.2.1       desc_1.4.2        
    ## [29] pkgload_1.3.2      brio_1.1.3         digest_0.6.31      stringi_1.7.12    
    ## [33] insight_0.19.1     grid_4.2.3         rprojroot_2.0.3    DMEA_0.1.0        
    ## [37] cli_3.6.1          tools_4.2.3        yulab.utils_0.0.6  magrittr_2.0.3    
    ## [41] patchwork_1.1.2    tibble_3.2.1       pkgconfig_2.0.3    ggplotify_0.1.0   
    ## [45] rmarkdown_2.21     rstudioapi_0.14    R6_2.5.1           compiler_4.2.3
