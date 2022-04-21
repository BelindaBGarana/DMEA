# DMEA
Drug Mechanism Enrichment Analysis (DMEA) measures the enrichment of drug mechanisms for Class A vs. Class B samples. Any type of ranking metric for drugs is appropriate, but we recommend using Pearson correlation estimates between Weighted Voting score and drug AUC. To replicate our analysis, you will need to input a list of gene weights, a gmt object for drug sets, as well as a gene expression matrix and drug sensitivity matrix which can be obtained from CCLE and PRISM databases, respectively. Inputs except your gene weights of interest and CCLE RNAseq data (due to size) are available in the "Inputs" folder. You can read more about DMEA here: https://www.biorxiv.org/content/10.1101/2022.03.15.484520v1.

This package contains 4 functions: 1 function which can perform everything described above (DMEA) and 3 others which comprise the 3 custom steps of the DMEA function such that each can be used as a ranking method in their own right.

DMEA usually takes just 1 minute with 1,351 drugs spanning 85 drug mechanisms-of-action sets from the PRISM dataset and 300+ cell line samples.

To install the DMEA R package, run the next 2 lines in your command line:

if (!require(devtools)){ install.packages(dev.tools) }

devtools::install_github('BelindaBGarana/DMEA')

Summary of functions:

-WV: used to rank samples based on molecular subtype. Inputs: expression dataframe and weight list (e.g., log2(fold-change) for RNA expression between two classes).

-rank.corr: used to run correlations & regressions between 1 rank list and 2+ other rank lists. Also makes scatter plots for correlations which pass your FDR threshold.

-drugSEA: used to measure enrichment of drug sets based on drug rank (typically drug AUC or correlation estimate with WV score). Also makes mountain plots for enrichments which pass your FDR threshold. Direction-adjustment is available if needed (e.g., if drugs are grouped by target or another tag which does not distinguish directionality like drug mechanism-of-action does).

-DMEA: performs the above 3 functions all-in-one for drug mechanism-of-action sets.

