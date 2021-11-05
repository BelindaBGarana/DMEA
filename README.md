# DMEA
Drug Mechanism Enrichment Analysis (DMEA) measures the enrichment of drug mechanisms for Class A vs. Class B samples. Any type of ranking metric for drugs is appropriate, but we recommend using Pearson correlation estimates between Weighted Voting score and drug AUC. 

This package contains 4 functions: 1 function which can perform everything described above (DMEA) and 3 others which comprise the 3 custom steps of the DMEA function such that each can be used as a ranking method in their own right.

DMEA can take 1-5 minutes for 500 genes or less, but it can take 40+ minutes for 1,000+ genes. Using CCLE RNA expression and PRISM drug sensitivity data, we have obtained consistent results even if only using the top 500 genes based on absolute rank. If you have a large set of gene weights (1,000+), we recommend running DMEA with your top 500 gene weights.

To install the DMEA R package, run the next 2 lines in your command line:

if (!require(devtools)){ install.packages(dev.tools) }

devtools::install_github('BelindaBGarana/DMEA')

Summary of functions:

-WV: used to rank samples based on molecular subtype. Inputs: expression dataframe and weight list (e.g., log2(fold-change) for RNA expression between two classes).

-rank.corr: used to run correlations & regressions between 1 rank list and 2+ other rank lists. Also makes scatter plots for correlations which pass your FDR threshold.

-drugSEA: used to measure enrichment of drug sets based on drug rank (typically drug AUC or correlation estimate with WV score). Also makes mountain plots for enrichments which pass your FDR threshold. Direction-adjustment is available if needed (e.g., if drugs are grouped by target or another tag which does not distinguish directionality like drug mechanism-of-action does).

-DMEA: performs the above 3 functions all-in-one for drug mechanism-of-action sets.

