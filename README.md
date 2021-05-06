# DMEA
Drug Mechanism Enrichment Analysis -4 functions for drug set enrichment analysis, direction-matching rank lists, correlations across a matrix, and weighted gene voting:

drugSEA -Runs GSEA to identify enriched drug sets -Optional: Adjusts correlation estimates for direction of drug activity -Produces mountain plots for significant enrichments

match.sets -Direction-matches rank lists and calculates average rank -Produces tables for top, bottom matches

rank.corr -Runs Pearson & Spearman correlations between a list and a matrix -Produces scatterplots for significant correlations

WGV -Uses gene weights and a gene expression matrix to score each sample (weighted gene voting) -Produces dataframes with weighted gene voting scores and intermediate matrix results
