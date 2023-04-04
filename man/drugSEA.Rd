\name{drugSEA}

\alias{drugSEA}

\title{Drug Set Enrichment Analysis
}

\description{Using drug rank metrics and drug set annotations, drugSEA calculates the enrichment of drug sets and produces plots. It can also generate a gmt object for future enrichment analyses.
}

\usage{
drugSEA(data, gmt=NULL, drug="Drug", rank.metric="Pearson.est", set.type="moa",
              direction.adjust=NULL, FDR=0.25, num.permutations=1000,
              stat.type="Weighted", min.per.set=6, sep = "[|]",
              exclusions = c("-666", "NA", "na", "NaN", "NULL"), descriptions=NULL,
              convert.synonyms=FALSE)
}

\arguments{
  \item{data}{Dataframe containing rank metric for each drug.}

  \item{gmt}{Optional: gmt object containing drug sets (e.g., output from GSA.read.gmt function from GSA package). If gmt is not provided, set.type parameter must match a column name in data. Default: NULL.}

  \item{drug}{Column name in data containing the drug names. Default: "Drug".}

  \item{rank.metric}{Column name in data containing the metric for ranking drugs. Default: "Pearson.est".}

  \item{set.type}{Column name in data containing the set membership information for each drug. Optional unless using direction.adjust feature or gmt is not provided. Default: "moa".}

\item{direction.adjust}{Vector of conditions for set.type to identify which, if any, estimates to multiply by -1 (e.g., c(" agonist", " activator")). Default: NULL.}

  \item{FDR}{False discovery rate cut-off value to produce mountain plots for significant enrichments. Default: 0.05.}

  \item{num.permutations}{Number of permutations to perform. Default: 1000.}

  \item{stat.type}{Kolmogorov-Smirnov statistic type. Must be either "Weighted" (weight = 1) or "Classic" (weight = 0). Default: "Weighted".}

  \item{min.per.set}{The minimum number of drugs required per set. Default: 6.}

  \item{sep}{Character(s) to use for splitting drug set names. Optional unless gmt is not provided. Default: "[|]".}

  \item{exclusions}{A vector of drug set names to be excluded from the generated gmt. Optional unless gmt is not provided. Default: c("-666", "NA", "na", "NaN", "NULL").}

  \item{descriptions}{Optional: column name in data containing descriptions for each drug set. If not provided and gmt is also not provided, set.type will be used for the descriptions instead. Default: NULL.}

  \item{convert.synonyms}{Logical operator describing whether or not to convert unannotated drug names to synonyms for drugs in the PRISM database. Default: FALSE.}
}

\value{
  \item{gmt }{gmt object containing list of elements (genesets), set names (geneset.names), and set descriptions (geneset.descriptions).}

  \item{result }{Dataframe containing enrichment analysis results for drug sets.}

  \item{mtn.plots }{ggplot object list of mountain plots for drug sets with enrichments which meet the FDR threshold.}

  \item{volcano.plot }{ggplot object for volcano plot summarizing enrichment of all drug sets evaluated.}

  \item{removed.sets }{Dataframe containing any drug sets which had too few elements to be evaluated and their number of elements.}

  \item{replaced.drugs }{Dataframe containing input drug names which were replaced based on PubChem synonyms and their PubChem CID.}

  \item{unannotated.drugs }{Dataframe containing input drug names which were not matched to a drug set and their ranks.}
}

\references{
For more information, please visit: https://belindabgarana.github.io/DMEA
}

\author{Belinda B. Garana, James H. Joly
}

\note{Mountain plots can be saved individually using for(i in 1:length(drugSEA.output$mtn.plots)){ggsave(filename, drugSEA.output$mtn.plots[[i]])}.
Output gmt can be saved for future use as an R object using saveRDS(DMEA.output$gmt, file).
}

\section{Warning }{Each drug must have 1 non-zero, numeric rank.metric value to be included in the analysis. NA and zero values will be omitted after coercing the rank.metric values to be of type "numeric". If no gmt is provided, then set.type must be provided in the input data. Either the gmt or set.type column must provide annotations for 2+ drug sets.

If convert.synonyms=TRUE, an internet connection is required and the drug synonym conversion is fastest if the input drug names are PubChem CIDs.
}

\examples{
## Example 1
# select example drug rank list
# this example has a column for moa annotations named "moa"
# otherwise, a gmt object would also be needed
selected.ex <- "Drug_rank_list/CMap/PRISM/Cell_lines_sensitive_to_HMGCR_inhibitor_lovastatin"
rank.data <- read.csv(file=paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Examples/",selected.ex,"/DMEA_input.csv"))

# perform drugSEA
drugSEA.test <- drugSEA(rank.data, drug = colnames(rank.data)[1],
                  rank.metric = colnames(rank.data)[2])

## Example 2
# select example drug rank list
selected.ex2 <- "Gene_signature/GSE31625_NSCLC_sensitive_vs_resistant_to_erlotinib"
rank.data2 <- read.csv(file=paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Examples/",selected.ex2,"/DMEA_correlation_results.csv"))

# import gmt containing PRISM moa annotations
gmt <- GSA::GSA.read.gmt(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/MOA_gmt_file_n6_no_special_chars.gmt")

# perform drugSEA
drugSEA.test2 <- drugSEA(rank.data2, gmt)
}
