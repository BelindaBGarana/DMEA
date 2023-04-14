\name{DMEA}

\alias{DMEA}

\title{Drug Mechanism Enrichment Analysis
}

\description{Ranks samples using Weighted Voting,
correlates Weighted Voting scores to drug sensitivity,
and then runs enrichment analysis using correlation estimates to
determine the enrichment of drug mechanisms-of-action.
}

\usage{
DMEA(drug.sensitivity, gmt=NULL, expression, weights, value="AUC",
    sample.names=colnames(expression)[1], gene.names=colnames(weights)[1],
    weight.values=colnames(weights)[2], rank.metric="Pearson.est", FDR=0.25,
    num.permutations=1000, stat.type="Weighted", drug.info=NULL, drug="Drug",
    set.type="moa", min.per.set=6, sep="[|]",
    exclusions=c("-666", "NA", "na", "NaN", "NULL"), descriptions=NULL,
    min.per.corr=3, scatter.plots=TRUE, scatter.plot.type="pearson",
    FDR.scatter.plots=0.05, xlab="Weighted Voting Score", ylab=value,
    position.x="min", position.y="min", se=TRUE)
}

\arguments{
  \item{drug.sensitivity}{Data frame containing drug sensitivity scores
  (e.g., AUC) with drug names as column names,
  except for one column which contains sample names.
  }

  \item{gmt}{Optional: gmt object containing drug sets (e.g., output from
  GSA.read.gmt function from GSA package). If gmt is not provided, drug.info
  parameter must be provided. Default: NULL.
  }

  \item{expression}{Data frame containing normalized gene expression with gene
  symbols as column names, except for one column which contains sample names.
  }

  \item{weights}{Data frame with gene.names and weight.values (e.g.,
  log2(fold-change) values).
  }

  \item{value}{Description of the values starting at column 3. Do not include
  spaces! Default: "AUC".
  }

  \item{sample.names}{Column name in drug sensitivity and expression data
  frames containing the sample names. Default: colnames(expression)[1].
  }

  \item{gene.names}{Column name in weights data frame containing gene symbols.
  Default: colnames(weights)[1].
  }

  \item{weight.values}{Column name in weights data frame containing weight
  values (e.g., log2(fold-change) values). Default: colnames(weights)[2].
  }

  \item{rank.metric}{Estimate for ranking drugs for enrichment analysis. Must
  be "Pearson.est", "Spearman.est", "Slope", "Rank.slope", or another column
  name in the output produced by the rank_corr function. Default: "Pearson.est".
  }

  \item{FDR}{False discovery rate cut-off value to produce mountain plots for
  significant enrichments. Default: 0.25.
  }

  \item{num.permutations}{Number of permutations to perform for enrichment
  analysis. Default: 1000.
  }

  \item{stat.type}{Statistic type for enrichment analysis. Must be either
  "Weighted" (score weight = 1) or "Classic" (score weight = 0). Default:
  "Weighted".
  }

  \item{drug.info}{Data frame with set membership information for each element.
  Optional unless gmt is not provided. Default: NULL.
  }

  \item{drug}{Column name in drug.info containing the drug names. Optional
  unless gmt is not provided. Default: "Drug".
  }

  \item{set.type}{Column name in drug.info containing the set membership
  information for each drug. Optional unless gmt is not provided. Default:
  "moa".
  }

  \item{min.per.set}{The minimum number of drugs required per set. Default: 6.
  }

  \item{sep}{Character(s) to use for splitting drug set names. Optional unless
  gmt is not provided. Default: "[|]".
  }

  \item{exclusions}{A vector of drug set names to be excluded from the generated
  gmt. Optional unless gmt is not provided. Default: c("-666", "NA", "na", "NaN",
  "NULL").
  }

  \item{descriptions}{Optional: column name in drug.info containing
  descriptions for each drug set. If not provided and gmt is also not provided,
  set.type will be used for the descriptions instead. Default: NULL.
  }

  \item{min.per.corr}{Minimum number of samples required to run correlations.
  Default: 3.
  }

  \item{scatter.plots}{Logical operator describing whether or not to produce
  scatter plots for correlations. Default: TRUE.
  }

  \item{scatter.plot.type}{Correlation type for plots. Must be specified as
  "pearson" or "spearman". Default: "pearson".
  }

  \item{FDR.scatter.plots}{False discovery rate cut-off value to produce
  scatter plots for significant correlations. Default: 0.05.
  }

  \item{xlab}{Text describing rank metric for x-axis label of scatter plots.
  Default: "Weighted Voting Score".
  }

  \item{ylab}{Text describing variable values for y-axis label of scatter
  plots. Default: value.
  }

  \item{position.x}{X-axis position of correlation estimate on scatter plots.
  Must be "min", "mid", "max", or numeric. Default: "mid".
  }

  \item{position.y}{Y-axis position of correlation estimate on scatter plots.
  Must be "min", "mid", "max", or numeric. Default: "max".
  }

    \item{se}{Logical operator describing whether or not to include the 95
    percent confidence interval limits of the regression line on the scatter
    plots. Default: TRUE.
    }
}

\value{
  \item{WV.scores}{Data frame with Weighted Voting scores for each sample.}

  \item{unused.weights}{Data frame with input weights not used for Weighted
  Voting.}

  \item{corr.result}{Data frame with Pearson & Spearman correlation
  estimates, p-values, and FDR q-values (Benjamini-Hochberg method) for each
  variable item.}

  \item{corr.scatter.plots}{ggplot object with a scatter plot for each
  significant correlation on each page.}

  \item{gmt}{gmt object containing list of elements (genesets), set names
  (geneset.names), and set descriptions (geneset.descriptions).}

  \item{result}{Data frame containing GSEA results for drug sets.}

  \item{mtn.plots}{ggplot object list of mountain plots for drug sets with
  enrichment which meet the FDR threshold.}

  \item{volcano.plot}{ggplot object for volcano plot summarizing enrichment of
  all drug sets evaluated.}

  \item{removed.sets}{Data frame containing any drug sets which had too few
  elements to be evaluated and their number of elements.}

  \item{unannotated.drugs}{Data frame containing input drug names which were
  not matched to a drug set and their ranks.}
}

\references{
For more information, please visit: https://belindabgarana.github.io/DMEA
}

\author{Belinda B. Garana, James H. Joly
}

\note{Mountain plots can be saved individually using :
for(i in seq_len(length(DMEA.output$mtn.plots))){
ggsave(filename, DMEA.output$mtn.plots[[i]])
}

Output gmt can be saved for future use as an R object using:
saveRDS(DMEA.output$gmt, file)
}

\section{Warning}{The column name for your samples must be the same in your
drug sensitivity data frame as in your expression data frame. You can set the
name of the column containing your samples using the parameter sample.names.

If no gmt is provided, then drug.info with set.type annotations must be
provided. Either the gmt or set.type column must provide annotations for 2+
drug sets.

Each drug must have only one numeric, nonzero rank.
}

\seealso{
  \link[DMEA]{WV}
  \link[DMEA]{rank_corr}
  \link[DMEA]{as_gmt}
  \link[DMEA]{drugSEA}
}

\examples{
## Step 1: prepare drug sensitivity data frame
# create list of sample names
Sample_ID <- seq(from = 1, to = 21)

drug.sensitivity <- as.data.frame(Sample_ID)

# create list of drug names
Drug <- paste0("Drug_", seq(from = 1, to = 24))

# give each drug values representative of AUC sensitivity scores
for(i in 1:length(Drug)){
  drug.sensitivity[,c(Drug[i])] <- rnorm(length(Sample_ID),
                                        mean = 0.83, sd = 0.166)
}

## Step 2: prepare drug info data frame
# alternatively, a gmt object could be provided
info <- as.data.frame(Drug)

# moa is our default name for the column containing drug set annotations
info$moa <- rep(paste("Set", LETTERS[seq(from = 1, to = 4)]), 6)

## Step 3: prepare gene weight data frame
# create list of gene symbols
Gene <- paste0("Gene_", seq(from = 1, to = 50))

# by default, gene symbols are found in
# the first column of your weights data frame
weights <- as.data.frame(Gene)

# give each gene a weight
# by default, gene weight values are found in
# the second column of your weights data frame
weights$Rank_metric <- rnorm(length(Gene))

## Step 4: prepare expression data frame
# by default, column 1 of your expression data frame is
# the column name from which sample names are gathered and
# the column containing sample names in your drug sensitivity
# data frame should have the same name
expr <- as.data.frame(Sample_ID)

# give each gene values representative of normalized RNA expression
# each gene is represented by a column in your expression data frame
for(i in 1:length(Gene)){
  expr[,c(Gene[i])] <- rnorm(length(Sample_ID), sd = 0.5)
}

## Step 5: perform DMEA and store results
DMEA.test <- DMEA::DMEA(drug.sensitivity, expression = expr,
                        weights = weights, drug.info = info)
}
