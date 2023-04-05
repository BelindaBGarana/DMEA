\name{WV}

\alias{WV}

\title{Weighted Voting
}

\description{Weights are multiplied by their expression values, then summed across all parameters (e.g., genes) for each sample or cell line to produce a score. Various types of expression data may be used (e.g., RNAseq, proteomics data).
}

\usage{
WV(expression, weights, sample.names=colnames(expression)[1],
          gene.names=colnames(weights)[1], weight.values=colnames(weights)[2])
}

\arguments{
  \item{expression}{Dataframe containing normalized gene expression with gene names as column names, except for one column which contains sample names.}

  \item{weights}{Dataframe with gene.names and weight.values (e.g., log2(fold-change) values). All gene.names must be present in the expression dataframe.}

  \item{sample.names}{Column name in expression dataframe containing the sample names. Default: colnames(expression)[1].}

  \item{gene.names}{Column name in weights dataframe containing the gene symbols. Default: colnames(weights)[1].}

  \item{weight.values}{Column name in weights dataframe containing the weight values (e.g., log2(fold-change) values). Default: colnames(weights)[2].}
}

\value{
  \item{scores}{Dataframe with Weighted Voting scores for each sample name.}

  \item{unused.weights}{Dataframe with input weights not used for Weighted Voting.}
}

\references{
For more information, please visit: https://belindabgarana.github.io/DMEA
}

\author{Belinda B. Garana
}

\section{Warning }{Make sure your weights have the same gene names as your expression dataframe, and that your expression dataframe is formatted with gene names along columns.}

\examples{
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

# run Weighted Voting to score each sample
WV.test <- DMEA::WV(expression, weights)
}

