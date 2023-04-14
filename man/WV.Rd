\name{WV}

\alias{WV}

\title{Weighted Voting
}

\description{Weights are multiplied by their expression values, then summed
across all parameters (e.g., genes) for each sample or cell line to produce a
score. Various types of expression data may be used (e.g., RNAseq, proteomics
data).
}

\usage{
WV(expression, weights, sample.names=colnames(expression)[1],
        gene.names=colnames(weights)[1], weight.values=colnames(weights)[2])
}

\arguments{
  \item{expression}{Data frame containing normalized gene expression with gene
  symbols as column names, except for one column which contains sample names.
  }

  \item{weights}{Data frame with gene.names and weight.values (e.g.,
  log2(fold-change) values). All gene.names must be present in the expression
  data frame.
  }

  \item{sample.names}{Column name in expression data frame containing the sample
  names. Default: colnames(expression)[1].
  }

  \item{gene.names}{Column name in weights data frame containing the gene
  symbols. Default: colnames(weights)[1].
  }

  \item{weight.values}{Column name in weights data frame containing the weight
  values (e.g., log2(fold-change) values). Default: colnames(weights)[2].
  }
}

\value{
  \item{scores}{Data frame with Weighted Voting scores for each sample name.}

  \item{unused.weights}{Data frame with input weights not used for Weighted
  Voting.}
}

\references{
For more information, please visit: https://belindabgarana.github.io/DMEA
}

\author{Belinda B. Garana
}

\section{Warning}{Make sure your weights have the same gene names as your
expression data frame, and that your expression data frame is formatted with gene
symbols along columns.}

\examples{
## Step 1: prepare gene weight data frame
# create list of gene symbols
Gene <- paste0("Gene_", seq(from = 1, to = 50))

# by default, gene symbols are found in
# the first column of your weights data frame
weights <- as.data.frame(Gene)

# give each gene a weight
# by default, gene weight values are found in
# the second column of your weights data frame
weights$Rank_metric <- rnorm(length(Gene))

## Step 2: prepare expression data frame
# by default, column 1 of your expression data frame is
# the column name from which sample names are gathered and
# the column containing sample names in your drug sensitivity
# data frame should have the same name

# create list of sample names
Sample_ID <- seq(from = 1, to = 21)

# by default, column 1 of your expression data frame is
# the column name from which sample names are gathered
expression <- as.data.frame(Sample_ID)

# give each gene values representative of normalized RNA expression
# each gene is represented by a column in your expression data frame
for(i in 1:length(Gene)){
  expression[,c(Gene[i])] <- rnorm(length(Sample_ID), sd = 0.5)
}

## Step 3: run Weighted Voting to score each sample
WV.test <- DMEA::WV(expression, weights)
}

