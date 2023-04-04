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

  \item{gene.names}{Column name in weights dataframe containing the gene names. Default: colnames(weights)[1].}

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
# download RNAseq data
download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin",
              destfile = paste0(getwd(),"/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"))
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin",
              destfile = paste0(getwd(),"/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"))
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
RNA.df <- rbind(RNA.first200, RNA.rest)

# select example gene signature
selected.ex <- "Gene_signature/GSE31625_NSCLC_sensitive_vs_resistant_to_erlotinib"
gene.sig <- read.csv(file=paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Examples/",selected.ex,"/Filtered_gene_signature_no_duplicates.csv"))

# run Weighted Voting to score each sample
WV.test <- WV(RNA.df, gene.sig)
}

