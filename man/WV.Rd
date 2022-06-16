\name{WV}
\alias{WV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Voting
%%  ~~function to do ... ~~
}
\description{Weights are multiplied by their expression values, then summed across all parameters (e.g., genes) for each sample or cell line to produce a score. Various types of expression data may be used (e.g., RNAseq, proteomics data).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
WV(expression, weights, sample.names=colnames(expression)[1],gene.names=colnames(weights)[1], weight.values=colnames(weights)[2])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression}{Dataframe containing normalized gene expression with gene names as column names, except for one column which contains sample names.
}
  \item{weights}{Dataframe with gene.names and weight.values (e.g., log2(fold-change) values). All gene.names must be present in the expression dataframe.
}
  \item{sample.names}{Column name in expression dataframe containing the sample names. Default: colnames(expression)[1].
}
  \item{gene.names}{Column name in weights dataframe containing the gene names. Default: colnames(weights)[1].
}
  \item{weight.values}{Column name in weights dataframe containing the weight values (e.g., log2(fold-change) values). Default: colnames(weights)[2].
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Dataframe with Weighted Voting scores for each sample name.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\references{
%% ~put references to the literature/web site here ~
}
\author{Belinda B. Garana
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{Make sure your weights have the same gene names as your expression dataframe, and that your expression dataframe is formatted with gene names along columns.}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# create sample ids
id <- LETTERS[seq(from = 1, to = 26)]

# create gene expression for each sample id
gene.1 <- as.numeric(seq(0,1, by = 0.04))
gene.2 <- as.numeric(seq(0,1, by = 0.04))
gene.3 <- as.numeric(seq(0,1, by = 0.04))

# prepare expression dataframe
expression <- as.data.frame(cbind(id, gene.1, gene.2, gene.3)) # example expression dataframe with 26 samples and 3 genes

# prepare weights dataframe
gene.names <- c("gene.1", "gene.2", "gene.3")
weights <- as.data.frame(gene.names)
weights$weight.value <- 1 # in this example, these genes all have a weight of 1

# run Weighted Voting to score each sample
WV(expression,weights)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
