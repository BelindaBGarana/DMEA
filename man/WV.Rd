\name{WV}
\alias{WV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Voting
%%  ~~function to do ... ~~
}
\description{Weights are multiplied by their gene expression values, then summed across all genes for each sample or cell line to produce a score. Other types of expression data may also be used (e.g., proteomics data).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
WV(expression, weights, sample.names=colnames(expression)[1],gene.names=colnames(weights)[1], weight.values=colnames(weights)[2])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression}{Dataframe containing normalized gene expression with gene names as column names, except for one column which contains sample names.
%%     ~~Describe \code{x} here~~
}
  \item{weights}{Dataframe with gene.names and weight.values (e.g., log2(fold-change) values). All gene.names must be present in the expression dataframe.
%%     ~~Describe \code{x} here~~
}
  \item{sample.names}{Column name in expression dataframe containing the sample names. Default: colnames(expression)[1].
%%     ~~Describe \code{x} here~~
}
  \item{gene.names}{Column name in weights dataframe containing the gene names. Default: colnames(weights)[1].
%%     ~~Describe \code{x} here~~
}
  \item{weight.values}{Column name in weights dataframe containing the weight values (e.g., log2(fold-change) values). Default: colnames(weights)[2].
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{scores }{Dataframe with Weighted Voting scores for each sample.}
  \item{matrix }{Dataframe with Weighted Voting values for each gene and sample.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
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
#create sample ids
id <- LETTERS[seq(from = 1, to = 26)]

#create gene expression for each sample id
gene.1 <- as.numeric(seq(0,1, by = 0.04))
gene.2 <- as.numeric(seq(0,1, by = 0.04))
gene.3 <- as.numeric(seq(0,1, by = 0.04))

#prepare expression dataframe
expression <- as.data.frame(cbind(id, gene.1, gene.2, gene.3)) #Example expression dataframe with 26 samples and 3 genes

#prepare weights dataframe
gene.names <- c("gene.1", "gene.2", "gene.3")
weights <- as.data.frame(gene.names)
weights$weight.value <- 1 #In this example, these genes all have a weight of 1

#run Weighted Voting to score each sample
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
