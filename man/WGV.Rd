\name{WGV}
\alias{WGV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Gene Voting
%%  ~~function to do ... ~~
}
\description{Gene weights are multiplied by their expression values, then summed across genes for each sample or cell line to produce a score.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
WGV(data, weights, gene.names=colnames(weights)[1], weight.values=colnames(weights)[2])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Dataframe of gene expression with samples or cell lines along rows and gene names along columns.
%%     ~~Describe \code{x} here~~
}
  \item{weights}{Dataframe with gene.names and weight.values (such as log2(fold-change) values). All gene.names must be present in the data.
%%     ~~Describe \code{x} here~~
}
  \item{gene.names}{Column name containing the gene names. Default: colnames(weights)[1].
%%     ~~Describe \code{x} here~~
}
  \item{weight.values}{Column name containing the weight values (such as log2(fold-change) values). Default: colnames(weights)[2].
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{scores }{Dataframe with Weighted Gene Voting scores for each sample or cell line.}
  \item{matrix }{Dataframe with Weighted Gene Voting values for each gene and sample or cell line.}
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
{Warning }{Make sure your weights have the same gene names as your expression data, and that your expression data is formatted with gene names along columns.}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#create id
id <- LETTERS[seq(from = 1, to = 26)]

#create rank metric for each id
rank.metric <- as.numeric(seq(1,2, by = 0.04))

#create gene expression for each id
gene.1 <- as.numeric(seq(0,1, by = 0.04))
gene.2 <- as.numeric(seq(0,1, by = 0.04))
gene.3 <- as.numeric(seq(0,1, by = 0.04))

#prepare data
data <- as.data.frame(cbind(id, rank.metric, gene.1, gene.2, gene.3)) #Example input with 26 id.var and 3 genes

#prepare weights
gene.names <- c("gene.1", "gene.2", "gene.3")
weights <- as.data.frame(gene.names)
weights$weight.value <- 1

#run Weighted Gene Voting to score each id
WGV(data,weights)
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
