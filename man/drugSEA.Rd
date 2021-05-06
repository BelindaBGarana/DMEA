\name{drugSEA}
\alias{drugSEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Drug Set Enrichment Analysis
%%  ~~function to do ... ~~
}
\description{Using correlation estimates with drug sensitivity and drug set annotations, drugSEA determines the enrichment of drug sets.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
drugSEA(data, gmt, drug="Drug", estimate="Spearman.est", set.type=NULL, direction.adjust=NULL,FDR=0.05, num.permutations=1000, stat.type="Weighted")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Dataframe containing correlation estimates and moa for each drug.
%%     ~~Describe \code{x} here~~
}
  \item{gmt}{R object read in as a .gmt file containing mechanism-of-action sets of drugs.
%%     ~~Describe \code{x} here~~
}
  \item{drug}{Column name in data containing the drug names. Default: "Drug".
%%     ~~Describe \code{x} here~~
}
  \item{estimate}{Column name in data containing the correlation estimates. Default: "Spearman.est".
%%     ~~Describe \code{x} here~~
}
  \item{set.type}{Column name in data containing the drug set descriptions. Optional unless using direction.adjust feature.
%%     ~~Describe \code{x} here~~
}
\item{direction.adjust}{List of conditions for set.type to identify which, if any, estimates to multiply by -1. Default: NULL.
%%     ~~Describe \code{x} here~~
}
  \item{FDR}{False discovery rate threshold. Default: 0.05.
%%     ~~Describe \code{x} here~~
}
  \item{num.permutations}{Number of permutations to perform. Default: 1000.
%%     ~~Describe \code{x} here~~
}
  \item{stat.type}{Character string set to "Weighted" (weight = 1) or "Classic" (score weight = 0). Default: "Weighted".
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{result }{Dataframe containing GSEA results for drug sets.}
  \item{mtn.plots }{ggplot object list of mountain plots for drug sets with enrichment which meet the FDR threshold.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Belinda B. Garana, James H. Joly
%%  ~~who you are~~
}
\note{Mountain plots can be saved individually out of a list i using ggsave(filename, drugSEA$mtn.plots[i]).
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Prepare dataframe containing correlation estimates and mechanism-of-action for each drug
Drug <- LETTERS[seq(from = 1, to = 26)]
Spearman.est <- seq(-1,1, by = 0.077)
moa <- LETTERS[seq(from = 1, to = 26)]
corr.results <- as.data.frame(cbind(Drug, Spearman.est, moa))

#Prepare .gmt object
drug.set.A <- c("A","B","C","D","E")
drug.set.B <- c("V","W","X","Y","Z")
drug.set.X <- c("J","Q","O","E","F")
drug.set.Y <- c("D","S","K","L","R")
drug.set.Z <- c("G","W","P","B","T")
names <- c("drug.set.A","drug.set.B","drug.set.X","drug.set.Y","drug.set.Z")
Letter.Sets <- list(drug.set.A,drug.set.B,drug.set.X,drug.set.Y,drug.set.Z)
names(Letter.Sets) <- NULL
Drug.Sets <- list(genesets = Letter.Sets, geneset.names = names, geneset.descriptions = names)
Drug.Sets <- list(genesets = Letter.Sets, geneset.names = names)

#Perform drugSEA
drugSEA(data=corr.results, gmt=Drug.Sets)
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
