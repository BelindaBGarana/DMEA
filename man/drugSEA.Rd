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
drugSEA(data, gmt, drug="Drug", estimate="Pearson.est", set.type="moa", direction.adjust=NULL,
                    FDR=0.25, num.permutations=1000, stat.type="Weighted")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Dataframe containing correlation estimates and moa for each drug.
}
  \item{gmt}{gmt list containing drug sets (e.g., output from GSA.read.gmt function from GSA package).
}
  \item{drug}{Column name in data containing the drug names. Default: "Drug".
}
  \item{estimate}{Column name in data containing the estimate for ranking drugs for enrichment analysis. Default: "Pearson.est".
}
  \item{set.type}{Column name in data containing the drug set descriptions. Optional unless using direction.adjust feature.
}
\item{direction.adjust}{List of conditions for set.type to identify which, if any, estimates to multiply by -1. Default: NULL.
}
  \item{FDR}{False discovery rate cut-off value to produce mountain plots for significant enrichments. Default: 0.05.
}
  \item{num.permutations}{Number of permutations to perform. Default: 1000.
}
  \item{stat.type}{Statistic type for enrichment analysis. Must be either "Weighted" (weight = 1) or "Classic" (score weight = 0). Default: "Weighted".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{result }{Dataframe containing enrichment analysis results for drug sets.}
  \item{mtn.plots }{ggplot object list of mountain plots for drug sets with enrichments which meet the FDR threshold.}
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
#Prepare dataframe containing correlation estimates for each drug
Drug <- LETTERS[seq(from = 1, to = 26)]
Pearson.est <- seq(-1,1, by = 0.077)
corr.results <- as.data.frame(cbind(Drug, Pearson.est))

#Prepare gmt object
drug.set.A <- c("A","B","C","D","E")
drug.set.B <- c("V","W","X","Y","Z")
drug.set.X <- c("J","Q","O","E","F")
drug.set.Y <- c("D","S","K","L","R")
drug.set.Z <- c("G","W","P","B","T")
names <- c("drug.set.A","drug.set.B","drug.set.X","drug.set.Y","drug.set.Z")
Drug.Sets <- list(drug.set.A,drug.set.B,drug.set.X,drug.set.Y,drug.set.Z)
names(Drug.Sets) <- NULL
gmt <- list(genesets = Drug.Sets, geneset.names = names, geneset.descriptions = names)

#Perform drugSEA
drugSEA(data=corr.results, gmt)
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
