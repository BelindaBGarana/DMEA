\name{drugSEA}
\alias{drugSEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Drug Set Enrichment Analysis
%%  ~~function to do ... ~~
}
\description{Using drug rank metrics and drug set annotations, drugSEA calculates the enrichment of drug sets and produces plots. It can also generate a gmt object for future enrichment analyses.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
drugSEA(data, gmt, drug="Drug", rank.metric="Pearson.est", set.type="moa", direction.adjust=NULL,
                    FDR=0.25, num.permutations=1000, stat.type="Weighted", min.per.set=6,
                    sep = "[|]", exclusions = c("-666", "NA", "NaN", "NULL"), descriptions=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Dataframe containing rank metric for each drug.
}
  \item{gmt}{Optional: gmt object containing drug sets (e.g., output from GSA.read.gmt function from GSA package). If gmt is not provided, set.type parameter must match a column name in data. Default: NULL.
}
  \item{drug}{Column name in data containing the drug names. Default: "Drug".
}
  \item{rank.metric}{Column name in data containing the metric for ranking drugs. Default: "Pearson.est".
}
  \item{set.type}{Column name in data containing the set membership information for each drug. Optional unless using direction.adjust feature or gmt is not provided. Default: "moa".
}
\item{direction.adjust}{Vector of conditions for set.type to identify which, if any, estimates to multiply by -1 (e.g., c(" agonist", " activator")). Default: NULL.
}
  \item{FDR}{False discovery rate cut-off value to produce mountain plots for significant enrichments. Default: 0.05.
}
  \item{num.permutations}{Number of permutations to perform. Default: 1000.
}
  \item{stat.type}{Kolmogorov-Smirnov statistic type. Must be either "Weighted" (weight = 1) or "Classic" (weight = 0). Default: "Weighted".
}
  \item{min.per.set}{The minimum number of drugs required per set. Default: 6.
}
  \item{sep}{Character(s) to use for splitting drug set names. Optional unless gmt is not provided. Default: "[|]".
}
  \item{exclusions}{A vector of drug set names to be excluded from the generated gmt. Optional unless gmt is not provided. Default: c("-666", "NA", "NaN", "NULL").
}
  \item{descriptions}{Optional: column name in data containing descriptions for each drug set. If not provided and gmt is also not provided, set.type will be used for the descriptions instead. Default: NULL.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{gmt }{gmt object containing list of elements (genesets), set names (geneset.names), and set descriptions (geneset.descriptions).}
  \item{result }{Dataframe containing enrichment analysis results for drug sets.}
  \item{mtn.plots }{ggplot object list of mountain plots for drug sets with enrichments which meet the FDR threshold.}
  \item{volcano.plot }{ggplot object for volcano plot summarizing enrichment of all drug sets evaluated.}
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
\note{Mountain plots can be saved individually using for(i in 1:length(drugSEA.output$mtn.plots)){ggsave(filename, drugSEA.output$mtn.plots[[i]])}.
Output gmt can be saved for future use as an R object using saveRDS(DMEA.output$gmt, file).
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# prepare dataframe containing rank metrics (e.g., Pearson.est) for each drug
Drug <- LETTERS[seq(from = 1, to = 26)]
Pearson.est <- seq(-1,1, by = 0.077)
data <- as.data.frame(cbind(Drug, Pearson.est))

# prepare gmt object
drug.set.A <- c("A","B","C","D","E")
drug.set.B <- c("V","W","X","Y","Z")
drug.set.X <- c("J","Q","O","E","F")
drug.set.Y <- c("D","S","K","L","R")
drug.set.Z <- c("G","W","P","B","T")
names <- c("drug.set.A","drug.set.B","drug.set.X","drug.set.Y","drug.set.Z")
Drug.Sets <- list(drug.set.A,drug.set.B,drug.set.X,drug.set.Y,drug.set.Z)
names(Drug.Sets) <- NULL
gmt <- list(genesets = Drug.Sets, geneset.names = names, geneset.descriptions = names)

# perform drugSEA
drugSEA(data, gmt)
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
