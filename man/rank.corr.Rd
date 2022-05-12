\name{rank.corr}
\alias{rank.corr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Correlations between a rank metric and a set of 2+ variables
%%  ~~function to do ... ~~
}
\description{Runs correlations (Pearson & Spearman) and regressions (value- & rank-based) between a rank metric and a set of 2+ variables as well as statistical analysis. Outputs a dataframe of results, as well as scatter plots for the correlation type and FDR-cutoff specified.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rank.corr(data, variable="Gene", value="Intensity",type="pearson", N.min=3, plots=TRUE, FDR=0.05, xlab=colnames(data)[2], ylab=value, position.x="mid", position.y="max", se=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Dataframe with sample names in column 1, the rank metric in column 2, and the set of variables (e.g. genes) as column names starting at column 3.
}
  \item{variable}{Description of the column names starting at column 3. Default: "Gene".
}
  \item{value}{Description of the values starting at column 3. Do not include spaces! Default: "Intensity".
}
  \item{type}{Correlation type for plots. Must be specified as "pearson" or "spearman". Default: "pearson".
}
  \item{N.min}{Minimum number of samples required to run correlations. Default: 3.
}
  \item{plots}{Set TRUE to produce scatter plots. Default: TRUE.
}
  \item{FDR}{False discovery rate cut-off value to produce scatter plots for significant correlations. Default: 0.05.
}
  \item{xlab}{Text describing rank metric for x-axis label of scatter plots. Default: colnames(data)[2].
}
  \item{ylab}{Text describing variable values for y-axis label of scatter plots. Default: value.
}
  \item{position.x}{X-axis position of correlation estimate on scatter plots. Must be "min", "mid", "max", or numeric. Default: "mid".
}
  \item{position.y}{Y-axis position of correlation estimate on scatter plots. Must be "min", "mid", "max", or numeric. Default: "max".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{result }{Dataframe with correlation and regression estimates, p-values, and FDR q-values (Benjamini-Hochberg method) for each variable item.}
  \item{scatter.plots }{ggplot object with a scatter plot for each significant correlation on each page.}
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
\note{Scatter plots can be saved collectively to 1 file using ggsave(filename, rank.corr$scatter.plots).
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{For this function to work, you need 2 or more variables in your set.
Do not name your variable as "variable.name" or your value as "value.name". Do not include spaces in "value" variable.}

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

#prepare input
data <- as.data.frame(cbind(id, rank.metric, gene.1, gene.2, gene.3)) #Example input with 26 id.var and 3 genes

#run correlations between rank metric and each gene's expression
rank.corr(data)
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
