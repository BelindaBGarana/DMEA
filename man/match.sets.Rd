\name{match.sets}
\alias{match.sets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Direction-match 2+ rank-lists
%%  ~~function to do ... ~~
}
\description{Combine 2+ rank-lists based on matched direction of association.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
match.sets(data, directions, ranks, id.var=colnames(data)[1],
n=10, rows.top=NULL, rows.bottom=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Dataframe containing a merge of 2+ rank-lists and their directions.
%%     ~~Describe \code{x} here~~
}
  \item{directions}{Vector containing column names of directions for each set.
%%     ~~Describe \code{x} here~~
}
  \item{ranks}{Vector containing column names of ranks for each set.
%%     ~~Describe \code{x} here~~
}
  \item{id.var}{Column name containing id variable. Default: colnames(data)[1].
%%     ~~Describe \code{x} here~~
}
  \item{n}{Number of id variables to include in top or bottom lists. Default: 10.
%%     ~~Describe \code{x} here~~
}
  \item{rows.top}{Row variable for tableGrob to output top.table. Default: NULL.
%%     ~~Describe \code{x} here~~
}
  \item{rows.bottom}{Row variable for tableGrob to output bottom.table. Default: NULL.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{result }{Dataframe with average direction and average rank of each id variable across sets.}
  \item{top }{Dataframe with top direction-matched id variables, including their average rank and direction.}
  \item{bottom }{Dataframe with bottom direction-matched id variables, including their average rank and direction.}
  \item{top.table }{tableGrob of top direction-matched id variables, including their average rank and direction.}
  \item{bottom.table }{tableGrob of bottom direction-matched id variables, including their average rank and direction.}
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
{Warning }{Be sure that your rank columns contain integers.}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#example input: correlation between rank metric and gene expression
id <- LETTERS[seq(from = 1, to = 26)]
rank.metric <- as.numeric(seq(1,2, by = 0.04))
gene.1 <- as.numeric(seq(0,1, by = 0.04))
gene.2 <- as.numeric(seq(0,1, by = 0.04))
gene.3 <- as.numeric(seq(0,1, by = 0.04))
data <- as.data.frame(cbind(id, rank.metric, gene.1, gene.2, gene.3)) #Example input with 26 id.var and 3 genes
correlations <- rank.corr(data=data,id.var="id",rank.var="rank.metric",variable="Gene",value="Intensity",type="spearman",plots=FALSE)$result
correlations$Rank.Pearson <- rank(desc(as.numeric(correlations$Pearson.est)))
correlations$Rank.Spearman <- rank(desc(as.numeric(correlations$Spearman.est)))

#example analysis: match Pearson and Spearman correlation results
match.sets(data=correlations, directions=c("Pearson.est","Spearman.est"),
  ranks=c("Rank.Pearson","Rank.Spearman"))
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
