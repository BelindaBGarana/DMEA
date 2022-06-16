\name{as.gmt}
\alias{as.gmt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gmt formatting
%%  ~~function to do ... ~~
}
\description{Formats sets into gmt object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
as.gmt(data, element.names = "Drug", set.names = "moa", min.per.set = 6, sep = "[|]", exclusions = c("-666", "NA", "NaN", "NULL"), descriptions = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Dataframe with set membership information for each element.
}
  \item{element.names}{Column name in data containing the element names. Default: "Drug".
}
  \item{set.names}{Column name in data containing the set names. Default: "moa".
}
  \item{min.per.set}{The minimum number of elements required per set for the output. Default: 6.
}
  \item{sep}{Character(s) to use for splitting set names. Default: "[|]".
}
  \item{exclusions}{A vector of set names to be excluded from the output. Default: c("-666", "NA", "NaN", "NULL").
}
  \item{descriptions}{Optional: column name in data containing descriptions for each set. If not provided, set names will be used instead. Default: NULL.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{gmt object containing list of elements (genesets), set names (geneset.names), and set descriptions (geneset.descriptions).
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
\note{The output can be used for enrichment analysis in R or saved for future use as an R object using saveRDS(gmt, file).
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# prepare input dataframe
Drug <- LETTERS[seq(from = 1, to = 18)]
moa <- c(rep(LETTERS[19], 6), rep(LETTERS[20], 6), rep(LETTERS[21], 6))
data <- as.data.frame(cbind(Drug, moa)) # this example input has 18 drugs and 3 sets

# produce gmt object
as.gmt(data)
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
