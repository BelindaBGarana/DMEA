\name{DMEA}
\alias{DMEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Drug Mechanism Enrichment Analysis
%%  ~~function to do ... ~~
}
\description{Ranks samples using Weighted Voting, correlates Weighted Voting scores to drug sensitivity, and then runs enrichment analysis using correlation estimates to determine the enrichment of drug mechanisms-of-action.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DMEA(drug.sensitivity, gmt=NULL, expression, weights, value="AUC", expr.sample.names=colnames(expression)[1],
                 gene.names=colnames(weights)[1], weight.values=colnames(weights)[2],
                 rank.metric="Pearson.est", FDR=0.25, num.permutations=1000, stat.type="Weighted",
                 drug.info=NULL, drug="Drug", set.type="moa", min.per.set=6, sep="[|]", exclusions=c("-666", "NA", "na", "NaN", "NULL"), descriptions=NULL,
                 min.per.corr=3, scatter.plots=TRUE, scatter.plot.type="pearson", FDR.scatter.plots=0.05,
                 xlab="Weighted Voting Score", ylab=value, position.x="min", position.y="min", se=TRUE)
}
\arguments{
  \item{drug.sensitivity}{Dataframe containing drug sensitivity scores (e.g., AUC) with drug names as column names, except for one column which contains sample names.
}
  \item{gmt}{Optional: gmt object containing drug sets (e.g., output from GSA.read.gmt function from GSA package). If gmt is not provided, drug.info parameter must be provided. Default: NULL.
}
  \item{expression}{Dataframe containing normalized gene expression with gene names as column names, except for one column which contains sample names.
}
  \item{weights}{Dataframe with gene.names and weight.values (e.g., log2(fold-change) values). All gene.names must be present in the expression dataframe.
}
  \item{value}{Description of the values starting at column 3. Do not include spaces! Default: "AUC".
}
  \item{sample.names}{Column name in drug sensitivity and expression dataframes containing the sample names. Default: colnames(expression)[1].
}
  \item{gene.names}{Column name in weights dataframe containing the gene names. Default: colnames(weights)[1].
}
  \item{weight.values}{Column name in weights dataframe containing the weight values (e.g., log2(fold-change) values). Default: colnames(weights)[2].
}
  \item{estimate}{Estimate for ranking drugs for enrichment analysis. Must be "Pearson.est", "Spearman.est", "Slope", "Rank.slope", or another column name in the output typically produced by the rank.corr function. Default: "Pearson.est".
}
  \item{FDR}{False discovery rate cut-off value to produce mountain plots for significant enrichments. Default: 0.25.
}
  \item{num.permutations}{Number of permutations to perform for enrichment analysis. Default: 1000.
}
  \item{stat.type}{Statistic type for enrichment analysis. Must be either "Weighted" (score weight = 1) or "Classic" (score weight = 0). Default: "Weighted".
}
  \item{drug.info}{Dataframe with set membership information for each element. Optional unless gmt is not provided. Default: NULL.
}
  \item{drug}{Column name in drug.info containing the drug names. Optional unless gmt is not provided. Default: "Drug".
}
  \item{set.type}{Column name in drug.info containing the set membership information for each drug. Optional unless gmt is not provided. Default: "moa".
}
  \item{min.per.set}{The minimum number of drugs required per set. Default: 6.
}
  \item{sep}{Character(s) to use for splitting drug set names. Optional unless gmt is not provided. Default: "[|]".
}
  \item{exclusions}{A vector of drug set names to be excluded from the generated gmt. Optional unless gmt is not provided. Default: c("-666", "NA", "na", "NaN", "NULL").
}
  \item{descriptions}{Optional: column name in drug.info containing descriptions for each drug set. If not provided and gmt is also not provided, set.type will be used for the descriptions instead. Default: NULL.
}
  \item{min.per.corr}{Minimum number of samples required to run correlations. Default: 3.
}
  \item{scatter.plots}{Set TRUE to produce scatter plots for correlations. Default: TRUE.
}
  \item{scatter.plot.type}{Correlation type for plots. Must be specified as "pearson" or "spearman". Default: "pearson".
}
  \item{FDR.scatter.plots}{False discovery rate cut-off value to produce scatter plots for significant correlations. Default: 0.05.
}
  \item{xlab}{Text describing rank metric for x-axis label of scatter plots. Default: "Weighted Voting Score".
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
  \item{WV.scores }{Dataframe with Weighted Voting scores for each sample.}
  \item{corr.result }{Dataframe with Pearson & Spearman correlation estimates, p-values, and FDR q-values (Benjamini-Hochberg method) for each variable item.}
  \item{corr.scatter.plots }{ggplot object with a scatter plot for each significant correlation on each page.}
  \item{gmt }{gmt object containing list of elements (genesets), set names (geneset.names), and set descriptions (geneset.descriptions).}
  \item{result }{Dataframe containing GSEA results for drug sets.}
  \item{mtn.plots }{ggplot object list of mountain plots for drug sets with enrichment which meet the FDR threshold.}
  \item{volcano.plot }{ggplot object for volcano plot summarizing enrichment of all drug sets evaluated.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Belinda B. Garana, James H. Joly
%%  ~~who you are~~
}
\note{Mountain plots can be saved individually using for(i in 1:length(DMEA.output$mtn.plots)){ggsave(filename, DMEA.output$mtn.plots[[i]])}.
Output gmt can be saved for future use as an R object using saveRDS(DMEA.output$gmt, file).
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{The column name for your samples must be the same in your drug sensitivity dataframe as in your expression dataframe. You can set the name of the column containing your samples using the parameter sample.names.
If no gmt is provided, then drug.info with set.type annotations must be provided. Either the gmt or set.type column must provide annotations for 2+ drug sets.}

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
expression <- as.data.frame(cbind(id, gene.1, gene.2, gene.3)) #Example expression dataframe with 26 samples and 3 genes

# prepare weights dataframe
gene.names <- c("gene.1", "gene.2", "gene.3")
weights <- as.data.frame(gene.names)
weights$weight.value <- 1 #In this example, these genes all have a weight of 1

# prepare dataframe containing drug sensitivity metric AUC for each drug
drug.sensitivity <- data.frame(id)
Drug <- LETTERS[seq(from = 1, to = 26)]
for (i in 1:length(Drug)){
drug.sensitivity[,c(Drug[i])] <- seq(0,1, by = 0.04)
}

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

# perform DMEA
DMEA(drug.sensitivity, gmt, expression, weights)
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
