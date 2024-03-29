DMEA <- function(drug.sensitivity, gmt = NULL, expression, weights,
                 value = "AUC", sample.names = colnames(expression)[1],
                 gene.names = colnames(weights)[1],
                 weight.values = colnames(weights)[2],
                 rank.metric = "Pearson.est", FDR = 0.25,
                 num.permutations = 1000, stat.type = "Weighted",
                 drug.info = NULL, drug = "Drug", set.type = "moa",
                 min.per.set = 6, sep = "[|]",
                 exclusions = c("-666", "NA", "na", "NaN", "NULL"),
                 descriptions = NULL, min.per.corr = 3, scatter.plots = TRUE,
                 scatter.plot.type = "pearson", FDR.scatter.plots = 0.05,
                 xlab = "Weighted Voting Score", ylab = value,
                 position.x = "min", position.y = "min", se = TRUE) {
  # check classes of expression and weight data frames are as expected
  testthat::expect_is(expression, "data.frame")
  testthat::expect_is(weights, "data.frame")

  # Weighted Voting (WV)
  WV.result <- WV(
    expression = expression, weights = weights,
    sample.names = sample.names, gene.names = gene.names,
    weight.values = weight.values
  )

  # check classes of drug sensitivity data and WV output are as expected
  testthat::expect_is(drug.sensitivity, "data.frame")
  testthat::expect_is(WV.result$scores, "data.frame")

  # merge drug sensitivity dataframe with WV
  WV.result.drug.sensitivity <- merge(WV.result$scores, drug.sensitivity,
    by = sample.names
  )

  # rank_corr
  corr.results <- rank_corr(
    data = WV.result.drug.sensitivity, variable = drug, value = value,
    type = scatter.plot.type, min.per.corr = min.per.corr,
    plots = scatter.plots, FDR = FDR.scatter.plots, xlab = xlab, ylab = ylab,
    position.x = position.x, position.y = position.y, se = se
  )

  # check class of rank_corr output is as expected
  testthat::expect_is(corr.results$result, "data.frame")

  # if necessary, generate gmt object
  if (is.null(gmt)) {
    if (is.null(drug.info)) {
      stop(paste(
        "drug.info dataframe containing set membership information",
        "must be provided as input if no gmt object is provided"
      ))
    } else {
      corr.output <- merge(corr.results$result, drug.info, by = drug)
      gmt <- as_gmt(
        data = corr.output, element.names = drug, set.names = set.type,
        min.per.set, sep, exclusions, descriptions
      )
    }
  } else {
    corr.output <- corr.results$result
  }

  # Drug Mechanism Enrichment Analysis (DMEA)
  DMEA.results <- drugSEA(
    corr.output, gmt, drug, rank.metric, set.type,
    direction.adjust = NULL, FDR, num.permutations,
    stat.type, min.per.set, sep, exclusions, 
    descriptions, convert.synonyms = FALSE
  )

  return(list(
    WV.scores = WV.result$scores,
    unused.weights = WV.result$unused.weights,
    corr.result = corr.results$result,
    corr.scatter.plots = corr.results$scatter.plots,
    gmt = gmt,
    result = DMEA.results$result,
    mtn.plots = DMEA.results$mtn.plots,
    volcano.plot = DMEA.results$volcano.plot,
    removed.sets = DMEA.results$removed.sets,
    unannotated.drugs = DMEA.results$unannotated.drugs
  ))
}
