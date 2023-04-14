WV <- function(expression, weights, sample.names = colnames(expression)[1],
               gene.names = colnames(weights)[1],
               weight.values = colnames(weights)[2]) {
  # expression: cell lines should be rows, genes should be columns
  # weights: one column should contain gene symbols,
  # another column should be gene weights
  message("Calculating Weighted Voting scores...")

  # check classes of inputs are as expected
  testthat::expect_is(expression, "data.frame")
  testthat::expect_is(weights, "data.frame")

  # remove gene names not found in expression
  unused.weights <- weights[!(weights[, c(gene.names)] %in%
                                colnames(expression)[2:ncol(expression)]), ]
  filtered.weights <- weights[weights[, c(gene.names)] %in%
                                colnames(expression)[2:ncol(expression)],
                              c(gene.names, weight.values)]
  filtered.weights <- na.omit(filtered.weights)

  # check expression for gene names
  filtered.expr <- expression %>%
    dplyr::select(c(tidyselect::all_of(sample.names),
                    tidyselect::all_of(filtered.weights[, c(gene.names)])))
  filtered.expr <- na.omit(filtered.expr)

  if (nrow(filtered.expr) > 0) {
    ## prep for matrix multiplication
    # store sample names in row names
    rownames(filtered.expr) <- filtered.expr[, c(sample.names)]

    # remove sample names from expression
    filtered.expr.data <- dplyr::select(filtered.expr,
                                        -c(tidyselect::all_of(sample.names)))

    # remove gene names from weights
    weight.data <- dplyr::select(filtered.weights,
                                 -c(tidyselect::all_of(gene.names)))

    # format data frames as matrices
    weight.matrix <- as.matrix(weight.data)
    expr.matrix <- as.matrix(filtered.expr.data)

    ## perform matrix multiplication
    scores <- as.data.frame(expr.matrix %*% weight.matrix)

    ## format dataframe for output
    colnames(scores)[1] <- "WV"
    scores[, c(sample.names)] <- rownames(scores)
    scores <- scores %>% relocate(colnames(scores)[2], .before = WV)

    # check classes of outputs are as expected
    testthat::expect_is(scores, "data.frame")
    testthat::expect_is(unused.weights, "data.frame")
  } else {
    stop(paste("expression data frame does not contain",
               "any of the gene symbols provided"))
  }
  return(list(scores = scores, unused.weights = unused.weights))
}
