as_gmt <- function(data, element.names = "Drug", set.names = "moa",
                   min.per.set = 6, sep = "[|]",
                   exclusions = c("-666", "NA", "na", "NaN", "NULL"),
                   descriptions = NULL) {
  message("Generating gmt object for enrichment analysis...")

  # check class of input is as expected
  testthat::expect_is(data, "data.frame")

  all.sets <- unique(data[, c(set.names)])
  if (length(all.sets) > 0) {
    all.sets <- all.sets[!(all.sets %in% exclusions)]
    if (length(all.sets) > 0) {
      elements <- list()

      # use parallel computing if possible
      if (requireNamespace("parallel") &
        requireNamespace("snow") &
        requireNamespace("doSNOW")) {
        cores <- parallel::detectCores() # number of cores available
        if (cores[1] > 1) {
          cl <- snow::makeCluster(cores[1] - 1) # cluster using all but 1 core
          doSNOW::registerDoSNOW(cl) # register cluster
        }
      }

      # fill sets with elements
      final.sets <- c()
      final.elements <- list()
      for (i in seq_len(length(all.sets))) {
        elements[[all.sets[i]]] <- list()
        for (j in seq_len(nrow(data))) {
          sets <- strsplit(as.character(data[j, c(set.names)]), sep)[[1]]
          if (all.sets[i] %in% sets) {
            elements[[all.sets[i]]] <- unique(c(
              elements[[all.sets[i]]],
              data[j, c(element.names)]
            ))
          }
        }
        # store sets with elements
        if (length(elements[[all.sets[i]]]) >= min.per.set) {
          final.sets <- c(final.sets, all.sets[i])
          final.elements[[all.sets[i]]] <- elements[[all.sets[i]]]
        }
      }
      rm(elements)
      rm(all.sets)

      if (requireNamespace("parallel") &
        requireNamespace("snow") &
        requireNamespace("doSNOW")) {
        if (cores[1] > 1) {
          snow::stopCluster(cl) # stop cluster
          rm(cl)
        }
      }

      # if there are no set descriptions, use set names
      if (!is.null(descriptions)) {
        final.set.info <- dplyr::distinct(data[
          data[, c(set.names)] %in%
            final.sets,
          c(set.names, descriptions)
        ])
        final.set.info <- final.set.info[, c(descriptions)]
      } else {
        final.set.info <- final.sets
      }

      # compile gmt information
      gmt <- list(
        genesets = final.elements,
        geneset.names = final.sets,
        geneset.descriptions = final.set.info
      )

      # check class of output is as expected
      testthat::expect_is(gmt, "list")

      return(gmt)
    } else {
      stop("no set annotations were left after exclusions were removed")
    }
  } else {
    stop("set annotations must be provided to generate a gmt object")
  }
}
