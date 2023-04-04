as_gmt <- function(data, element.names = "Drug", set.names = "moa", min.per.set=6,
                   sep = "[|]", exclusions = c("-666", "NA", "na", "NaN", "NULL"), descriptions = NULL){
  print("Generating gmt object for enrichment analysis...")

  all.sets <- unique(data[,c(set.names)])
  if(length(all.sets) > 0){
    all.sets <- all.sets[all.sets %in% exclusions == FALSE]
    if(length(all.sets) > 0){
      elements <- list()

      # use parallel computing if possible
      if(require(parallel) &
         require(snow) &
         require(doSNOW)){
        cores <- parallel::detectCores() # number of cores available
        if(cores[1] > 1){
          cl <- snow::makeCluster(cores[1]-1) # cluster using all but 1 core
          doSNOW::registerDoSNOW(cl) # register cluster
        }
      }

      # fill sets with elements
      final.sets <- c()
      final.elements <- list()
      for(i in 1:length(all.sets)){
        elements[[all.sets[i]]] <- list()
        for(j in 1:nrow(data)){
          sets <- strsplit(as.character(data[j, c(set.names)]), sep)[[1]]
          if(all.sets[i] %in% sets){elements[[all.sets[i]]] <- unique(c(elements[[all.sets[i]]], data[j, c(element.names)]))}
        }
        # store sets with elements
        if(length(elements[[all.sets[i]]])>=min.per.set){
          final.sets <- c(final.sets, all.sets[i])
          final.elements[[all.sets[i]]] <- elements[[all.sets[i]]]
        }
      }
      rm(elements)
      rm(all.sets)

      if(require(parallel) &
         require(snow) &
         require(doSNOW)){
        if(cores[1] > 1){
          snow::stopCluster(cl) #stop cluster
          rm(cl)
        }
      }

      # if there are no set descriptions, use set names
      if(!is.null(descriptions)){
        final.set.info <- dplyr::distinct(data[data[,c(set.names)] %in% final.sets, c(set.names, descriptions)])
        final.set.info <- final.set.info[, c(descriptions)]
      }else{
        final.set.info <- final.sets
      }

      gmt <- list(genesets=final.elements, geneset.names=final.sets, geneset.descriptions=final.set.info)
      return(gmt)
    }else{
      stop("no set annotations were left after exclusions were removed")
    }
  }else{
    stop("set annotations must be provided to generate a gmt object")
  }
}
