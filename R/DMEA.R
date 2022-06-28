#DMEA
#Author: Belinda B. Garana (BG), Date: 2020-12-03; last edit: BG 2022-06-28
#Note: drugSEA co-authored with JJ (GSEA_custom by JJ & revised by BG for drugSEA; gsea_mountain_plot by JJ & revised by BG)
#Note: thanks to NG for ng.theme (used in rank.corr)

rm(list=ls(all=TRUE))

WV <- function(expression, weights, sample.names=colnames(expression)[1],
                gene.names=colnames(weights)[1], weight.values=colnames(weights)[2]){
  # data: cells should be rows, genes should be columns; weights: one column should be gene names, another column should be gene weights
  print("Calculating Weighted Voting scores...")
  library(dplyr);library(tidyselect);
  
  # check expression for gene names
  filtered.expr <- expression %>% dplyr::select(c(tidyselect::all_of(sample.names),tidyselect::all_of(weights[,c(gene.names)])))
  filtered.expr <- na.omit(filtered.expr)
  if(nrow(filtered.expr)>0){
    # prep for matrix multiplication
    rownames(filtered.expr) <- filtered.expr[,c(sample.names)] #store sample names in row names
    filtered.expr.data <- dplyr::select(filtered.expr, -c(tidyselect::all_of(sample.names))) #remove sample names from expression
    filtered.weights <- weights[weights[,c(gene.names)] %in% colnames(filtered.expr.data),c(gene.names, weight.values)] #reduce weights for overlap with expression
    weight.data <- dplyr::select(filtered.weights, -c(tidyselect::all_of(gene.names))) #remove gene names from weights
    weight.matrix <- as.matrix(weight.data)
    expr.matrix <- as.matrix(filtered.expr.data)
    
    # perform matrix multiplication
    scores <- as.data.frame(expr.matrix %*% weight.matrix)
    
    # format dataframe for output
    colnames(scores)[1] <- "WV"
    scores[,c(sample.names)] <- rownames(scores)
    scores <- scores %>% relocate(colnames(scores)[2], .before=WV)
   }else{
    stop("expression dataframe does not contain data for the gene names provided")
  }
  return(scores)
}

rank.corr <- function(data, variable="Drug", value="AUC",type="pearson", min.per.corr=3, plots=TRUE, FDR=0.05, xlab=colnames(data)[2], ylab=value, position.x="mid", position.y="max", se=TRUE){
  print("Running correlations and regressions...")
  library(dplyr);library(qvalue);library(ggplot2);library(stats);library(reshape2);library(gridExtra);
  
  cores <- parallel::detectCores() # number of cores available
  if(cores[1] > 1){
    cl <- snow::makeCluster(cores[1]-1) # cluster using all but 1 core
    doSNOW::registerDoSNOW(cl) # register cluster
  }
  
  # correlations and regressions
  not_all_na <- function(x) any(!is.na(x))
  all.data.corr <- data
  all.data.corr <- all.data.corr %>% select_if(not_all_na)
  variable.set <- colnames(all.data.corr[,3:ncol(all.data.corr)])
  corr <- data.frame(variable.set)
  colnames(corr)[1] <- variable
  corr[,c("Pearson.est","Pearson.p","Pearson.q",
          "Spearman.est","Spearman.p","Spearman.q",
          "Slope","Intercept","R.squared",
          "Rank.slope","Rank.intercept","Rank.R.squared","N")] <- NA
  data.corr <- list()
  for (j in 3:ncol(all.data.corr)){
    data.corr[[j-2]] <- na.omit(all.data.corr[,c(2,j)])
    a <- nrow(data.corr[[j-2]]) # number of samples in the correlation (N)
    if (a >= min.per.corr) {corr$N[j-2] <- a
    x <- as.numeric(data.corr[[j-2]][,1]) # list of rank metric
    y <- as.numeric(data.corr[[j-2]][,2]) # list of gene expression (for each gene j)

    Regression <- stats::lm(y ~ x)
    corr$Slope[j-2] <- Regression$coeff[[2]]
    corr$Intercept[j-2] <- Regression$coeff[[1]]
    corr$R.squared[j-2] <- summary(Regression)$r.squared

    Rank.regression <- stats::lm(rank(y) ~ rank(x))
    corr$Rank.slope[j-2] <- Rank.regression$coeff[[2]]
    corr$Rank.intercept[j-2] <- Rank.regression$coeff[[1]]
    corr$Rank.R.squared[j-2] <- summary(Rank.regression)$r.squared

    Pearson <- stats::cor.test(x, y, method="pearson")
    corr$Pearson.est[j-2] <- Pearson$estimate
    corr$Pearson.p[j-2] <- Pearson$p.value

    Spearman <- stats::cor.test(x, y, method="spearman", exact=FALSE)
    corr$Spearman.est[j-2] <- Spearman$estimate
    corr$Spearman.p[j-2] <- Spearman$p.value}
  }
  corr.no.na <- corr[!is.na(corr$N), ]
  qobj.pearson <- qvalue(p=corr.no.na$Pearson.p,pi0=1)
  qobj.spearman <- qvalue(p=corr.no.na$Spearman.p,pi0=1)
  corr.no.na$Pearson.q <- qobj.pearson$qvalues
  corr.no.na$Spearman.q <- qobj.spearman$qvalues

  # scatter plots for significant results
  if (plots==TRUE){
    # load themes for plots
    ng.theme <- ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.border = element_rect(fill=NA), panel.background = element_blank(),
                      axis.line = element_line(colour = "black"), axis.text.x = element_text(colour = "black"),
                      axis.text.y = element_text(colour = "black"), axis.ticks.x = element_line(colour="black"),
                      axis.ticks.y = element_line(colour="black"), legend.title = element_blank(),
                      axis.title.y = element_text(size=8, colour="black"))

    bg.theme <- ggplot2::theme(legend.background = element_rect(), legend.position="top", legend.text = element_text(size=14), legend.key = element_blank(),
                      axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                      axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                      plot.title = element_text(lineheight=.8, face="bold", size=36))
    # prepare dataframe for plots
    id.var <- colnames(data)[1]
    rank.var <- colnames(data)[2]
    plot.data <- reshape2::melt(all.data.corr, id=c(id.var,rank.var), variable.name = variable, value.name = value, na.rm=TRUE)
    if (type=="pearson"){
      results <- unique(corr.no.na[which(corr.no.na$Pearson.q<=FDR),1])
      if (length(results)>0){
        a <- list()
        for (i in 1:length(results)) {
          sig.data <- plot.data[plot.data[,3]==results[i],]
          min.x <- min(sig.data[,c(rank.var)])
          max.x <- max(sig.data[,c(rank.var)])
          mid.x <- 0.5*(min.x+max.x)
          min.y <- min(sig.data[,c(value)])
          max.y <- max(sig.data[,c(value)])
          mid.y <- 0.5*(min.y+max.y)
          if (position.x=="min"){pos.x=min.x} else if(position.x=="mid"){pos.x=mid.x
          } else if(position.x=="max"){pos.x=max.x} else if(is.numeric(position.x)==TRUE){pos.x==position.x}
          if (position.y=="min"){pos.y=min.y} else if(position.y=="mid"){pos.y=mid.y
          } else if(position.y=="max"){pos.y=max.y} else if(is.numeric(position.y)==TRUE){pos.y==position.y}

          Pearson.est <- corr.no.na[corr.no.na[,1]==results[i],]$Pearson.est
          Pearson.p <- corr.no.na[corr.no.na[,1]==results[i],]$Pearson.p
          stats_pearson <- substitute(r == est*","~~"p"~"="~p,
                                      list(est = format(Pearson.est, digits = 3),
                                           p = format(Pearson.p, digits = 3)))

          a[[i]] <- ggplot2::ggplot(data=sig.data, aes_string(x=rank.var, y=value)) +
            ggplot2::geom_point() + ggplot2::labs(x = xlab, y = ylab) + ggplot2::ggtitle(results[i]) + 
            ggplot2::geom_smooth(method="lm",size = 1.5,linetype = 'solid',color="blue",se = se,na.rm = TRUE) +
            ggplot2::geom_text(x=pos.x,y=pos.y,vjust="inward",hjust="inward", colour="blue", parse=TRUE,
                       label = as.character(as.expression(stats_pearson)), size = 8) + ng.theme + bg.theme
        }
        scatter.plots <- gridExtra::marrangeGrob(a, nrow=1, ncol=1)
      }else{
        scatter.plots <- NA
        warning("No correlations met the FDR cut-off to produce scatter plots")
      }
    } else if (type=="spearman"){
      results <- unique(corr.no.na[which(corr.no.na$Spearman.q<=FDR),1])
      if (length(results)>0){
        a <- list()
        for (i in 1:length(results)) {
          sig.data <- plot.data[plot.data[,3]==results[i],]
          min.x <- 1
          mid.x <- 0.5*length(sig.data[,c(rank.var)])
          max.x <- length(sig.data[,c(rank.var)])
          min.y <- 1
          mid.y <- 0.5*length(sig.data[,c(value)])
          max.y <- length(sig.data[,c(value)])
          if (position.x=="min"){pos.x=min.x} else if(position.x=="mid"){pos.x=mid.x
          } else if(position.x=="max"){pos.x=max.x} else if(is.numeric(position.x)==TRUE){pos.x==position.x}
          if (position.y=="min"){pos.y=min.y} else if(position.y=="mid"){pos.y=mid.y
          } else if(position.y=="max"){pos.y=max.y} else if(is.numeric(position.y)==TRUE){pos.y==position.y}

          Spearman.est <- corr.no.na[corr.no.na[,1]==results[i],]$Spearman.est
          Spearman.p <- corr.no.na[corr.no.na[,1]==results[i],]$Spearman.p
          stats_spearman <- substitute(rho == est*","~~"p"~"="~p,
                                       list(est = format(Spearman.est, digits = 3),
                                            p = format(Spearman.p, digits = 3)))

          a[[i]] <- ggplot(data=sig.data, aes_string(x=rank(sig.data[,c(rank.var)]), y=rank(sig.data[,c(value)]))) +
            geom_point() + labs(x = paste0(xlab," Rank"), y = paste0(ylab," Rank")) + ggtitle(results[i]) +
            geom_smooth(method="lm",size = 1.5,linetype = 'solid',color="blue",se = se,na.rm = TRUE) +
            geom_text(x = pos.x, y = pos.y, label = as.character(as.expression(stats_spearman)), colour="blue",
                      size = 8, parse=TRUE) + ng.theme + bg.theme
        }
        scatter.plots <- gridExtra::marrangeGrob(a, nrow=1, ncol=1)
      }else{
        warning("No correlations met the FDR cut-off to produce scatter plots")
        scatter.plots <- NA
      }
    }else{
      warning("Type must be specified as either spearman or pearson to produce scatter plots")
      scatter.plots <- NA
    }
  }else{scatter.plots <- NA}
  
  if(cores[1] > 1){
    snow::stopCluster(cl) # stop cluster
    rm(cl)
  }
  
  outputs <- list(result=corr.no.na, scatter.plots = scatter.plots)
  return(outputs)
}

as.gmt <- function(data, element.names = "Drug", set.names = "moa", min.per.set=6, 
                   sep = "[|]", exclusions = c("-666", "NA", "na", "NaN", "NULL"), descriptions = NULL){
  print("Generating gmt object for enrichment analysis...")
  library(dplyr);library(parallel);library(snow);library(doSNOW);
  
  all.sets <- unique(data[,c(set.names)])
  if(length(all.sets) > 0){
    all.sets <- all.sets[all.sets %in% exclusions == FALSE]
    if(length(all.sets) > 0){
      elements <- list()
      
      cores <- parallel::detectCores() # number of cores available
      if(cores[1] > 1){
        cl <- snow::makeCluster(cores[1]-1) # cluster using all but 1 core
        doSNOW::registerDoSNOW(cl) # register cluster
      }
      
      final.sets <- c()
      final.elements <- list()
      for(i in 1:length(all.sets)){
        elements[[all.sets[i]]] <- list()
        for(j in 1:nrow(data)){ 
          sets <- strsplit(as.character(data[j, c(set.names)]), sep)[[1]]
          if(all.sets[i] %in% sets){elements[[all.sets[i]]] <- c(elements[[all.sets[i]]], data[j, c(element.names)])}
        }
        # store sets with elements
        if(length(elements[[all.sets[i]]])>=min.per.set){
          final.sets <- c(final.sets, all.sets[i])
          final.elements[[all.sets[i]]] <- elements[[all.sets[i]]]
        }
      }
      rm(elements)
      rm(all.sets)
      if(cores[1] > 1){
        snow::stopCluster(cl) #stop cluster
        rm(cl)
      }
      
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

drugSEA <- function(data, gmt=NULL, drug="Drug", rank.metric="Pearson.est", set.type="moa", direction.adjust=NULL,
                    FDR=0.25, num.permutations=1000, stat.type="Weighted", min.per.set=6,
                    sep = "[|]", exclusions = c("-666", "NA", "na", "NaN", "NULL"), descriptions=NULL){
  library(dplyr);library(parallel);library(snow);library(doSNOW);library(stats);
  library(sjmisc);library(testthat);library(foreach);
  library(ggplot2);library(cowplot);library(aplot);library(ggrepel);
  
  # load GSEA_custom function
  GSEA_custom <- function(input.df, gmt.list,
                          num.permutations = 1000,
                           stat.type = "Weighted", min.per.set){
    nperm = num.permutations # number of permutations
    if(stat.type == "Classic"){
      score.weight = 0
    }else if(stat.type == "Weighted"){
      score.weight = 1
    }

    # Read in gene expression data
    # drug names should be first column
    # rank metric name should be second column
    data_in <- input.df

    # find biggest drug set
    biggest.set <- min.per.set
    for(i in 1:length(gmt.list$genesets)){
      if(length(gmt.list$genesets[[i]]) > biggest.set){
        biggest.set <- length(gmt.list$genesets[[i]])
      }
    }
    
    Drug.Sets <- matrix(data = "", nrow = biggest.set, ncol = length(gmt.list$genesets))
    colnames(Drug.Sets) <- gmt.list$geneset.names
    for(i in 1:ncol(Drug.Sets)){
      for(j in 1:length(gmt.list$genesets[[i]])){
        Drug.Sets[j,i] <- gmt.list$genesets[[i]][[j]][1]
      }
    }
    Drug.Sets <- as.data.frame(Drug.Sets)

    testthat::expect_is(data_in, "data.frame")
    testthat::expect_is(Drug.Sets, "data.frame")

    GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = score.weight, correl.vector = NULL){
      tag.indicators <- sign(match(gene.list, gene.set, nomatch = 0))
      no.tag.indicator <- 1 - tag.indicators
      N <- length(gene.list)
      Nh <- numhits_pathway
      Nm <- N - Nh
      if (weighted.score.type == 0){
        correl.vector <- rep(1,N)
      }
      alpha <- weighted.score.type
      correl.vector <- abs(correl.vector**alpha)
      sum.correl.tag <- sum(correl.vector[tag.indicators == 1])
      norm.tag <- 1.0/sum.correl.tag
      norm.no.tag <- 1.0/Nm
      RES <- cumsum(tag.indicators * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
      max.ES <- max(RES)
      min.ES <- min(RES)
      if (max.ES > - min.ES) {
        ES <- signif(max.ES, digits = 5)
        arg.ES <- which.max(RES)
      } else {
        ES <- signif(min.ES, digits=5)
        arg.ES <- which.min(RES)
      }
      return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicators))
    } # for real ES

    GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = score.weight, correl.vector = NULL) {

      N <- length(gene.list)
      Nh <- numhits_pathway
      Nm <-  N - Nh

      loc.vector <- vector(length=N, mode="numeric")
      peak.res.vector <- vector(length=Nh, mode="numeric")
      valley.res.vector <- vector(length=Nh, mode="numeric")
      tag.correl.vector <- vector(length=Nh, mode="numeric")
      tag.diff.vector <- vector(length=Nh, mode="numeric")
      tag.loc.vector <- vector(length=Nh, mode="numeric")

      loc.vector[gene.list] <- seq(1, N)
      tag.loc.vector <- loc.vector[gene.set]

      tag.loc.vector <- sort(tag.loc.vector, decreasing = F)

      if(weighted.score.type == 0){
        tag.correl.vector <- rep(1, Nh)
      }else if(weighted.score.type == 1){
         tag.correl.vector <- correl.vector[tag.loc.vector]
         tag.correl.vector <- abs(tag.correl.vector)
      }else if(weighted.score.type == 2){
         tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
         tag.correl.vector <- abs(tag.correl.vector)
       }else{
         tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
         tag.correl.vector <- abs(tag.correl.vector)
      }

      norm.tag <- 1.0/sum(tag.correl.vector)
      tag.correl.vector <- tag.correl.vector * norm.tag
      norm.no.tag <- 1.0/Nm
      tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
      tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
      tag.diff.vector <- tag.diff.vector * norm.no.tag
      peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
      valley.res.vector <- peak.res.vector - tag.correl.vector
      max.ES <- max(peak.res.vector)
      min.ES <- min(valley.res.vector)
      ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
      return(ES)
    } #for permutation ES

    Samples <- colnames(data_in)[2]
    Drug.Sets.All <- colnames(Drug.Sets)

    annotations <- matrix(data = 0, nrow = nrow(data_in), ncol = length(Drug.Sets.All))
    colnames(annotations) <- Drug.Sets.All

    annotations <- as.data.frame(annotations)

    annotations <- cbind(data_in[,1],annotations)
    colnames(annotations) <- c(colnames(data_in)[1], Drug.Sets.All)
    annotations <- as.matrix(annotations)

    num.hits.pathways <- list()

    ### Annotate drug sets
    for (j in 1:length(Drug.Sets.All)){
      temp.pathway <- Drug.Sets[,Drug.Sets.All[j]]
      for (i in 1:nrow(annotations)){
        if (annotations[i,1] %in% temp.pathway){
          annotations[i,j+1] = "X";
        }
      }
      num.hits.pathways[[Drug.Sets.All[j]]] <- sum(annotations[,Drug.Sets.All[j]] == "X")
    }

    num.hits.pathways.df <- matrix(unlist(num.hits.pathways))
    row.names(num.hits.pathways.df) = Drug.Sets.All
    Drug.Sets.under <- which(num.hits.pathways.df < min.per.set)
    Drug.Sets.over <- which(num.hits.pathways.df >= min.per.set)
    Drug.Sets.All <- Drug.Sets.All[Drug.Sets.over] # only keep drug sets with >= min.per.set drugs
    if(length(Drug.Sets.under) > 0){
      warning(paste0("Removing drug sets with less than ",min.per.set," drugs observed in data set"))
      annotations <- annotations[,c(colnames(data_in)[1],Drug.Sets.All)]
    }
    if(length(Drug.Sets.over) < 2){stop("annotations for 2+ drug sets are required")} # BG 2022-06-25
    annotations <- as.data.frame(annotations)
    data_in <- merge(data_in, annotations, by = colnames(data_in)[1])
    data_in <- stats::na.omit(data_in)
    rm(annotations)
    
    # Find out how many cores are available (if you don't already know)
    cores<-parallel::detectCores()
    if(cores[1] > 1){
      cl <- snow::makeCluster(cores[1]-1) # cluster using all but 1 core
      doSNOW::registerDoSNOW(cl) # register cluster
    }

    GSEA.Results <- matrix(data = NA, nrow = length(Drug.Sets.All), ncol = 7)
    colnames(GSEA.Results) <- c("Rank_metric","Drug_set","ES","NES",
                                "p_value","Position_at_max","FDR_q_value")
    GSEA.Results <- as.data.frame(GSEA.Results)
    GSEA.Results$Drug_set <- Drug.Sets.All
    GSEA.Results$Rank_metric <- Samples

    ## Assuming first two columns in data table are drug names and rank metric (e.g. Foldchange, SNR)
    data_in[,Samples] <- as.numeric(as.character(data_in[,Samples]))
    data_in <- data_in[order(-data_in[,Samples]),] # sort by descending order for the rank metric
    rownames(data_in) <- 1:nrow(data_in) # reorder row indices for counting in for loop below
    ions <- nrow(data_in)
    
    #for plotting
    ks_results_plot <- list()
    positions.of.hits <- list()
    
    gene.list <- 1:ions
    rank_metric <- data_in[,Samples] #Save the rank metric
    
    pos_gene_set <- array(data = 0, dim = nrow(data_in), dimnames = NULL);
    
    ## Calculate Real KS Statistic
    for (i in 1:length(Drug.Sets.All)){
      data_in3 <- data_in[,Drug.Sets.All[i]]
      numhits_pathway <- sum(data_in3 == "X"); # check to see if there is anything in the column (e.g. X)
      if (numhits_pathway > 1){
        pos_gene_set <- which(data_in[,Drug.Sets.All[i]] %in% c("X"))
        KS_real <- GSEA.EnrichmentScore(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
        GSEA.Results[GSEA.Results$Drug_set == Drug.Sets.All[i],]$ES <- KS_real$ES;
        GSEA.Results[GSEA.Results$Drug_set == Drug.Sets.All[i],]$Position_at_max <- KS_real$arg.ES;
        ks_results_plot[[Drug.Sets.All[i]]] = KS_real$RES
        positions.of.hits[[Drug.Sets.All[i]]] = pos_gene_set
      }
    }
    Mountain.Plot.Info <- list(MountainPlot = ks_results_plot, Position.of.hits = positions.of.hits)
    rm(pos_gene_set)
    rm(numhits_pathway)
    rm(data_in3)
    rm(KS_real)
    
    KSRandomArray <- matrix(data = NA, nrow = nperm, ncol = length(Drug.Sets.All))
    num.Drug.Sets.all <- length(Drug.Sets.All)
    `%dopar%` <- foreach::`%dopar%`
    
    KSRandomArray <- foreach::foreach(L = 1:nperm, .combine = "rbind") %dopar% {
      temp.KSRandomArray <- matrix(data = NA, nrow = 1, ncol = num.Drug.Sets.all)
      for(i in 1:length(Drug.Sets.All)){
        numhits_pathway <- length(positions.of.hits[[Drug.Sets.All[i]]])
        pos_gene_set <- sample(1:ions,numhits_pathway)
        temp.KSRandomArray[,i] <- GSEA.EnrichmentScore2(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
      }
      temp.KSRandomArray
    }
    colnames(KSRandomArray) <- Drug.Sets.All
    
    KSRandomArray <- data.frame(matrix(unlist(KSRandomArray), nrow = nperm, byrow = T))
    colnames(KSRandomArray) <- Drug.Sets.All
    KSRandomArray <- stats::na.omit(KSRandomArray)
    
    KSRandomArray <- as.data.frame(KSRandomArray)
    ### normalize the GSEA distribution
    KSRandomArray.Norm <- matrix(data = NA, nrow = nrow(KSRandomArray), ncol = ncol(KSRandomArray))
    colnames(KSRandomArray.Norm) <- colnames(KSRandomArray)
    avg <- 0
    KSRandomArray.temp <- 0
    for (i in 1:ncol(KSRandomArray.Norm)){
      avg <- 0
      KSRandomArray.temp <- KSRandomArray[,i]
      pos.temp <- KSRandomArray.temp[which(KSRandomArray.temp >= 0)]
      neg.temp <- KSRandomArray.temp[which(KSRandomArray.temp <= 0)] #BG OR EQUAL TO
      
      avg.pos <- mean(pos.temp)
      avg.neg <- mean(neg.temp)
      
      norm.pos.temp <- pos.temp / avg.pos
      norm.neg.temp <- neg.temp / avg.neg * -1
      
      norm.perms <- c(norm.pos.temp,norm.neg.temp)
      
      KSRandomArray.Norm[,i] <- norm.perms
      
    }
    GSEA.NES.perms <- as.vector(KSRandomArray.Norm)
    rm(KSRandomArray.Norm)
    GSEA.NES.perms.pos <- GSEA.NES.perms[which(GSEA.NES.perms >= 0)]
    GSEA.NES.perms.neg <- GSEA.NES.perms[which(GSEA.NES.perms <= 0)] #BG OR EQUAL TO
    rm(GSEA.NES.perms)
    percent.pos.GSEA <- sum(GSEA.Results$ES >= 0) / length(GSEA.Results$ES) #BG OR EQUAL TO
    percent.neg.GSEA <- sum(GSEA.Results$ES <= 0) / length(GSEA.Results$ES) #BG OR EQUAL TO
    
    # Calculate GSEA NES and p-value
    for (i in 1:length(Drug.Sets.All)){
      temp.gene.set <- Drug.Sets.All[i]
      temp.ES <- GSEA.Results[GSEA.Results$Drug_set == temp.gene.set,]$ES
      if (temp.ES >= 0){ # BG OR EQUAL TO
        pos.perms <- KSRandomArray[,temp.gene.set]
        pos.perms <- pos.perms[which(pos.perms >= 0)] # BG OR EQUAL TO
        #p-val
        GSEA.Results[GSEA.Results$Drug_set == temp.gene.set,]$p_value = signif(sum(pos.perms >= temp.ES) / length(pos.perms),digits = 3) #BG OR EQUAL TO
        #NES
        GSEA.Results[GSEA.Results$Drug_set == temp.gene.set,]$NES = signif(temp.ES / mean(pos.perms), digits = 3)
      } else if (temp.ES <= 0){ # BG OR EQUAL TO
        neg.perms <- KSRandomArray[,temp.gene.set]
        neg.perms <- neg.perms[which(neg.perms <= 0)] # BG OR EQUAL TO
        #p-val
        GSEA.Results[GSEA.Results$Drug_set == temp.gene.set,]$p_value = signif(sum(neg.perms <= temp.ES) / length(neg.perms),digits = 3) #BG OR EQUAL TO
        #NES
        GSEA.Results[GSEA.Results$Drug_set == temp.gene.set,]$NES = signif(temp.ES / mean(neg.perms) * -1, digits = 3)
      }
    }
    
    # Calculate GSEA FDR
    for(i in 1:length(Drug.Sets.All)){
      temp.gene.set <- Drug.Sets.All[i]
      temp.NES <- GSEA.Results[GSEA.Results$Drug_set == temp.gene.set,]$NES
      if(temp.NES >= 0){ # BG OR EQUAL TO
        percent.temp <- sum(GSEA.NES.perms.pos >= GSEA.Results[GSEA.Results$Drug_set == temp.gene.set,]$NES) / length(GSEA.NES.perms.pos) #BG OR EQUAL TO
        percent.pos.stronger <- sum(GSEA.Results$NES >= GSEA.Results[GSEA.Results$Drug_set == temp.gene.set,]$NES) / sum(GSEA.Results$NES >= 0) #BG
        GSEA.Results[GSEA.Results$Drug_set == temp.gene.set,]$FDR_q_value = ifelse(signif(percent.temp / percent.pos.stronger, digits = 3) < 1, signif(percent.temp / percent.pos.stronger, digits = 3), 1) #BG
      }else if(temp.NES <= 0){ # BG OR EQUAL TO
        percent.temp <- sum(GSEA.NES.perms.neg <= temp.NES) / length(GSEA.NES.perms.neg) # BG OR EQUAL TO
        percent.neg.stronger <- sum(GSEA.Results$NES <= temp.NES) / sum(GSEA.Results$NES <= 0) # BG
        GSEA.Results[GSEA.Results$Drug_set == temp.gene.set,]$FDR_q_value = ifelse(signif(percent.temp / percent.neg.stronger, digits = 3) < 1, signif(percent.temp / percent.neg.stronger, digits = 3), 1) #BG
      }
    }
    if(cores[1] > 1){
      snow::stopCluster(cl) # stop cluster
      rm(cl)
    }

    return(list(GSEA.Results = GSEA.Results,
                Mountain.Plot.Info = Mountain.Plot.Info,
                ranking.metric = rank_metric))
  }
  
  # load mountain plot function
  gsea_mountain_plot <- function(GSEA.list, Sample.Name, Gene.Set.A, color = TRUE){
    library(ggplot2);
    if(color == TRUE){
      color.palette <- c("red")
    }else if(color == FALSE){
      color.palette <- c("black")
    }

    GSEA.list <- GSEA.list
    plotting.info <- GSEA.list$Mountain.Plot.Info
    ES_profiles <- plotting.info$MountainPlot
    position.of.hits <- plotting.info$Position.of.hits

    ranking.metric <- as.numeric(GSEA.list$ranking.metric)

    ES_GeneSetA <- ES_profiles[[Gene.Set.A]]

    pos.GeneSetA <- position.of.hits[[Gene.Set.A]]

    myName <- Gene.Set.A

    GSEA.results <- GSEA.list$GSEA.Results
    GSEA.results <- GSEA.results[GSEA.results$Rank_metric == Sample.Name,]
    GSEA.results$p_value <- as.numeric(as.character(GSEA.results$p_value))
    GSEA.results$FDR_q_value <- as.numeric(as.character(GSEA.results$FDR_q_value))

    gsea.normalized.enrichment.score <- signif(GSEA.results[GSEA.results$Drug_set == myName,]$NES, digits =3)
    gsea.p.value <- signif(GSEA.results[GSEA.results$Drug_set == myName,]$p_value, digits = 3)
    gsea.q.val <- signif(GSEA.results[GSEA.results$Drug_set == myName,]$FDR_q_value, digits = 3)

    results.as.text <- paste("NES",gsea.normalized.enrichment.score,"\np-value:", gsea.p.value,
                             "\nq-value:",gsea.q.val)

    gsea.hit.indices.A <- pos.GeneSetA

    gsea.es.profile.A <- ES_GeneSetA[gsea.hit.indices.A]

    ngenes <- length(ES_GeneSetA)
    combined.ES <- c(ES_GeneSetA)
    max.ES <- max(combined.ES)
    min.ES <- min(combined.ES)

    mtn.plot <- NULL
    mtn.plot <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x=c(0,pos.GeneSetA,ngenes), y = c(0,gsea.es.profile.A,0)), color = color.palette[1], size = 1) +
      ggplot2::labs(x = NULL, y = "ES", title = myName) + ggplot2::scale_x_continuous(expand = c(0,0)) +
      ggplot2::theme_bw() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.length = ggplot2::unit(0, "pt"),
        axis.title.x = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0,0,0,0),"cm"),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

    if(!is.null(gsea.normalized.enrichment.score) & length(gsea.normalized.enrichment.score) > 0){
      if(gsea.normalized.enrichment.score > 0){
        mtn.plot <- mtn.plot + ggplot2::annotate("text",label = results.as.text, x = ngenes * 0.90, y= max.ES * 0.70, color = "black")
      }else if(gsea.normalized.enrichment.score < 0){
        mtn.plot <- mtn.plot + ggplot2::annotate("text",label = results.as.text, x = ngenes * 0.10, y= min.ES * 0.30, color = "black")
      }
    }

    hit.markers <- ggplot2::ggplot() +
      cowplot::theme_nothing() +
      ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA)
      ) +
      ggplot2::geom_col(ggplot2::aes(x = pos.GeneSetA, y = 1), color = color.palette[1], fill = color.palette[1], width = 0.75) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      aplot::xlim2(mtn.plot)

    show.rank.metric <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = 1:ngenes, y = ranking.metric), color = "grey50") +
      ggplot2::geom_area(ggplot2::aes(x = 1:ngenes, y = ranking.metric), fill = "grey50") +
      ggplot2::labs(x = "Rank", y = "Rank Metric") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0,0,0,0),"cm")
      ) +
      aplot::xlim2(mtn.plot)
    
    assembled <- mtn.plot + hit.markers + show.rank.metric + patchwork::plot_layout(ncol = 1, heights = c(2.5,0.5,2))
    return(assembled)
  }
  
  ## remove any drugs with rank.metric = 0, NA (signed ranks are required for enrichment analysis)
  data[,c(rank.metric)] <- as.numeric(data[,c(rank.metric)])
  data <- data[data[,c(rank.metric)]!=0 & !is.na(data[,c(rank.metric)]),]
  
  ## check for duplicate entries for each drug (only 1 rank metric allowed per drug)
  drug.names <- unique(data[,c(drug)])
  if(length(drug.names) < nrow(data)){
    stop("each drug must have only 1 rank.metric value")
  }
  
  ## if necessary, generate gmt object
  if(is.null(gmt)){
    gmt <- as.gmt(data, element.names=drug, set.names=set.type, min.per.set, sep, exclusions, descriptions)
  }
  
  ## if necessary, adjust rank metrics
  if(length(direction.adjust)>0){
    info <- dplyr::distinct(data[,c(drug, rank.metric, set.type)])
    info$adjust <- 1
    info[grep(direction.adjust, info[,c(set.type)], value = T), "adjust"] <- -1
    info$adjusted.est <- as.numeric(info[,c(rank.metric)])*as.numeric(info$adjust)
    input <- info[,c(drug,"adjusted.est")]
    est.name <- paste0(rank.metric, " direction-adjusted")
  }else{
    input <- data[,c(drug, rank.metric)]
    est.name <- rank.metric
  }
  colnames(input) <- c(drug, est.name)
  
  ## perform enrichment analysis using GSEA_custom
  print("Running enrichment analysis...")
  EA <- GSEA_custom(input, gmt, num.permutations, stat.type, min.per.set)
  EA.Results <- EA$GSEA.Results
  
  ## select results which meet FDR threshold and produce mountain plots
  significant.hits <- EA.Results[which(EA.Results$FDR_q_value < FDR),]
  temp.plot <- list()
  if (nrow(significant.hits) > 0){
    for (i in 1:nrow(significant.hits)){
      temp <- gsea_mountain_plot(GSEA.list = EA, Sample.Name = est.name, Gene.Set.A = significant.hits$Drug_set[i])
      temp.plot[significant.hits$Drug_set[i]] <- list(temp)
    }
  }else{warning("No enrichments met the FDR cut-off to produce mountain plots")}
  
  ## produce volcano plot
  plot.data <- EA.Results
  if(nrow(plot.data[plot.data$p_value==0,])>0){plot.data[plot.data$p_value==0,]$p_value <- 0.00099}
  # define x, y limits
  limit.x <- ceiling(max(abs(as.numeric(plot.data$NES)), na.rm=TRUE))
  limit.y <- ceiling(max(-as.numeric(log(plot.data$p_value, 10)), na.rm=TRUE))
  # categorize data by significance level if there are significant hits
  if(nrow(significant.hits) > 0){
    plot.data$Significance <- paste0("FDR > ", FDR)
    plot.data[plot.data$FDR_q_value < FDR,]$Significance <- paste0("FDR < ", FDR)
    plot.data$Significance <- factor(plot.data$Significance,levels=c(paste0("FDR < ", FDR), paste0("FDR > ", FDR)))
    volc <- ggplot2::ggplot(data = plot.data, aes(x = NES, y = -log(p_value,10), color=Significance)) + ggplot2::geom_point(size = 4) + 
      ggrepel::geom_text_repel(data=subset(plot.data,Significance==paste0("FDR < ", FDR)),mapping=aes(label=Drug_set,size=I(4)),nudge_y = 0.25) +
      ggplot2::scale_color_manual(values=c("red","azure4"),name="Significance",breaks=c(paste0("FDR < ", FDR), paste0("FDR > ", FDR))) +
      ggplot2::xlim(-limit.x,limit.x) + ylim(0,limit.y) + ggplot2::xlab("Normalized Enrichment Score") + ggplot2::ylab("-Log(p-value)") +
      ggplot2::geom_vline(xintercept=0,linetype="solid",color="grey",size=0.5) +
      ggplot2::theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = 'black', size = 0.65),
            legend.text=element_text(size=10),axis.text=element_text(size=10),axis.title=element_text(size=20,face="bold"),
            panel.background = element_rect(fill="white", colour="white", size=0.5,linetype="solid", color="black"), text = element_text(size = 10),
            legend.position = "bottom", legend.key = element_blank())
  }else{
    volc <- ggplot2::ggplot(data = plot.data, aes(x = NES, y = -log(p_value,10))) + ggplot2::geom_point(size = 4, color="azure4") + 
      ggplot2::xlim(-limit.x,limit.x) + ggplot2::ylim(0,limit.y) + ggplot2::xlab("Normalized Enrichment Score") + ggplot2::ylab("-Log(p-value)") +
      ggplot2::geom_vline(xintercept=0,linetype="solid",color="grey",size=0.5) +
      ggplot2::theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = 'black', size = 0.65),
            axis.text=element_text(size=10),axis.title=element_text(size=20,face="bold"),
            panel.background = element_rect(fill="white", colour="white", size=0.5,linetype="solid", color="black"), text = element_text(size = 10))
  }
  
  outputs <- list(gmt = gmt, result = EA.Results, mtn.plots = temp.plot, volcano.plot = volc)
  return(outputs)
}

DMEA <- function(drug.sensitivity, gmt=NULL, expression, weights, value="AUC", sample.names=colnames(expression)[1],
                 gene.names=colnames(weights)[1], weight.values=colnames(weights)[2],
                 rank.metric="Pearson.est", FDR=0.25, num.permutations=1000, stat.type="Weighted", 
                 drug.info=NULL, drug="Drug", set.type="moa", min.per.set=6, sep="[|]", exclusions=c("-666", "NA", "na","NaN", "NULL"), descriptions=NULL,
                 min.per.corr=3, scatter.plots=TRUE, scatter.plot.type="pearson",FDR.scatter.plots=0.05, 
                 xlab="Weighted Voting Score", ylab=value, position.x="min", position.y="min", se=TRUE){
  # Weighted Voting (WV)
  WV.result <- WV(expression=expression,weights=weights,sample.names=sample.names,gene.names=gene.names,weight.values=weight.values)

  # merge drug sensitivity dataframe with WV
  WV.result.drug.sensitivity <- merge(WV.result,drug.sensitivity,by=sample.names)

  # rank.corr 
  corr.results <- rank.corr(data=WV.result.drug.sensitivity,variable="Drug",value=value,type=scatter.plot.type,min.per.corr=min.per.corr,
                            plots=scatter.plots,FDR=FDR.scatter.plots,xlab=xlab,ylab=ylab,position.x=position.x,position.y=position.y,se=se)

  # if necessary, generate gmt object
  if(is.null(gmt)){
    if(is.null(drug.info)){
      stop("drug.info dataframe containing set membership information must be provided as input if no gmt object is provided")
      }else{
      corr.output <- merge(corr.results$result, drug.info, by=drug)
      gmt <- as.gmt(data=corr.output, element.names=drug, set.names=set.type, min.per.set, sep, exclusions, descriptions)
    }
  }else{
    corr.output <- corr.results$result
  }
  
  # Drug Mechanism Enrichment Analysis (DMEA)
  DMEA.results <- drugSEA(data=corr.output, gmt, rank.metric=rank.metric, stat.type=stat.type, num.permutations=num.permutations, FDR=FDR)

  return(list(WV.scores = WV.result,
              corr.result = corr.results$result,
              corr.scatter.plots = corr.results$scatter.plots,
              gmt = gmt,
              DMEA.result = DMEA.results$result,
              DMEA.mtn.plots = DMEA.results$mtn.plots,
              DMEA.volcano.plot = DMEA.results$volcano.plot))
}
