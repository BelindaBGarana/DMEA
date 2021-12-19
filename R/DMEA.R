#DMEA
#BG 20201203; last edit: BG 20211218
#Note: drugSEA co-authored with JJ (GSEA_custom by JJ & revised by BG for drugSEA; gsea_mountain_plot by JJ & revised by BG)
#Note: thanks to NG for ng.theme (used in rank.corr)

rm(list=ls(all=TRUE))

WV <- function(expression, weights, sample.names=colnames(expression)[1],
                gene.names=colnames(weights)[1], weight.values=colnames(weights)[2]){#data: cells should be rows, genes should be columns; weights: one column should be gene names, another column should be gene weights
  library(dplyr);

  cores <- parallel::detectCores() #number of cores available
  cl <- snow::makeCluster(cores[1]-1) #cluster using all but 1 core
  doSNOW::registerDoSNOW(cl) #register cluster

  rownames(expression) <- expression[,c(sample.names)]
  input.df <- expression %>% select(weights[,c(gene.names)])
  WV.matrix <- input.df #make a dataframe to store resulting WV score matrix
  cells <- unique(rownames(WV.matrix))
  genes <- unique(colnames(WV.matrix))
  WV.scores <- as.data.frame(cells) #make a dataframe to store resulting WV scores
  WV.scores$WV <- NA
  for(m in 1:length(cells)) {
    for(j in 1:length(genes)) {
      temp.gene <- input.df %>% select(genes[j]) #select expression for particular gene j
      temp.weight <- weights[weights[,c(gene.names)]==genes[j],c(weight.values)]
      temp.gene$cells <- rownames(temp.gene) #add column with cell
      temp.gene.cell <- temp.gene[which(temp.gene$cells==cells[m]),] #get expression for cell line m for gene j
      temp.gene.cell$cells <- NULL #remove cell column
      temp.WV <- as.numeric(temp.gene.cell[1,1])*as.numeric(temp.weight) #multiply expression value by log2FC
      WV.matrix[m,j] <- temp.WV #store WV value for each gene, cell combo
    }
    WV.scores[WV.scores$cells==cells[m],]$WV <- sum(as.numeric(WV.matrix[m,])) #sum WV values across genes for each cell line
  }
  colnames(WV.scores) [1] <- sample.names

  snow::stopCluster(cl) #stop cluster
  rm(cl)

  outputs <- list(scores = WV.scores, matrix = WV.matrix)
  return(outputs)
}

rank.corr <- function(data, variable="Gene", value="Intensity",type="pearson", N.min=3, plots=TRUE, FDR=0.05, xlab=rank.var, ylab=value, position.x="mid", position.y="max", se=TRUE){
  library(dplyr);library(qvalue);library(ggplot2);
  #Correlations
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
    a <- nrow(data.corr[[j-2]]) #number of samples in the correlation (N)
    if (a >= N.min) {corr$N[j-2] <- a
    x <- as.numeric(data.corr[[j-2]][,1]) #list of rank metric
    y <- as.numeric(data.corr[[j-2]][,2]) #list of gene expression (for each gene j)

    Regression <- lm(y ~ x)
    corr$Slope[j-2] <- Regression$coeff[[2]]
    corr$Intercept[j-2] <- Regression$coeff[[1]]
    corr$R.squared[j-2] <- summary(Regression)$r.squared

    Rank.regression <- lm(rank(y) ~ rank(x))
    corr$Rank.slope[j-2] <- Rank.regression$coeff[[2]]
    corr$Rank.intercept[j-2] <- Rank.regression$coeff[[1]]
    corr$Rank.R.squared[j-2] <- summary(Rank.regression)$r.squared

    Pearson <- cor.test(x, y, method="pearson")
    corr$Pearson.est[j-2] <- Pearson$estimate
    corr$Pearson.p[j-2] <- Pearson$p.value

    Spearman <- cor.test(x, y, method="spearman", exact=FALSE)
    corr$Spearman.est[j-2] <- Spearman$estimate
    corr$Spearman.p[j-2] <- Spearman$p.value}
  }
  corr.no.na <- corr[!is.na(corr$N), ]
  qobj.pearson <- qvalue(p=corr.no.na$Pearson.p,pi0=1)
  qobj.spearman <- qvalue(p=corr.no.na$Spearman.p,pi0=1)
  corr.no.na$Pearson.q <- qobj.pearson$qvalues
  corr.no.na$Spearman.q <- qobj.spearman$qvalues

  #Scatter plots for significant results
  if (plots==TRUE){
    library(reshape2); library(ggplot2); library(gridExtra);
    #Load themes for plots
    ng.theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.border = element_rect(fill=NA), panel.background = element_blank(),
                      axis.line = element_line(colour = "black"), axis.text.x = element_text(colour = "black"),
                      axis.text.y = element_text(colour = "black"), axis.ticks.x = element_line(colour="black"),
                      axis.ticks.y = element_line(colour="black"), legend.title = element_blank(),
                      axis.title.y = element_text(size=8, colour="black"))

    bg.theme <- theme(legend.background = element_rect(), legend.position="top", legend.text = element_text(size=14), legend.key = element_blank(),
                      axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                      axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                      plot.title = element_text(lineheight=.8, face="bold", size=36))
    #Prepare dataframe for plots
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

          a[[i]] <- ggplot(data=sig.data, aes_string(x=rank.var, y=value)) +
            geom_point() + labs(x = xlab, y = ylab) + ggtitle(results[i]) +
            geom_smooth(method="lm",size = 1.5,linetype = 'solid',color="blue",se = se,na.rm = TRUE) +
            geom_text(x=pos.x,y=pos.y,vjust="inward",hjust="inward", colour="blue", parse=TRUE,
                       label = as.character(as.expression(stats_pearson)), size = 8) + ng.theme + bg.theme
        }
        scatter.plots <- marrangeGrob(a, nrow=1, ncol=1)
      }else{
        scatter.plots <- NA
        print("No correlations met the FDR cut-off to produce scatter plots")
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
        scatter.plots <- marrangeGrob(a, nrow=1, ncol=1)
      }else{
        scatter.plots <- NA
        print("No correlations met the FDR cut-off to produce scatter plots")
      }
    }else{
      print("type must be specified as either spearman or pearson to produce scatter plots")
    }
  }else{scatter.plots <- NA}
  outputs <- list(result=corr.no.na, scatter.plots = scatter.plots)
  return(outputs)
}

drugSEA <- function(data, gmt, drug="Drug", estimate="Pearson.est", set.type="moa", direction.adjust=NULL,
                    FDR=0.25, num.permutations=1000, stat.type="Weighted"){
  library(sjmisc);library(dplyr);
  #adjust estimates
  if(length(direction.adjust)>0){
    info <- distinct(data[,c(drug, estimate, set.type)])
    info$adjust <- 1
    info[grep(direction.adjust, info[,c(set.type)], value = T), "adjust"] <- -1
    info$adjusted.est <- as.numeric(info[,c(estimate)])*as.numeric(info$adjust)
    corr.adj <- info[,c(drug,"adjusted.est")]
    GSEA.sample.name <- "Adjusted Estimate"
    colnames(corr.adj) <- c("Gene",GSEA.sample.name)
    input <- corr.adj
  }else{
    GSEA.sample.name <- "Estimate"
    input <- data[,c(drug,estimate)]
    colnames(input) <- c("Gene",GSEA.sample.name)
  }

  #load GSEA_custom function
  GSEA_custom <- function(input.df, gmt.list,
                          num.permutations = 1000,
                           stat.type = "Weighted"){
    loop.time <- Sys.time()

    nperm = num.permutations #number of permutations
    if (stat.type == "Classic"){
      score.weight = 0
    }
    if (stat.type == "Weighted"){
      score.weight = 1
    }


    #Read in gene expression data
    #Genes should be first column, named "Gene"
    #Samples should be columns 2:N
    data_in <- input.df

    gmt.for.reformat <- gmt.list
    Gene.Sets <- t(plyr::ldply(gmt.for.reformat$genesets, rbind)) #reformat gmt list to desired format
    colnames(Gene.Sets) <- gmt.for.reformat$geneset.names

    Gene.Sets <- as.data.frame(Gene.Sets)

    testthat::expect_is(data_in, "data.frame")
    testthat::expect_is(Gene.Sets, "data.frame")

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
        #      ES <- max.ES
        ES <- signif(max.ES, digits = 5)
        arg.ES <- which.max(RES)
      } else {
        #      ES <- min.ES
        ES <- signif(min.ES, digits=5)
        arg.ES <- which.min(RES)
      }
      return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicators))
    } #for real ES

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

      if (weighted.score.type == 0) {
        tag.correl.vector <- rep(1, Nh)
      } else if (weighted.score.type == 1) {
         tag.correl.vector <- correl.vector[tag.loc.vector]
         tag.correl.vector <- abs(tag.correl.vector)
      } else if (weighted.score.type == 2) {
         tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
         tag.correl.vector <- abs(tag.correl.vector)
       } else {
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

    Samples <- colnames(data_in)
    if (Samples[1] != "Gene"){
      stop("Please ensure that your data frame is organized with the first column to be named 'Gene'")
    }
    Samples <- Samples[-1]

    Gene.Sets.All <- colnames(Gene.Sets)

    annotations <- matrix(data = 0, nrow = nrow(data_in), ncol = length(Gene.Sets.All))
    colnames(annotations) <- Gene.Sets.All

    annotations <- as.data.frame(annotations)

    annotations <- cbind(data_in$Gene,annotations)
    colnames(annotations) <- c("Gene", Gene.Sets.All)
    annotations <- as.matrix(annotations)

    num.hits.pathways <- list()

    ### Annotate gene sets
    for (j in 1:length(Gene.Sets.All)){
      temp.pathway <- Gene.Sets[,Gene.Sets.All[j]]
      for (i in 1:nrow(annotations)){
        if (annotations[i,"Gene"] %in% temp.pathway){
          annotations[i,j+1] = "X";
        }
      }
      num.hits.pathways[[Gene.Sets.All[j]]] <- sum(annotations[,Gene.Sets.All[j]] == "X")
    }

    num.hits.pathways.df <- matrix(unlist(num.hits.pathways))
    row.names(num.hits.pathways.df) = Gene.Sets.All
    num.gene.sets.under.5 <- which(num.hits.pathways.df < 5)
    num.gene.sets.over.4 <- which(num.hits.pathways.df > 4)
    Gene.Sets.All <- Gene.Sets.All[num.gene.sets.over.4] # write over gene sets all to keep > 4
    if (length(num.gene.sets.under.5) > 1){
      print("Warning: Removing drug sets with less than 5 drugs observed in data set.")
      #gene.sets.to.remove <- Gene.Sets.All[num.gene.sets.under.5]
      annotations <- annotations[,c("Gene",Gene.Sets.All)]
      #annotations[,which(colnames(annotations) %in% gene.sets.to.remove)] <- NULL
    }
    annotations <- as.data.frame(annotations)
    data_in <- merge(data_in, annotations, by = "Gene")

    data_in <- stats::na.omit(data_in)

    GSEA.Results.All.Samples <- matrix(data = NA, nrow = 0, ncol = 7)
    colnames(GSEA.Results.All.Samples) <- c("Sample","Gene.Set","KS","KS_Normalized",
                                            "p-value","Position at Max",
                                            "FDR q-value")
    Mountain.Plot.Info.All.Samples <- list()
    rank_metric.All.Samples <- list()

    #Find out how many cores are available (if you don't already know)
    cores<-parallel::detectCores()
    #Create cluster with desired number of cores, leave one open for the machine
    #core processes
    cl <- snow::makeCluster(cores[1]-1)
    #Register cluster
    doSNOW::registerDoSNOW(cl)

    rm(annotations)
    data_in2 <- array(data = NA)
    for (u in 1:length(Samples)){

      data_in2 <- cbind(subset(data_in, select = Gene.Sets.All),
                        dplyr::select(data_in, Samples[u]))  #select one Sample type and the genes and Gene.Sets.A.and.B
      data_in2[,Samples[u]] <- as.numeric(as.character(data_in2[,Samples[u]]))
      data_in2 <- data_in2[order(-data_in2[,Samples[u]]),] #sort by descending order for the rank metric
      rownames(data_in2) <- 1:nrow(data_in2) #reorder row indices for counting in for loop below

      ## Assuming first two columns in data table are Genes and Rank Metric (e.g. Foldchange, SNR)

      GSEA.Results <- matrix(data = NA, nrow = length(Gene.Sets.All), ncol = 7)
      colnames(GSEA.Results) <- c("Sample","Gene.Set","KS","KS_Normalized",
                                  "p_value","Position_at_max",
                                  "FDR_q_value")
      GSEA.Results <- as.data.frame(GSEA.Results)
      GSEA.Results$Gene.Set <- Gene.Sets.All
      GSEA.Results$Sample <- Samples[u]

      ions <- nrow(data_in2)

      #for plotting
      ks_results_plot <- list()
      positions.of.hits <- list()

      #ks_results_plot <- as.data.frame(ks_results_plot)
      gene.list <- 1:ions
      rank_metric <- data_in2[,Samples[u]] #Save the rank metric

      pos_gene_set <- array(data = 0, dim = nrow(data_in2), dimnames = NULL);

      ## Calculate Real KS Statistic
      for (i in 1:length(Gene.Sets.All)){
        data_in3 <- data_in2[,Gene.Sets.All[i]]
        numhits_pathway <- sum(data_in3 == "X"); #check to see if there is anything in the column (e.g. X)
        if (numhits_pathway > 1){
          pos_gene_set <- which(data_in2[,Gene.Sets.All[i]] %in% c("X"))
          KS_real <- GSEA.EnrichmentScore(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
          GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$KS <- KS_real$ES;
          GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$Position_at_max <- KS_real$arg.ES;
          ks_results_plot[[Gene.Sets.All[i]]] = KS_real$RES
          positions.of.hits[[Gene.Sets.All[i]]] = pos_gene_set
        }
      }

      Mountain.Plot.Info <- list(MountainPlot = ks_results_plot, Position.of.hits = positions.of.hits)
      rm(pos_gene_set)
      rm(numhits_pathway)
      rm(data_in3)
      rm(KS_real)

      #print("Calculating permutations...")

      #pb <- utils::txtProgressBar(max = num.permutations, style = 3)
      #progress <- function(n) utils::setTxtProgressBar(pb, n)
      #opts <- list(progress = progress)

      KSRandomArray <- matrix(data = NA, nrow = nperm, ncol = length(Gene.Sets.All))
      num.gene.sets.all <- length(Gene.Sets.All)
      `%dopar%` <- foreach::`%dopar%`

      #KSRandomArray <- foreach::foreach(L = 1:nperm, .combine = "rbind",.options.snow = opts) %dopar% {
      KSRandomArray <- foreach::foreach(L = 1:nperm, .combine = "rbind") %dopar% {
        temp.KSRandomArray <- matrix(data = NA, nrow = 1, ncol = num.gene.sets.all)
        for(i in 1:length(Gene.Sets.All)){
          numhits_pathway <- length(positions.of.hits[[Gene.Sets.All[i]]])
          pos_gene_set <- sample(1:ions,numhits_pathway)
          temp.KSRandomArray[,i] <- GSEA.EnrichmentScore2(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
        }
        temp.KSRandomArray
      }
      colnames(KSRandomArray) <- Gene.Sets.All

      #rm(opts)
      #rm(pb)
      KSRandomArray <- data.frame(matrix(unlist(KSRandomArray), nrow = nperm, byrow = T))
      colnames(KSRandomArray) <- Gene.Sets.All
      KSRandomArray <- stats::na.omit(KSRandomArray)

      #print("Normalizing enrichment scores...")
      KSRandomArray <- as.data.frame(KSRandomArray)
      ###normalize the GSEA distribution
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
      percent.pos.GSEA <- sum(GSEA.Results$KS >= 0) / length(GSEA.Results$KS) #BG OR EQUAL TO
      percent.neg.GSEA <- sum(GSEA.Results$KS <= 0) / length(GSEA.Results$KS) #BG OR EQUAL TO

      # Calculate GSEA NES and p-value
      #print("Calculating GSEA NES and p-value...")
      for (i in 1:length(Gene.Sets.All)){
        temp.gene.set <- Gene.Sets.All[i]
        temp.KS <- GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS
        if (temp.KS >= 0){ #BG OR EQUAL TO
          pos.perms <- KSRandomArray[,temp.gene.set]
          pos.perms <- pos.perms[which(pos.perms >= 0)] #BG OR EQUAL TO
          #p-val
          GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$p_value = signif(sum(pos.perms >= temp.KS) / length(pos.perms),digits = 3) #BG OR EQUAL TO
          #NES
          GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized = signif(temp.KS / mean(pos.perms), digits = 3)
        } else if (temp.KS <= 0){ #BG OR EQUAL TO
          neg.perms <- KSRandomArray[,temp.gene.set]
          neg.perms <- neg.perms[which(neg.perms <= 0)] #BG OR EQUAL TO
          #p-val
          GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$p_value = signif(sum(neg.perms <= temp.KS) / length(neg.perms),digits = 3) #BG OR EQUAL TO
          #NES
          GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized = signif(temp.KS / mean(neg.perms) * -1, digits = 3)
        }
      }

      # Calculate GSEA FDR
      #print("Calculating GSEA FDR...")
      for (i in 1:length(Gene.Sets.All)){
        temp.gene.set <- Gene.Sets.All[i]
        temp.NES <- GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized
        if (temp.NES >= 0){ #BG OR EQUAL TO
          #FDR
          percent.temp <- sum(GSEA.NES.perms.pos >= GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized) / length(GSEA.NES.perms.pos) #BG OR EQUAL TO
          percent.pos.stronger <- sum(GSEA.Results$KS_Normalized >= GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized) / sum(GSEA.Results$KS_Normalized >= 0) #BG
          GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$FDR_q_value = ifelse(signif(percent.temp / percent.pos.stronger, digits = 3) < 1, signif(percent.temp / percent.pos.stronger, digits = 3), 1) #BG
        } else if (temp.NES <= 0){ #BG OR EQUAL TO
          #FDR
          percent.temp <- sum(GSEA.NES.perms.neg <= temp.NES) / length(GSEA.NES.perms.neg) #BG OR EQUAL TO
          percent.neg.stronger <- sum(GSEA.Results$KS_Normalized <= temp.NES) / sum(GSEA.Results$KS_Normalized <= 0) #BG
          GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$FDR_q_value = ifelse(signif(percent.temp / percent.neg.stronger, digits = 3) < 1, signif(percent.temp / percent.neg.stronger, digits = 3), 1) #BG
        }
      }

      GSEA.Results.All.Samples <- rbind(GSEA.Results.All.Samples,GSEA.Results)
      Mountain.Plot.Info.All.Samples <- c(Mountain.Plot.Info.All.Samples,Mountain.Plot.Info)
      rank_metric.All.Samples <- c(rank_metric.All.Samples,rank_metric)

      #print(paste("Sample #: ", u))

      #end.loop.time <- Sys.time()
      #total.loop.time <- signif(end.loop.time - loop.time, digits = 3)
      #print(paste("Time per Sample:" , total.loop.time))
    }

    snow::stopCluster(cl)
    rm(cl)

    return(list(GSEA.Results = GSEA.Results.All.Samples,
                Mountain.Plot.Info = Mountain.Plot.Info.All.Samples,
                ranking.metric = rank_metric.All.Samples))
  }
  #load mountain plot function
  gsea_mountain_plot <- function(GSEA.list, Sample.Name, Gene.Set.A, color = TRUE){
    library(ggplot2);
    if (color == TRUE){
      color.palette <- c("red")
    } else if (color == FALSE){
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
    GSEA.results <- GSEA.results[GSEA.results$Sample == Sample.Name,]
    GSEA.results$p_value <- as.numeric(as.character(GSEA.results$p_value))
    GSEA.results$FDR_q_value <- as.numeric(as.character(GSEA.results$FDR_q_value))

    gsea.normalized.enrichment.score <- signif(GSEA.results[GSEA.results$Gene.Set == myName,]$KS_Normalized, digits =3)
    gsea.p.value <- signif(GSEA.results[GSEA.results$Gene.Set == myName,]$p_value, digits = 3)
    gsea.q.val <- signif(GSEA.results[GSEA.results$Gene.Set == myName,]$FDR_q_value, digits = 3)

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
      #cowplot::theme_nothing() +
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

    if (!is.null(gsea.normalized.enrichment.score) & length(gsea.normalized.enrichment.score) > 0){
      if (gsea.normalized.enrichment.score > 0){
        mtn.plot <- mtn.plot + ggplot2::annotate("text",label = results.as.text, x = ngenes * 0.90, y= max.ES * 0.70, color = "black")
      } else if (gsea.normalized.enrichment.score < 0){
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
  #perform GSEA using GSEA_custom
  GSEA <- GSEA_custom(input,gmt,num.permutations=num.permutations,stat.type=stat.type)
  GSEA.Results <- GSEA$GSEA.Results
  #select results which meet FDR threshold and produce mountain plots
  significant.hits <- GSEA.Results[which(GSEA.Results$FDR_q_value<=FDR),]
  temp.plot <- list()
  if (nrow(significant.hits) > 0){
    for (i in 1:nrow(significant.hits)){
      temp <- gsea_mountain_plot(GSEA.list = GSEA, Sample.Name = GSEA.sample.name ,Gene.Set.A = significant.hits$Gene.Set[i])
      temp.plot[i] <- list(temp)
    }
    colnames(GSEA.Results)[colnames(GSEA.Results) == "Gene.Set"] <- "Drug.Set"
  } else {print("No enrichments met the FDR cut-off to produce mountain plots")
    colnames(GSEA.Results)[colnames(GSEA.Results) == "Gene.Set"] <- "Drug.Set"}
  outputs <- list(result=GSEA.Results, mtn.plots=temp.plot)
  return(outputs)
}

DMEA <- function(drug.sensitivity, gmt, expression, weights, value="AUC", expr.sample.names=colnames(expression)[1],
                 gene.names=colnames(weights)[1], weight.values=colnames(weights)[2],
                 estimate="Pearson.est", FDR=0.25, num.permutations=1000, stat.type="Weighted", N.min=3,
                 scatter.plots=TRUE, FDR.scatter.plots=0.05, xlab="Weighted Voting Score", ylab=value, position.x="min", position.y="min", se=TRUE){
  print("Calculating Weighted Voting scores...")
  #WV
  WV.result <- WV(expression=expression,weights=weights,sample.names=expr.sample.names,gene.names=gene.names,weight.values=weight.values)$scores

  #merge PRISM with WV
  WV.result.drug.sensitivity <- merge(WV.result,drug.sensitivity,by=expr.sample.names)

  print("Running correlations and regressions...")
  #rank.corr
  corr.results <- rank.corr(data=WV.result.drug.sensitivity,variable="Drug",value=value,N.min=N.min,plots=scatter.plots,FDR=FDR.scatter.plots,
                            xlab=xlab,ylab=ylab,position.x=position.x,position.y=position.y,se=se)

  print("Running enrichment analysis...")
  #drugSEA
  DMEA.results <- drugSEA(data=corr.results$result,gmt=gmt,estimate=estimate,stat.type=stat.type,num.permutations = num.permutations,FDR=FDR)

  return(list(WV.scores = WV.result,
              corr.result = corr.results$result,
              corr.scatter.plots = corr.results$scatter.plots,
              DMEA.result = DMEA.results$result,
              DMEA.mtn.plots = DMEA.results$mtn.plots))
}
