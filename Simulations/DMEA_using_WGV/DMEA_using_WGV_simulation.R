# Authors: Belinda B. Garana (BG), James H. Joly (JJ)
# JJ sent 12-17-20; BG adapted for DMEA; last edit: BG 2021-10-02
### We want to know the true positive / false negative rates of enrichment
### So we will create synthetic data with enrichments for both gene WGV and drug AUC.
### How?

### We will perturb gene expression and drug AUC in the same directions, then calculate:
### 1) WGV
### 2) Correlation coefficients for each drug
### 3) DMEA

### How to build this?

### Step 1: Simulate gene expression matrix for 200 cell lines representative of CCLE normalized RNAseq version 19Q4, then perturb 1 set of 25 genes +/- X
###         Note: need to name the cell lines such that it's easy to track, we will want to perturb drug AUC in the same direction
### Step 2: Calculate WGV scores for synthetic gene set
### Step 3: Simulate drug AUC for 200 cell lines representative of PRISM, then perturb 1 drug set of 10 drugs +/- Y
### Step 4: Run correlations (drug AUC vs. WGV)
### Step 5: Run DMEA for each different variant (+/- X gene expression, +/- Y drug AUC)
### Step 6: Re-run for 50 replicates

rm(list=ls(all=TRUE))
library(DMEA);library(dplyr);library(GSA);library(reshape2);library(data.table);
library(ggplot2);library(stringr);library(dostats);

### Step 0: Prep rows, columns, values
# Import pathways
drug.sets <- GSA.read.gmt(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Simulations/DMEA_using_WGV/MOA_gmt_file_n6_w_synthetic_set.gmt")
pathway.names <- drug.sets$geneset.names
num.pathways <- length(pathway.names)

# Import CCLE RNAseq data for gene list
download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin", 
              destfile = "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin", 
              destfile = "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
CCLE.RNAseq <- rbind(RNA.first200, RNA.rest)
gene.names <- unique(colnames(CCLE.RNAseq[,2:ncol(CCLE.RNAseq)]))
num.genes <- length(gene.names) #19144

# Make a synthetic gene set of size 25 genes
gene.set.size <- 25
#synthetic.gene.set <- sample(gene.names,gene.set.size) #we did this once and now we want to use the same gene set for each time, so lets save it separately
synthetic.gene.set <-  c("HDAC11","FSBP","PGF","NBEAL2","SMIM15","PRSS33","USP29","ZCCHC7","PIK3C2B","ZFAND2B","CLMP","APOL4","TXNIP",
                         "LRRC41","ARHGEF38","KRT85","TMEM121","DPYSL2","DAZ4","C1QA","CILP","AFF2","CHAC2","BLVRA","RNF114")

#Import PRISM data for drug list
prism <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv",header=T)
prism$X <- NULL

#Make a synthetic drug set of size 10 drugs - already done in "Metabolic Pathways and the drugs that target them - filtered - gmt.gmt"
drug.set.size <- 10
#synthetic.drug.set <- sample(drug.names,drug.set.size) # we did this once and now we want to use the same gene set for each time, so lets save it separately

synthetic.drug.set <-  c("blebbistatin....","BI.D1870","PF.05212384","oligomycin.a","levonorgestrel","norfloxacin","poziotinib",
                         "estradiol","CGP.54626","epirubicin")

# Make synthetic cell line names
synthetic.cell.names <- seq(from = 1, to = 200, by = 1)
synthetic.cell.names <- paste("Cell_No_", synthetic.cell.names, sep = "")

# How far do we want to vary values here?
values.to.vary <- seq(from = 0, to = 0.1, by = 0.01)

# Import drug info
drug.info <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv",header=T) #481 cell lines
colnames(drug.info)[colnames(drug.info)=="name"] <- "Drug"
drug.info$Drug <- gsub("[[:punct:][:blank:]]",".",drug.info$Drug)
drug.moa <- na.omit(distinct(drug.info[,c("Drug","moa")]))
Drug <- unique(drug.moa$Drug)
moa.drugs <- data.frame(Drug)
Drug <- unique(colnames(prism[,2:ncol(prism)]))
prism.drugs <- data.frame(Drug)
overlap.drugs <- merge(moa.drugs,prism.drugs,by="Drug")

# reduce PRISM data to drugs with moa
prism <- prism %>% select("CCLE_ID",overlap.drugs$Drug)
drug.names <- unique(colnames(prism[,2:ncol(prism)]))
num.drugs <- length(drug.names) #1351

# set drug AUC parameters
frac.low.AUC <- 0.72
num.drugs.low <- round(frac.low.AUC*num.drugs,digits=0) #973
num.drugs.high <- round((1-frac.low.AUC)*num.drugs,digits=0) #378
num.synth.drugs.low <- round(frac.low.AUC*drug.set.size,digits=0) #7
num.synth.drugs.high <- round((1-frac.low.AUC)*drug.set.size,digits=0) #3
low.mean.AUC <- 0.83
low.sd.AUC <- 0.1657471
high.mean.AUC <- 1.31
high.sd.AUC <- 0.1580992
half.low.sd.AUC <- 0.08
half.high.sd.AUC <- 0.08

# Run simulations (N replicates)
N.max <- 50
for (N in 1:N.max){
  loop.start.time <- Sys.time()

  ### Step 1: Generate gene expression matrix for 200 cell lines and 25 genes in synthetic gene set
  # Step 1a: make gene expression matrix
  expr.matrix <- matrix(data = NA, nrow = gene.set.size, ncol = length(synthetic.cell.names)+1)
  row.names(expr.matrix) <- synthetic.gene.set
  colnames(expr.matrix) <- c("Gene",synthetic.cell.names)
  expr.matrix[,"Gene"] <- synthetic.gene.set
  expr.matrix <- as.data.frame(expr.matrix)

  # Step 1b: make list of varied value expression matrices
  expr.matrices <- list()
  varied.value.names.all <- paste("value_added_",values.to.vary, sep = "")

  # Step 1c: replace values of synthetic gene set with varied values (perturbation); Need to do 0 separately
  varied.value <- values.to.vary[1]
  varied.value.name <- paste("value_added_",varied.value, sep = "")
  for (k in 1:length(synthetic.cell.names)){
    expr.matrix[expr.matrix$Gene %in% synthetic.gene.set, synthetic.cell.names[k]] <- rnorm(length(synthetic.gene.set), mean = 0, sd = 0.5)
  }
  expr.matrices[[varied.value.names.all[1]]] <- expr.matrix
  
  for(j in 2:length(values.to.vary)){
    varied.value <- values.to.vary[j]
    varied.value.name <- paste("value_added_",varied.value, sep = "")

    range.to.vary <- seq(from = -varied.value, to = varied.value, by = 2*varied.value / (length(synthetic.cell.names)-1))
    for (k in 1:length(synthetic.cell.names)){
        expr.matrix[expr.matrix$Gene %in% synthetic.gene.set, synthetic.cell.names[k]] <- rnorm(length(synthetic.gene.set), mean = range.to.vary[k], sd = 0.5)
    }
    expr.matrices[[varied.value.names.all[j]]] <- expr.matrix
  }


  ### Step 2: Calculate WGV scores for synthetic gene set
  # generate weight values by dividing +X/-X expression values for each gene in the synthetic gene set
  print("Calculating weights for synthetic gene set")
  diff.all.var <- list()
  for(i in 1:length(expr.matrices)){
    var.name <- names(expr.matrices[i])
    temp.expr <- expr.matrices[[var.name]] # select expression variation
    temp.expr.hits <- temp.expr[temp.expr$Gene %in% synthetic.gene.set, ] # select expression for synthetic gene set
    temp.expr.hits$Gene <- NULL #remove "Gene" column
    input.df <- as.data.frame(t(temp.expr.hits)) # Cells should be rows, genes should be columns
    genes <- unique(colnames(input.df))
    diff.values <- as.data.frame(genes) # make a dataframe to store resulting WGV scores
    diff.values$diff <- NA
    for(j in 1:length(genes)) {
      temp.gene <- input.df %>% select(genes[j]) #select expression for particular gene j
      temp.gene[,1] <- as.numeric(temp.gene[,1])
      # calculate diff based on top vs. bottom 20 cell lines' expression
      top20 <- temp.gene[181:200,]
      mean.top20 <- mean(as.numeric(top20))
      bottom20 <- temp.gene[1:20,]
      mean.bottom20 <- mean(as.numeric(bottom20))
      diff.values$diff[j] <- mean.top20-mean.bottom20
      }
    diff.all.var[[var.name]] <- diff.values
  }

  print("Calculating WGV")
  WGV.scores.all.var <- list()
  WGV.start.time <- Sys.time()
  for(i in 1:length(expr.matrices)){
    var.name <- names(expr.matrices[i])
    temp.diff.values <- na.omit(diff.all.var[[var.name]])
    temp.expr <- expr.matrices[[var.name]] # select expression variation
    temp.expr.hits <- temp.expr[temp.expr$Gene %in% temp.diff.values$genes, ] # select expression for synthetic gene set
    temp.expr.hits$Gene <- NULL # remove "Gene" column
    input.df <- as.data.frame(t(temp.expr.hits)) # Cells should be rows, genes should be columns
    WGV.matrix <- input.df # make a dataframe to store resulting WGV score matrix
    cells <- unique(rownames(WGV.matrix))
    genes <- unique(colnames(WGV.matrix))
    WGV.scores <- as.data.frame(cells) # make a dataframe to store resulting WGV scores
    WGV.scores$WGV <- NA
    for(m in 1:length(cells)) {
      for(j in 1:length(genes)) {
        temp.gene <- input.df %>% select(genes[j]) # select expression for particular gene j
        temp.diff <- temp.diff.values[temp.diff.values$genes==genes[j],]$diff
        temp.gene$cells <- rownames(temp.gene) # add column with cell
        temp.gene.cell <- temp.gene[which(temp.gene$cells==cells[m]),] # get expression for cell line m for gene j
        temp.gene.cell$cells <- NULL # remove cell column
        temp.WGV <- as.numeric(temp.gene.cell[1,1])*as.numeric(temp.diff) # multiply expression value by diff
        WGV.matrix[m,j] <- temp.WGV # store WGV value for each gene, cell combo
      }
      WGV.scores[WGV.scores$cells==cells[m],]$WGV <- sum(as.numeric(WGV.matrix[m,])) #sum WGV values across genes for each cell line
    }
    WGV.scores.all.var[[var.name]] <- WGV.scores
  }

  ### Step 3: Generate drug AUC that have a normal distribution for each cell line, then perturb synthetic drug set
AUC.matrix <- matrix(data = NA, nrow = length(drug.names), ncol = length(synthetic.cell.names)+1)
row.names(AUC.matrix) <- drug.names
colnames(AUC.matrix) <- c("Drug",synthetic.cell.names)
AUC.matrix[,"Drug"] <- drug.names
AUC.matrix <- as.data.frame(AUC.matrix)
# Fill out matrix with normal distribution of rank values
# And replace synthetic drug set with a value of X

for(i in 1:length(synthetic.cell.names)){
  temp.cell <- synthetic.cell.names[i]
  normal.low <- rnorm(num.drugs.low, mean=low.mean.AUC, sd=half.low.sd.AUC)
  normal.high<- rnorm(num.drugs.high, mean=high.mean.AUC, sd=half.high.sd.AUC)
  normal <- c(normal.low,normal.high)
  AUC.matrix[,temp.cell] <- normal
}


AUC.matrices <- list()

## Need to do 0 separately
varied.value <- values.to.vary[1]
varied.value.name <- paste("value_added_",varied.value, sep = "")
for (k in 1:length(synthetic.cell.names)){
  normal.low <- rnorm(num.synth.drugs.low, mean=low.mean.AUC, sd=half.low.sd.AUC)
  normal.high<- rnorm(num.synth.drugs.high, mean=high.mean.AUC, sd=half.high.sd.AUC)
  normal <- c(normal.low,normal.high)
  AUC.matrix[synthetic.drug.set, synthetic.cell.names[k]] <- normal
}
AUC.matrices[[varied.value.name]] <- AUC.matrix

for(j in 2:length(values.to.vary)){
  varied.value <- values.to.vary[j]
  varied.value.name <- paste("value_added_",varied.value, sep = "")

  range.to.vary <- seq(from = -varied.value, to = varied.value, by = 2*varied.value / (length(synthetic.cell.names)-1)) #from0 BG 20210913
  for (k in 1:length(synthetic.cell.names)){
    normal.low <- rnorm(num.synth.drugs.low, mean=low.mean.AUC+range.to.vary[k], sd=half.low.sd.AUC)
    normal.high<- rnorm(num.synth.drugs.high, mean=high.mean.AUC+range.to.vary[k], sd=half.high.sd.AUC)
    normal <- c(normal.low,normal.high)
    AUC.matrix[synthetic.drug.set, synthetic.cell.names[k]] <- normal
  }
  AUC.matrices[[varied.value.name]] <- AUC.matrix
}

### Step 4: Generate correlation coefficients (merge WGV & AUC beforehand!); took a minute or so per N
expr.matrix.names <- names(expr.matrices)
AUC.matrix.names <- names(AUC.matrices)

correlation.matrices <- list()
print("Generating Correlation Coefficients")
for(i in 1:length(expr.matrix.names)){
  for (j in 1:length(AUC.matrix.names)){
    expr.value <- paste("expr_",expr.matrix.names[i],sep = "")
    AUC.value <- paste("AUC_",AUC.matrix.names[j], sep = "")

    corr.name <- paste(expr.value,"WGV", AUC.value,"correlation", sep = "_")

    temp.WGV.results <- WGV.scores.all.var[[expr.matrix.names[i]]]
    temp.AUC.matrix <- AUC.matrices[[AUC.matrix.names[j]]]
    temp.AUC.matrix$Drug <- NULL
    temp.AUC.matrix <- t(temp.AUC.matrix)
    input.df <- cbind(temp.WGV.results,temp.AUC.matrix)

    temp.corr <- rank.corr(data=input.df,variable="Drug",value="AUC",plots=F)$result
    correlation.matrices[[corr.name]] <- temp.corr
  }
}

##### Step 5: Run DMEA
print("Running DMEA")
result.matrices <- list()
Correlation <- names(correlation.matrices)
synthetic.results <- as.data.frame(Correlation)
synthetic.results$ES <- NA
synthetic.results$NES <- NA
synthetic.results$p <- NA
synthetic.results$q <- NA
drugSEA.start.time <- Sys.time()
for(i in 1:length(correlation.matrices)){
  corr.name <- names(correlation.matrices[i])
  temp.corr <- correlation.matrices[[corr.name]]
  temp.result <- drugSEA(data = temp.corr,gmt = drug.sets, estimate = "Pearson.est", FDR=0)$result
  result.matrices[[corr.name]] <- temp.result
  # store synth drug set results in a separate df and save that
  synthetic.results[synthetic.results$Correlation==corr.name,]$ES <- temp.result[temp.result$Drug.Set=="Synthetic_Drug_Set",]$ES
  synthetic.results[synthetic.results$Correlation==corr.name,]$NES <- temp.result[temp.result$Drug.Set=="Synthetic_Drug_Set",]$NES
  synthetic.results[synthetic.results$Correlation==corr.name,]$p <- temp.result[temp.result$Drug.Set=="Synthetic_Drug_Set",]$p_value
  synthetic.results[synthetic.results$Correlation==corr.name,]$q <- temp.result[temp.result$Drug.Set=="Synthetic_Drug_Set",]$FDR_q_value
}
all.results <- rbindlist(result.matrices,use.names=TRUE,idcol="correlation",fill=TRUE)
write.csv(all.results,file=paste0(N,"_", drugSEA.start.time,"_all_drugSEA_results.csv"))
write.csv(synthetic.results,file=paste0(N,"_", drugSEA.start.time,"_Synthetic_Drug_Set_drugSEA_results.csv"))
}

## Generate tile plots
# Load themes for plots
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

# import results
sim.results <- import_list(dir(pattern = "Synthetic_Drug_Set_drugSEA_results.csv$"))

# rbind all dataframes in sim.resuls; make sure you only use N.max replicates
plot.data <- rbindlist(sim.results,use.names=TRUE,idcol="file",fill=TRUE)

# formatting: make columns for replicate number N, expr.value, and AUC.value
expr.values <- seq(0,0.1,0.01)
AUC.values <- seq(0,0.1,0.01)
for (i in 1:nrow(plot.data)){
  plot.data$N[i] <- substring(plot.data$file[i],1,2)
  plot.data$expr.string[i] <- substring(plot.data$Correlation[i],1,23)
  for (j in 1:length(expr.values)){
    if (grepl(expr.values[j], plot.data$expr.string[i], fixed=TRUE)==TRUE) {
      plot.data$expr.value[i] <- expr.values[j]
    }
  }
  plot.data$AUC.string[i] <- substring(plot.data$Correlation[i],24,50)
  for (k in 1:length(AUC.values)){
    if (grepl(AUC.values[k], plot.data$AUC.string[i], fixed=TRUE)==TRUE){
      plot.data$AUC.value[i] <- AUC.values[k]
    }
  }
}

#For each correlation, average NES and count nrow q<0.05 and divide by N.max for percent
plot.data[,c("Average.NES","SD.NES",
             "Average.q","SD.q",
             "percent.q0.05","percent.q0.25")] <- NA
Correlations <- unique(plot.data$Correlation)
for(i in 1:length(Correlations)){
  temp.data <- plot.data[plot.data$Correlation==Correlations[i],]
  reps.q0.05 <- nrow(temp.data[temp.data$q<0.05,])
  reps.q0.25 <- nrow(temp.data[temp.data$q<0.25,])
  for(j in 1:nrow(plot.data)){
    if(plot.data$Correlation[j]==Correlations[i]){
      plot.data$Average.NES[j] <- mean(temp.data$NES)
      plot.data$SD.NES[j] <- sd(temp.data$NES)
      plot.data$Average.q[j] <- mean(temp.data$q)
      plot.data$SD.q[j] <- sd(temp.data$q)
      plot.data$percent.q0.05[j] <- as.numeric(reps.q0.05*100/N.max)
      plot.data$percent.q0.25[j] <- as.numeric(reps.q0.25*100/N.max)
    }
  }
}
write.csv(plot.data,file="WGV_drugSEA_sim_N50.csv")
plot.data <- read.csv(file = "WGV_drugSEA_sim_N50.csv")
plot.data$X <- NULL

q0.05.tile.plot <- ggplot(plot.data, aes(x = expr.value, y = AUC.value, fill = percent.q0.05)) +
  geom_tile() + ng.theme + theme_light() + scale_fill_gradient2()+
  labs(x = "Expression Value Added", y = "AUC Value Added", fill = "% q<0.05")
ggsave(a0.05.tile.plot, file = paste0("WGV_drugSEA_sim_q0.05_N50_",Sys.Date(),".pdf"))

q0.25.tile.plot <- ggplot(plot.data, aes(x = expr.value, y = AUC.value, fill = percent.q0.25)) +
  geom_tile() + ng.theme + theme_light() + scale_fill_gradient2() + 
  labs(x = "Expression Value Added", y = "AUC Value Added", fill = "% q<0.25")
ggsave(q0.25.tile.plot, file = paste0("WGV_drugSEA_sim_q0.25_N50_",Sys.Date(),".pdf"))

NES.tile.plot <- ggplot(plot.data, aes(x = expr.value, y = AUC.value, fill = Average.NES)) +
  geom_tile() + ng.theme + theme_light() + scale_fill_gradient2() +
  labs(x = "Expression Value Added", y = "AUC Value Added", fill = "Mean NES")
ggsave(NES.tile.plot, file = paste0("WGV_drugSEA_sim_NES_N50_",Sys.Date(),".pdf"))
