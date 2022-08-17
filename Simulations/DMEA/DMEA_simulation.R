# Authors: Belinda B. Garana (BG), James H. Joly (JJ)
# JJ sent 12-17-20; BG adapted for DMEA; last edit: BG 2022-06-20
### We want to know the true positive / false negative rates of enrichment
### So we will create synthetic data with enrichments for a drug set
### We will vary drug rank metric for a set of X drugs all in one direction, then run DMEA

### How to build this?

### Step 1: Generate drug rank metrics that have a normal distribution for each cell line, then perturb 1 drug set of X drugs
### Step 2: Run DMEA for each different variant (+/- Y drug AUC)
### Step 3: Re-run for 50 replicates
### Step 4: Repeat with different values of X (5, 10, 20, 40)

rm(list=ls(all=TRUE))
if(!require(devtools)){install.packages("dev.tools")}
if(!require(DMEA)){devtools::install_github('BelindaBGarana/DMEA')}
library(DMEA);library(dplyr);library(GSA);library(reshape2);library(data.table);library(ggplot2);

### Step 0: Prep rows, columns, values
#Import PRISM data for drug list
prism <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/main/Inputs/PRISM_drug_mean_AUC_6-23-21.csv",header=T) #481 cell lines
prism$X <- NULL

#Import drug info
drug.info <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/main/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv",header=T) #481 cell lines
drug.info$X <- NULL
colnames(drug.info)[colnames(drug.info)=="name"] <- "Drug"
drug.info$Drug <- gsub("[[:punct:][:blank:]]",".",drug.info$Drug)
drug.moa <- na.omit(distinct(drug.info[,c("Drug","moa")]))
Drug <- unique(drug.moa$Drug)
moa.drugs <- data.frame(Drug)
Drug <- unique(colnames(prism[,2:ncol(prism)]))
prism.drugs <- data.frame(Drug)
overlap.drugs <- merge(moa.drugs,prism.drugs,by="Drug")

#reduce PRISM data to drugs with moa
prism <- prism %>% select("CCLE_ID",overlap.drugs$Drug)
drug.names <- unique(colnames(prism[,2:ncol(prism)]))
num.drugs <- length(drug.names) #1351

# Make synthetic cell line name - only need 1!
synthetic.cell.names <- "Cell_No_1"

# How far do we want to vary values here?
min.vary <- -1
max.vary <- 1
values.to.vary <- seq(from = min.vary, to = max.vary, by = 0.25)

# Import MOA set information
gmt <- GSA.read.gmt(file="https://raw.github.com/BelindaBGarana/DMEA/main/Inputs/MOA_gmt_file_n6_no_special_chars.gmt")

path.sim <- paste0("Sim_results_",Sys.Date())
dir.create(path.sim)
# Change path to where you want files saved
setwd(path.sim)

all.synth.sim.results <- list()
all.sim.results <- list()
all.perc.sig <- list()
synth.drug.sets <- data.frame(N_drugs=integer(), synthetic.drug.set=character())
drug.set.sizes <- c(5,6,10,15,20,25,30,35,40,45,50)
for(i in 1:(length(drug.set.sizes))){
  synthetic.drug.set <- sample(drug.names,drug.set.sizes[i])
  synthetic.drug.df <- as.data.frame(synthetic.drug.set)
  synthetic.drug.df$N_drugs <- drug.set.sizes[i]
  synth.drug.sets <- rbind(synth.drug.sets, synthetic.drug.df)
 
  # make new gmt with synthetic drug set
  new.drug.sets <- gmt$genesets
  new.drug.set.names <- gmt$geneset.names
  new.drug.set.descr <- gmt$geneset.descriptions
  new.drug.set.names[[length(gmt$genesets)+1]] <- paste0(drug.set.sizes[i],"_random_drugs")
  new.drug.set.descr[[length(gmt$genesets)+1]] <- paste0(drug.set.sizes[i],"_random_drugs")
  new.drug.sets[[length(gmt$genesets)+1]] <- synthetic.drug.set
  new.gmt <- list(genesets=new.drug.sets, geneset.names=new.drug.set.names, geneset.descriptions=new.drug.set.descr)
  
  synth.sim.results <- list()
  sim.results <- list()
  for(N in 1:50){
    ### Step 3: Generate drug AUC that have a normal distribution for each cell line, then perturb synthetic drug set
    rank.matrix <- matrix(data = NA, nrow = length(drug.names), ncol = length(synthetic.cell.names)+1)
    row.names(rank.matrix) <- drug.names
    colnames(rank.matrix) <- c("Drug",synthetic.cell.names)
    rank.matrix[,"Drug"] <- drug.names
    rank.matrix <- as.data.frame(rank.matrix)
    rank.matrix[,synthetic.cell.names] <- rnorm(num.drugs, mean=0, sd=0.5)
    
    synthetic.results <- as.data.frame(values.to.vary)
    synthetic.results[,c("ES","NES","p","pos_at_max","q")] <- NA
    results <- list()
    for(j in 1:length(values.to.vary)){
      # perturb synthetic drug set
      rank.matrix[synthetic.drug.set,synthetic.cell.names] <- rnorm(drug.set.sizes[i], mean=values.to.vary[j], sd=0.5)
      
      # run DMEA
      DMEA.result <- drugSEA(data = rank.matrix, gmt = new.gmt, rank.metric=synthetic.cell.names, FDR=0, min.per.set=5)$result
      synthetic.results[j, 2:6] <- DMEA.result[DMEA.result$Drug_set==paste0(drug.set.sizes[i],"_random_drugs"), 3:7]
      results[[as.character(values.to.vary[j])]] <- DMEA.result
    }
    synth.sim.results[[N]] <- synthetic.results 
    sim.results[[N]] <- rbindlist(results, use.names=TRUE, idcol="Value Added", fill=TRUE)
  }
  all.synth.sim.results[[as.character(drug.set.sizes[i])]] <- rbindlist(synth.sim.results, use.names=TRUE, idcol="Simulation Number", fill=TRUE)
  all.sim.results[[as.character(drug.set.sizes[i])]] <- rbindlist(sim.results, use.names=TRUE, idcol="Simulation Number", fill=TRUE)
  write.csv(all.sim.results[[as.character(drug.set.sizes[i])]],file=paste0(drug.set.sizes[i],"_drugs_",N,"_replicates_all_drugSEA_results_",Sys.Date(),".csv"))
  write.csv(all.synth.sim.results[[as.character(drug.set.sizes[i])]],file=paste0(drug.set.sizes[i],"_drugs_",N,"_replicates_Synthetic_Drug_Set_drugSEA_results_",Sys.Date(),".csv"))
  all.synth.sim.results.df <- all.synth.sim.results[[as.character(drug.set.sizes[i])]]
  
  perc.sig <- as.data.frame(values.to.vary)
  perc.sig[,c("N.replicates","avg.NES","sd.NES","perc.q.0.25","perc.q.0.05","perc.q.0.25.pos","perc.q.0.25.neg","perc.q.0.05.pos","perc.q.0.05.neg")] <- NA
  for(k in 1:length(values.to.vary)){
    temp.data <- all.synth.sim.results.df[all.synth.sim.results.df$values.to.vary==values.to.vary[k],]
    
    sig.temp.data <- temp.data[temp.data$p<0.05 & temp.data$q<0.25,]
    pos.sig.temp.data <- sig.temp.data[sig.temp.data$NES>0,]
    neg.sig.temp.data <- sig.temp.data[sig.temp.data$NES<0,]
    
    more.sig.temp.data <- sig.temp.data[sig.temp.data$q<0.05,]
    pos.more.sig.temp.data <- more.sig.temp.data[more.sig.temp.data$NES>0,]
    neg.more.sig.temp.data <- more.sig.temp.data[more.sig.temp.data$NES<0,]
    
    perc.sig$N.replicates[k] <- nrow(temp.data)
    perc.sig$avg.NES[k] <- mean(as.numeric(temp.data$NES))
    perc.sig$sd.NES[k] <- sd(as.numeric(temp.data$NES))
    perc.sig$perc.q.0.25[k] <- 100*nrow(sig.temp.data)/nrow(temp.data)
    perc.sig$perc.q.0.05[k] <- 100*nrow(more.sig.temp.data)/nrow(temp.data)
    perc.sig$perc.q.0.25.pos[k] <- 100*nrow(pos.sig.temp.data)/nrow(temp.data)
    perc.sig$perc.q.0.05.pos[k] <- 100*nrow(pos.more.sig.temp.data)/nrow(temp.data)
    perc.sig$perc.q.0.25.neg[k] <- 100*nrow(neg.sig.temp.data)/nrow(temp.data)
    perc.sig$perc.q.0.05.neg[k] <- 100*nrow(neg.more.sig.temp.data)/nrow(temp.data)
  }
  write.csv(perc.sig,file=paste0(drug.set.sizes[i],"_drugs_",N,"_replicates_percent_significant_synthetic_drug_set_enrichments.csv"))
  all.perc.sig[[as.character(drug.set.sizes[i])]] <- perc.sig
}
complete.synth.sim.results <- rbindlist(all.synth.sim.results, use.names=TRUE, idcol="Number of Drugs in Set", fill=TRUE)
complete.sim.results <- rbindlist(all.sim.results, use.names=TRUE, idcol="Number of Drugs in Set", fill=TRUE)
write.csv(complete.sim.results,file=paste0(N,"_replicates_all_drugSEA_results_simulation_",Sys.Date(),".csv"), row.names=FALSE)
write.csv(complete.synth.sim.results,file=paste0(N,"_replicates_synthetic_Drug_Set_drugSEA_results_",Sys.Date(),".csv"), row.names=FALSE)
complete.perc.sig <- rbindlist(all.perc.sig, use.names=TRUE, idcol="Number of Drugs in Set", fill=TRUE)
write.csv(complete.perc.sig,file=paste0(N,"_replicates_percent_significant_synthetic_drug_set_enrichments_",Sys.Date(),".csv"), row.names=FALSE)

## Generate tile plots
# prepare collapsed results
complete.perc.sig$N_drugs <- as.numeric(complete.perc.sig$`Number of Drugs in Set`)
complete.perc.sig$values.to.vary <- as.numeric(complete.perc.sig$values.to.vary)
complete.perc.sig$perc.q.0.05 <- as.numeric(complete.perc.sig$perc.q.0.05)

# load theme for plot
ng.theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.border = element_rect(fill=NA), panel.background = element_blank(),
                  axis.line = element_line(colour = "black"), axis.text.x = element_text(colour = "black"),
                  axis.text.y = element_text(colour = "black"), axis.ticks.x = element_line(colour="black"),
                  axis.ticks.y = element_line(colour="black"), legend.title = element_blank(),
                  axis.title.y = element_text(size=8, colour="black"))

# produce tile plots
q0.05.tile.plot <- ggplot(complete.perc.sig, aes(x = as.factor(N_drugs), y = values.to.vary, fill = perc.q.0.05)) +
  geom_tile(height=0.25) + ng.theme + theme_light() + scale_fill_gradient2() +
  scale_x_discrete(breaks = c(5,6,10,15,20,25,30,35,40,45,50,55)) + 
  labs(x="Number of Drugs in Set", y = "Value Added", fill = "% q<0.05")
ggsave(q0.05.tile.plot, file = paste0("drugSEA_sim_q0.05_N50_",Sys.Date(),".pdf"))

q0.25.tile.plot <- ggplot(complete.perc.sig, aes(x = as.factor(N_drugs), y = values.to.vary, fill = perc.q.0.25)) +
  geom_tile(height=0.25) + ng.theme + theme_light() + scale_fill_gradient2() +
  scale_x_discrete(breaks = c(5,6,10,15,20,25,30,35,40,45,50,55)) + 
  labs(x="Number of Drugs in Set", y = "Value Added", fill = "% q<0.25")
ggsave(q0.25.tile.plot, file = paste0("drugSEA_sim_q0.25_N50_",Sys.Date(),".pdf"))

NES.tile.plot <- ggplot(complete.perc.sig, aes(x = as.factor(N_drugs), y = values.to.vary, fill=avg.NES)) +
  geom_tile(height=0.25) + ng.theme + theme_light() + scale_fill_gradient2() +
  scale_x_discrete(breaks = c(5,6,10,15,20,25,30,35,40,45,50,55)) +  
  labs(x="Number of Drugs in Set", y = "Value Added", fill = "Mean NES")
ggsave(NES.tile.plot, file = paste0("drugSEA_sim_NES_N50_",Sys.Date(),".pdf"))
