# DMEA using gene signatures
# Author: Belinda B. Garana (BG), Date: 2022-05-05; last edit: BG 2022-06-30

# How to recreate this? Change path in line 12 and press run through the whole script
# Steps for each dataset:
# a - generate gene signature using limma::eBayes
# b - run DMEA::DMEA without any samples used to define the gene signature

rm(list=ls(all=TRUE))
library(plyr);library(dplyr);library(GSA);library(DMEA);library(GEOquery);library(limma);

path <- "/Users/belindagarana/Documents/Graham\ Lab/DMEA_2022-05-05/Gene_signatures/"

##### Part 1: Import data (RNAseq, PRISM drug AUC, drug moa gmt) #####
illegal.chars <- c("#","<","%",">","!","`","&","'","=","}","/",":","@") # source: https://www.mtu.edu/umc/services/websites/writing/characters-avoid/; would be nice to protect against " and \ too
illegal.chars.need.brackets <- c("$","+","*","|","{","?") # these characters need brackets for gsub, but not for str_contains

#### cell line info (CCLE 19Q4)
cell.line.info <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/main/Inputs/CCLE_sample_info.csv")
cell.line.info$X <- NULL

#### drug info (PRISM)
drug.info <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/main/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv")
colnames(drug.info)[colnames(drug.info)=="name"] <- "Drug"
drug.info$Drug <- gsub("[[:punct:][:blank:]]",".",drug.info$Drug)
drug.moa <- na.omit(distinct(drug.info[,c("Drug","moa")]))

#### RNAseq (CCLE 19Q4)
download.file("https://raw.github.com/BelindaBGarana/DMEA/main/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin", 
              destfile = paste0(getwd(),"/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"))
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
download.file("https://raw.github.com/BelindaBGarana/DMEA/main/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin", 
              destfile = paste0(getwd(),"/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"))
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
RNA.df <- rbind(RNA.first200, RNA.rest)

#### PRISM drug AUC
AUC.df <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/main/Inputs/PRISM_drug_mean_AUC_6-23-21.csv") #481 cell lines
AUC.df$X <- NULL

### only use drugs with moa annotations
rownames(AUC.df) <- AUC.df$CCLE_ID
AUC.df <- AUC.df[,colnames(AUC.df) %in% drug.moa$Drug]
AUC.df$CCLE_ID <- rownames(AUC.df)
AUC.df <- AUC.df[,c("CCLE_ID",colnames(AUC.df)[1:(ncol(AUC.df)-1)])] # for DMEA, the sample names must be in the first column for both the drug sensitivity and gene expression dataframes

#### drug moa sets
gmt <- GSA.read.gmt(file="https://raw.github.com/BelindaBGarana/DMEA/main/Inputs/MOA_gmt_file_n6_no_special_chars.gmt")

##### Part 2: Coldren et al (EGFR inhibitor) #####
dir.create("Coldren_et_al_NSCLC_sensitive_vs_resistant_to_gefitinib")
setwd(paste0(path,"Coldren_et_al_NSCLC_sensitive_vs_resistant_to_gefitinib"))
#### Part 2a - generate gene signature using limma::eBayes ####
Sys.setenv("VROOM_CONNECTION_SIZE"=300000) # increases size of the connection buffer

### Identify sensitive (SENS) or resistant (RES) cells based on Coldren et al, 2006 paper
SENS.cells <- c("NCIH358","NCIH322C","CALU3","NCIH1334","NCIH1648",
                "HCC827","HCC78","NCIH2126","HCC193","HCC95",
                "NCIH3255","HCC4006","NCIH2009","NCIH1650","NCIH820","NCIH1975") # not in CCLE RNAseq: CALU3
RES.cells <- c("NCIH125","NCIH1703","A549","NCIH157","NCIH460","NCIH520",
               "HCC44","HCC15","NCIH157","NCIH1299","HCC2279")
SENS.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% SENS.cells,] # 5 cell lines: NCIH1650, HCC95, NCIH1975, NCIH1648, NCIH2126
SENS.CCLE_ID <- SENS.cell.info$CCLE_ID
RES.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% RES.cells,] # 7 cell lines: NCIH520, NCIH460, NCIH1299, HCC44, A549, NCIH1703, HCC15
RES.CCLE_ID <- RES.cell.info$CCLE_ID

### Get data for SENS and RES cells
SENS.data <- merge(SENS.cell.info, RNA.df, by="CCLE_ID")
SENS.data$stripped_cell_line_name <- NULL
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL
RES.data <- merge(RES.cell.info, RNA.df, by="CCLE_ID")
RES.data$stripped_cell_line_name <- NULL
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL

### assign samples to groups 
data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
SENS.RES.data <- cbind(data.SENS,data.RES)
gset <- ExpressionSet(assayData = SENS.RES.data)
sml <- c(0,0,0,0,0,1,1,1,1,1,1,1)

### set up design matrix
gs <- factor(sml)
groups <- make.names(c("test","test2"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

### set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

### compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
tT$Gene <- rownames(tT) # 19,144 genes
tT.noblanks <- tT[tT$Gene!="",] # 19,144 genes
colnames(tT.noblanks)[1] <- "Log2FC"
weights <- na.omit(tT.noblanks[tT.noblanks$Gene %in% colnames(RNA.df.nooverlap),]) # 19,111 genes
#write.csv(weights, file="Full_gene_signature.csv", row.names=FALSE)
filtered.weights.q0.05 <- weights[weights$adj.P.Val<0.05,] # 31 genes
avg.filtered.weights.q0.05 <- ddply(filtered.weights.q0.05, .(Gene), summarize,
                                    Log2FC = mean(Log2FC, na.rm=TRUE),
                                    sd.Log2FC = sd(Log2FC, na.rm=TRUE)) # prevent duplicate gene names; 31 genes
if(nrow(avg.filtered.weights.q0.05) < nrow(filtered.weights.q0.05)){
  filtered.weights <- avg.filtered.weights.q0.05
}else{filtered.weights <- filtered.weights.q0.05}
write.csv(filtered.weights, file="Filtered_gene_signature_no_duplicates.csv", row.names=FALSE)

#### Part 2b - run DMEA::DMEA without any samples used to define the gene signature ####
### remove cell lines used to define gene weights
overlap <- c(SENS.CCLE_ID,RES.CCLE_ID)
RNA.df.nooverlap <- RNA.df[!RNA.df$CCLE_ID %in% overlap, ] # 602 cell lines out of 614
rownames(RNA.df.nooverlap) <- RNA.df.nooverlap$CCLE_ID

### run DMEA with top 500 genes based on abs(Log2FC)
if(nrow(filtered.weights) > 500){
  filtered.weights <- filtered.weights %>% slice_max(abs(Log2FC), n=500)
}
DMEA.results <- DMEA(AUC.df, gmt, RNA.df.nooverlap, filtered.weights, gene.names = "Gene", weight.values = "Log2FC")
write.csv(DMEA.results$result, file="DMEA_results.csv", row.names=FALSE)
ggsave(file="DMEA_volcano_plot.pdf", DMEA.results$volcano.plot, device="pdf")
write.csv(DMEA.output$WV.scores, "DMEA_WGV_scores.csv", row.names = FALSE)
write.csv(DMEA.output$corr.result, "DMEA_correlation_results.csv", row.names = FALSE)
ggsave("DMEA_correlation_plots.pdf", DMEA.output$corr.scatter.plots, device="pdf")
if(length(DMEA.results$mtn.plots)>0){
  for(j in 1:length(DMEA.results$mtn.plots)){
    # if moa name contains illegal filename characters, replace these characters with "-"
    moa.name <- selected.moa
    if(str_contains(moa.name, illegal.chars, logic="or")){
      for(k in 1:length(illegal.chars)){
        moa.name <- gsub(illegal.chars[k], "-", moa.name)
      }
    }else if(str_contains(moa.name, illegal.chars.need.brackets, logic="or")){
      for(k in 1:length(illegal.chars.need.brackets)){
        moa.name <- gsub(paste0("[",illegal.chars.need.brackets[k],"]"), "-", moa.name)
      }
    }
    ggsave(file=paste0("DMEA_mountain_plot_",moa.name,".pdf"),DMEA.results$mtn.plots[[j]], device="pdf")
  }
}

##### Part 3: GSE31625 (EGFR inhibitor) #####
dir.create("GSE31625_NSCLC_sensitive_vs_resistant_to_erlotinib")
setwd(paste0(path,"GSE31625_NSCLC_sensitive_vs_resistant_to_erlotinib"))

#### Part 3a - generate gene signature using limma::eBayes ####
Sys.setenv("VROOM_CONNECTION_SIZE"=300000) # increases size of the connection buffer

# load series and platform data from GEO
gset <- getGEO("GSE31625", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("X11111111XXXXXXXXX0000000XXXXXXXX0000X0000001111") # 0 is phenotype of numerator in FC; 1 of denominator; X: excluded
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }
SENS.expr <- as.data.frame(exprs(gset[,which(sml==0)]))
RES.expr <- as.data.frame(exprs(gset[,which(sml==1)]))

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("test","test2"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2$genes))

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title")) # 22,215 genes
tT.noblanks <- tT[tT$Gene.symbol!="",] # 21,146 genes
colnames(tT.noblanks)[6] <- "Log2FC"
colnames(tT.noblanks)[7] <- "Gene"
weights <- na.omit(tT.noblanks[tT.noblanks$Gene %in% colnames(RNA.df.nooverlap),]) # 18,980 genes
#write.csv(weights, file="Full_gene_signature.csv", row.names=FALSE)
filtered.weights.q0.05 <- weights[weights$adj.P.Val<0.05,] # 4,266 genes
avg.filtered.weights.q0.05 <- ddply(filtered.weights.q0.05, .(Gene), summarize,
                                    Log2FC = mean(Log2FC, na.rm=TRUE),
                                    sd.Log2FC = sd(Log2FC, na.rm=TRUE)) #prevent duplicate gene names; 3,336 genes
if(nrow(avg.filtered.weights.q0.05) < nrow(filtered.weights.q0.05)){
  filtered.weights <- avg.filtered.weights.q0.05
}else{filtered.weights <- filtered.weights.q0.05}
write.csv(filtered.weights, file="Filtered_gene_signature_no_duplicates.csv", row.names=FALSE)

#### Part 3b - run DMEA::DMEA without any samples used to define the gene signature ####
### remove cell lines used to define gene signature
SENS.cells <- c("NCIH1650") # this is the only sensitive cell line from this dataset in the CCLE RNAseq version 19Q4
RES.cells <- c("A549") # this is the only resistant cell line from this dataset in the CCLE RNAseq version 19Q4
overlap <- c(SENS.cells,RES.cells)
SENS.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% SENS.cells,]
SENS.CCLE_ID <- SENS.cell.info$CCLE_ID
RES.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% RES.cells,]
RES.CCLE_ID <- RES.cell.info$CCLE_ID
overlap <- c(SENS.CCLE_ID,RES.CCLE_ID)

SENS.data <- merge(SENS.cell.info, RNA.df, by="CCLE_ID")
SENS.data$stripped_cell_line_name <- NULL
RES.data <- merge(RES.cell.info, RNA.df, by="CCLE_ID")
RES.data$stripped_cell_line_name <- NULL

overlap <- c(SENS.CCLE_ID,RES.CCLE_ID)
RNA.df.nooverlap <- RNA.df[!RNA.df$CCLE_ID %in% overlap, ]
rownames(RNA.df.nooverlap) <- RNA.df.nooverlap$CCLE_ID

### run DMEA with top 500 genes based on abs(Log2FC)
if(nrow(filtered.weights) > 500){
  filtered.weights <- filtered.weights %>% slice_max(abs(Log2FC), n=500)
}
DMEA.results <- DMEA(AUC.df, gmt, RNA.df.nooverlap, filtered.weights, gene.names = "Gene", weight.values = "Log2FC")
ggsave(file="DMEA_volcano_plot.pdf", DMEA.results$volcano.plot, device="pdf")
write.csv(DMEA.output$WV.scores, "DMEA_WGV_scores.csv", row.names = FALSE)
write.csv(DMEA.output$corr.result, "DMEA_correlation_results.csv", row.names = FALSE)
ggsave("DMEA_correlation_plots.pdf", DMEA.output$corr.scatter.plots, device="pdf")
write.csv(DMEA.results$result, file="DMEA_results.csv", row.names=FALSE)
if(length(DMEA.results$mtn.plots)>0){
  for(j in 1:length(DMEA.results$mtn.plots)){
    # if moa name contains illegal filename characters, replace these characters with "-"
    moa.name <- selected.moa
    if(str_contains(moa.name, illegal.chars, logic="or")){
      for(k in 1:length(illegal.chars)){
        moa.name <- gsub(illegal.chars[k], "-", moa.name)
      }
    }else if(str_contains(moa.name, illegal.chars.need.brackets, logic="or")){
      for(k in 1:length(illegal.chars.need.brackets)){
        moa.name <- gsub(paste0("[",illegal.chars.need.brackets[k],"]"), "-", moa.name)
      }
    }
    ggsave(file=paste0("DMEA_mountain_plot_",moa.name,".pdf"),DMEA.results$mtn.plots[[j]], device="pdf")
  }
}

##### Part 4: GSE12790 (EGFR inhibitor) #####
dir.create("GSE12790_BRCA_sensitive_vs_resistant_to_erlotinib")
setwd(paste0(path,"GSE12790_BRCA_sensitive_vs_resistant_to_erlotinib"))

#### Part 4a - generate gene signature using limma::eBayes ####
Sys.setenv("VROOM_CONNECTION_SIZE"=300000) # increases size of the connection buffer

# load series and platform data from GEO
gset <- getGEO("GSE12790", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX1XX",
               "X0X11X1XX1X1X1XX1XX1X1110X1XX1XXX11X1XX0XX1XXXXX") # 0 is phenotype of numerator in FC; 1 of denominator; X: excluded
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("test","test2"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2$genes))

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title")) # 54,675 genes
tT.noblanks <- tT[tT$Gene.symbol!="",] #45,118 genes
colnames(tT.noblanks)[6] <- "Log2FC"
colnames(tT.noblanks)[7] <- "Gene"
weights <- na.omit(tT.noblanks[tT.noblanks$Gene %in% colnames(RNA.df.nooverlap),]) # 37,107 genes
#write.csv(weights, file="Full_gene_signature.csv", row.names=FALSE)
filtered.weights.q0.05 <- weights[weights$adj.P.Val<0.05,] #54
avg.filtered.weights.q0.05 <- ddply(filtered.weights.q0.05, .(Gene), summarize,
                                    Log2FC = mean(Log2FC, na.rm=TRUE),
                                    sd.Log2FC = sd(Log2FC, na.rm=TRUE)) #prevent duplicate gene names; 49 genes
if(nrow(avg.filtered.weights.q0.05) < nrow(filtered.weights.q0.05)){
  filtered.weights <- avg.filtered.weights.q0.05
}else{filtered.weights <- filtered.weights.q0.05}
write.csv(filtered.weights, file="Filtered_gene_signature_no_duplicates.csv", row.names=FALSE)

#### Part 4b - run DMEA::DMEA without any samples used to define the gene signature ####
### remove cell lines used to define gene signature
SENS.cells <- c("HCC1806") #HDQP1, CAL851 not in CCLE RNAseq
RES.cells <- c("CAL51", "CAL120", "MDAMB231", "HCC1954",
               "MDAMB453", "HCC1428", "T47D", "ZR751",
               "KPL1", "MDAMB415", "CAMA1", "MCF7") #BT20, HCC1569, EFM192A, BT474, BT483, HCC1500 not in CCLE RNAseq
overlap <- c(SENS.cells,RES.cells)
SENS.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% SENS.cells,] # 1 cell line
SENS.CCLE_ID <- SENS.cell.info$CCLE_ID
RES.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% RES.cells,] # 12 cell lines
RES.CCLE_ID <- RES.cell.info$CCLE_ID
overlap <- c(SENS.CCLE_ID,RES.CCLE_ID)

RNA.df.nooverlap <- RNA.df[!RNA.df$CCLE_ID %in% overlap, ] #601 cell lines out of 614
rownames(RNA.df.nooverlap) <- RNA.df.nooverlap$CCLE_ID

### run DMEA with top 500 genes based on abs(Log2FC)
if(nrow(filtered.weights) > 500){
  filtered.weights <- filtered.weights %>% slice_max(abs(Log2FC), n=500)
}
DMEA.results <- DMEA(AUC.df, gmt, RNA.df.nooverlap, filtered.weights, gene.names = "Gene", weight.values = "Log2FC")
write.csv(DMEA.results$result, file="DMEA_results.csv", row.names=FALSE)
ggsave(file="DMEA_volcano_plot.pdf", DMEA.results$volcano.plot, device="pdf")
write.csv(DMEA.output$WV.scores, "DMEA_WGV_scores.csv", row.names = FALSE)
write.csv(DMEA.output$corr.result, "DMEA_correlation_results.csv", row.names = FALSE)
ggsave("DMEA_correlation_plots.pdf", DMEA.output$corr.scatter.plots, device="pdf")
if(length(DMEA.results$mtn.plots)>0){
  for(j in 1:length(DMEA.results$mtn.plots)){
    # if moa name contains illegal filename characters, replace these characters with "-"
    moa.name <- selected.moa
    if(str_contains(moa.name, illegal.chars, logic="or")){
      for(k in 1:length(illegal.chars)){
        moa.name <- gsub(illegal.chars[k], "-", moa.name)
      }
    }else if(str_contains(moa.name, illegal.chars.need.brackets, logic="or")){
      for(k in 1:length(illegal.chars.need.brackets)){
        moa.name <- gsub(paste0("[",illegal.chars.need.brackets[k],"]"), "-", moa.name)
      }
    }
    ggsave(file=paste0("DMEA_mountain_plot_",moa.name,".pdf"),DMEA.results$mtn.plots[[j]], device="pdf")
  }
}

##### Part 5: GSE66539 (RAF inhibitor) #####
#### Part 5a - generate gene signature using limma::eBayes ####
dir.create("GSE66539_SKCM_sensitive_vs_resistant_to_dabrafenib")
setwd(paste0(path,"GSE66539_SKCM_sensitive_vs_resistant_to_dabrafenib"))

#### Using paired biopsy samples (3 patients before and after resistance) so no samples overlap with CCLE version 19Q4
### Using eBayes differential expression parameters
paired.ex <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/main/Examples/Gene_signature/GSE66539_SKCM_sensitive_vs_resistant_to_dabrafenib/Paired_log2FC_RAFi_GSE66539.csv")
paired.ex.biopsy <- paired.ex[,c("Gene.ID","M005","M019","M026")]
rownames(paired.ex.biopsy) <- paired.ex.biopsy$Gene.ID
paired.ex.biopsy.matrix <- as.matrix(paired.ex.biopsy[,-1])
gene.symbols <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/main/Examples/Gene_signature/GSE66539_SKCM_sensitive_vs_resistant_to_dabrafenib/ENSID_gene_names.csv")

paired.biopsy.gset <- ExpressionSet(assayData = paired.ex.biopsy.matrix)
fit <- lmFit(paired.biopsy.gset)  # fit linear model

# compute statistics and table of top significant genes
fit <- eBayes(fit, 0.01)
tT <- topTable(fit, adjust="fdr", sort.by="B", number=nrow(fit))
tT$ENSID <- rownames(tT)
tT.wsymbol <- merge(tT, gene.symbols, by="ENSID") #62,074 genes
tT.noblanks <- tT.wsymbol[tT.wsymbol$Gene!="",] #50,500 genes
colnames(tT.noblanks)[2] <- "Log2FC"
weights <- na.omit(tT.noblanks[tT.noblanks$Gene %in% colnames(RNA.df),]) #15,264 genes
#write.csv(weights, file="Full_gene_signature.csv", row.names=FALSE)
filtered.weights.q0.05 <- weights[weights$adj.P.Val<0.05,] #748 genes
avg.filtered.weights.q0.05 <- ddply(filtered.weights.q0.05, .(Gene), summarize,
                                    Log2FC = mean(Log2FC, na.rm=TRUE),
                                    sd.Log2FC = sd(Log2FC, na.rm=TRUE)) #prevent duplicate gene names; 748 genes
if(nrow(avg.filtered.weights.q0.05) < nrow(filtered.weights.q0.05)){
  filtered.weights <- avg.filtered.weights.q0.05
}else{filtered.weights <- filtered.weights.q0.05}
write.csv(filtered.weights, file="Filtered_gene_signature_no_duplicates.csv", row.names=FALSE)

#### Part 5b - run DMEA::DMEA without any samples used to define the gene signature ####
### run DMEA with top 500 genes based on abs(Log2FC)
if(nrow(filtered.weights) > 500){
  filtered.weights <- filtered.weights %>% slice_max(abs(Log2FC), n=500)
}
DMEA.results <- DMEA(AUC.df, gmt, RNA.df, filtered.weights, gene.names = "Gene", weight.values = "Log2FC")
write.csv(DMEA.results$result, file="DMEA_results.csv", row.names=FALSE)
ggsave(file="DMEA_volcano_plot.pdf", DMEA.results$volcano.plot, device="pdf")
write.csv(DMEA.output$WV.scores, "DMEA_WGV_scores.csv", row.names = FALSE)
write.csv(DMEA.output$corr.result, "DMEA_correlation_results.csv", row.names = FALSE)
ggsave("DMEA_correlation_plots.pdf", DMEA.output$corr.scatter.plots, device="pdf")
if(length(DMEA.results$mtn.plots)>0){
  for(j in 1:length(DMEA.results$mtn.plots)){
    # if moa name contains illegal filename characters, replace these characters with "-"
    moa.name <- selected.moa
    if(str_contains(moa.name, illegal.chars, logic="or")){
      for(k in 1:length(illegal.chars)){
        moa.name <- gsub(illegal.chars[k], "-", moa.name)
      }
    }else if(str_contains(moa.name, illegal.chars.need.brackets, logic="or")){
      for(k in 1:length(illegal.chars.need.brackets)){
        moa.name <- gsub(paste0("[",illegal.chars.need.brackets[k],"]"), "-", moa.name)
      }
    }
    ggsave(file=paste0("DMEA_mountain_plot_",moa.name,".pdf"),DMEA.results$mtn.plots[[j]], device="pdf")
  }
}


