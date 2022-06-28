#Using CMap drug rank lists for DMEA
#Author: Belinda B. Garana (BG), Date: 2022-05-12; last edit: BG 2022-06-25
rm(list=ls(all=TRUE))

library(plyr);library(dplyr);library(ggplot2);
library(sjmisc);library(DMEA);

##### To replicate: change path and create folders with CMap query output for each query type and case #####
path <- "/Users/belindagarana/Documents/Graham\ Lab/drugSEA/CMap/"
cases.L1000 <- c("GSE14003_Proteasome_inhibitor_bortezomib_treated_JEKO1_cells_10H_vs_untreated",
                 "GSE28896_Glucocorticoid_agonist_dexamethasone_treated_CD34_cells_24H_vs_untreated",
                 "GSE32547_Pitavastatin_treated_HUVEC_cells_1_uM_at_4H_vs_DMSO_treated",
                 "GSE33643_PI3K_MTOR_inhibitor_BEZ235_treated_A2058_cells_3_doses_at_24H_vs_DMSO_treated",
                 "GSE35230_MEK_inhibitor_GSK212_treated_A375_clones_30nM_at_24H_vs_DMSO_treated")
cases.PRISM <- c("Cell_lines_sensitive_to_HMGCR_inhibitor_lovastatin",
                 "Cell_lines_with_EGFR_activating_mutation_ELREA746DEL",
                 "Cell_lines_with_high_gene_expression_of_PDGFRA")
query.types <- c("L1000", "PRISM")
for(k in 1:length(query.types)){
  # set rank.metric and cases for each query type
  if(query.types[k]=="L1000"){
    cases <- cases.L1000
    rank.metric <- "norm_cs"
  }else if(query.types[k]=="PRISM"){
    cases <- cases.PRISM
    rank.metric <- "query"
  }
  for(i in 1:length(cases)){
    setwd(paste0(path,query.types[k],"/",cases[i],"/"))
    # need CMap output to replicate lines 31-48 (could not upload to GitHub without interfering with DMEA package)
    # # import query result
    # if(query.types[k]=="L1000"){
    #   setwd(paste0(path,query.types[k],"/",cases[i],"/arfs/TAG/"))
    #   input.df <- read.delim(file="query_result.gct", skip = 2)
    #   input.df <- input.df[input.df$qc_pass==1,] 
    # }else if(query.types[k]=="PRISM"){
    #   input.df <- read.delim(file="ncs.gct", skip = 2)
    # }
    # # remove any drugs without moa annotations, duplicates, or NA
    # input.df <- input.df[input.df$pert_iname!="na",] # removes empty first row in CMap query output
    # input.df <- input.df[input.df$moa!="-666",] # remove drugs without known moa (moa="-666")
    # input.df[,c(rank.metric)] <- as.numeric(input.df[,c(rank.metric)])
    # if(query.types[k]=="L1000"){
    #   input.df <- ddply(input.df, .(pert_iname, moa), summarize, norm_cs = mean(norm_cs, na.rm=TRUE)) # average across cell lines
    #   setwd(paste0(path,query.types[k],"/",cases[i],"/"))
    # }
    # input.df <- distinct(na.omit(input.df[!duplicated(input.df[,c("pert_iname")]),c("pert_iname", rank.metric, "moa")])) # remove any NA, duplicates
    # write.csv(input.df, file="DMEA_input.csv")
    input.df <- read.csv(file="DMEA_input.csv")
    
    # run DMEA
    DMEA.results <- drugSEA(input.df, drug="pert_iname", rank.metric=rank.metric)
    saveRDS(DMEA.results$gmt, file="DMEA_gmt.rds")
    write.csv(DMEA.results$result, file="DMEA_results.csv")
    ggsave(file="DMEA_volcano_plot.pdf", DMEA.results$volcano.plot)
    if(length(DMEA.results$mtn.plots)>0){
      for(j in 1:length(DMEA.results$mtn.plots)){
        if(str_contains(names(DMEA.results$mtn.plots)[[j]], "/")){
          moa.name <- gsub("/", "-", names(DMEA.results$mtn.plots)[[j]])
          ggsave(file=paste0("DMEA_mountain_plot_",moa.name,".pdf"),DMEA.results$mtn.plots[[j]])
        }else{
          ggsave(file=paste0("DMEA_mountain_plot_",names(DMEA.results$mtn.plots)[[j]],".pdf"),DMEA.results$mtn.plots[[j]])
        }
      }
    }
  }
}
