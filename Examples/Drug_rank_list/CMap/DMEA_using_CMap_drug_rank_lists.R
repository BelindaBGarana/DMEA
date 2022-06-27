#Using CMap drug rank lists for DMEA
#Author: Belinda B. Garana (BG), Date: 2022-05-12; last edit: BG 2022-06-25
rm(list=ls(all=TRUE))

library(plyr);library(dplyr);library(ggplot2);
library(sjmisc);library(DMEA);

##### To replicate: change path to download location, unzip .gct files for L1000 queries, and create folders with CMap query output for each query type and case #####
path <- "/Users/belindagarana/Documents/Graham\ Lab/drugSEA/CMap/"
cases.L1000 <- c("CMap_L1000_query_GSE14003_Proteasome_inhibitor_bortezomib_treated_JEKO1_cells_10H_vs_untreated_accessed_2022-05-11",
                 "CMap_L1000_query_GSE28896_Glucocorticoid_agonist_dexamethasone_treated_CD34_cells_24H_vs_untreated_accessed_2022-05-11",
                 "CMap_L1000_query_GSE32547_Pitavastatin_treated_HUVEC_cells_1_uM_at_4H_vs_DMSO_treated_accessed_2022-05-11",
                 "CMap_L1000_query_GSE33643_PI3K_MTOR_inhibitor_BEZ235_treated_A2058_cells_3_doses_at_24H_vs_DMSO_treated_accessed_2022-05-11",
                 "CMap_L1000_query_GSE35230_MEK_inhibitor_GSK212_treated_A375_clones_30nM_at_24H_vs_DMSO_treated_accessed_2022-05-11")
cases.PRISM <- c("CMap_PRISM_query_Cell_lines_sensitive_to_HMGCR_inhibitor_lovastatin_accessed_2022-05-11",
                 "CMap_PRISM_query_Cell_lines_with_EGFR_activating_mutation_ELREA746DEL_accessed_2022-05-11",
                 "CMap_PRISM_query_Cell_lines_with_high_gene_expression_of_PDGFRA_accessed_2022-05-11")
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
    # import query result
    if(query.types[k]=="L1000"){
      setwd(paste0(path,query.types[k],"/",cases[i],"/arfs/TAG/"))
      input.df <- read.delim(file="query_result.gct", skip = 2)
      input.df <- input.df[input.df$qc_pass==1,] 
    }else if(query.types[k]=="PRISM"){
      setwd(paste0(path,query.types[k],"/",cases[i],"/"))
      input.df <- read.delim(file="ncs.gct", skip = 2)
    }
    # remove any drugs without moa annotations, duplicates, or NA
    input.df <- input.df[input.df$pert_iname!="na",] # removes empty first row in CMap query output
    input.df <- input.df[input.df$moa!="-666",] # remove drugs without known moa (moa="-666")
    input.df[,c(rank.metric)] <- as.numeric(input.df[,c(rank.metric)])
    if(query.types[k]=="L1000"){
      input.df <- ddply(input.df, .(pert_iname, moa), summarize, norm_cs = mean(norm_cs, na.rm=TRUE)) # average across cell lines
      setwd(paste0(path,query.types[k],"/",cases[i],"/"))
    }
    input.df <- distinct(na.omit(input.df[!duplicated(input.df[,c("pert_iname")]),c("pert_iname", rank.metric, "moa")])) # remove any NA, duplicates
    output.folder <- paste0("DMEA_analysis_",Sys.Date())
    dir.create(output.folder)
    setwd(paste0(path,query.types[k],"/",cases[i],"/",output.folder))
    write.csv(input.df, file=paste0("DMEA_input_",cases[i],".csv"))
    
    # run DMEA
    DMEA.results <- drugSEA(input.df, drug="pert_iname", rank.metric=rank.metric)
    saveRDS(DMEA.results$gmt, file=paste0("gmt_DMEA_",cases[i],".rds"))
    write.csv(DMEA.results$result, file=paste0("DMEA_results_",cases[i],".csv"))
    ggsave(file=paste0("DMEA_volcano_plot_",cases[i],".pdf"), DMEA.results$volcano.plot)
    if(length(DMEA.results$mtn.plots)>0){
      for(j in 1:length(DMEA.results$mtn.plots)){
        if(str_contains(names(DMEA.results$mtn.plots)[[j]], "/")){
          moa.name <- gsub("/", "-", names(DMEA.results$mtn.plots)[[j]])
          ggsave(file=paste0("DMEA_mountain_plot_",moa.name,"_",cases[i],".pdf"),DMEA.results$mtn.plots[[j]])
        }else{
          ggsave(file=paste0("DMEA_mountain_plot_",names(DMEA.results$mtn.plots)[[j]],"_",cases[i],".pdf"),DMEA.results$mtn.plots[[j]])
        }
      }
    }
  }
}
