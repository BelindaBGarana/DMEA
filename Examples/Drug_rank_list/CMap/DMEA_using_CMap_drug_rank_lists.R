#Using CMap drug rank lists for DMEA
#Author: Belinda B. Garana (BG), Date: 2022-05-12; last edit: BG 2022-07-02
rm(list=ls(all=TRUE))

library(plyr);library(dplyr);library(ggplot2);
library(sjmisc);library(DMEA);

illegal.chars <- c("#","<","%",">","!","`","&","'","=","}","/",":","@") # source: https://www.mtu.edu/umc/services/websites/writing/characters-avoid/; would be nice to protect against " and \ too
illegal.chars.need.brackets <- c("$","+","*","|","{","?") # these characters need brackets for gsub, but not for str_contains

## To replicate: change path and run through whole script
# If recreating from CMap direct output, create folders with matching CMap output gct file for each query type and case
# and uncomment lines 41-58 and comment out line 59 before running whole script
path <- "/Users/belindagarana/Documents/Graham\ Lab/drugSEA/CMap/"
cases.L1000 <- c("Delfarah_et_al_HMEC_senescence",
                 "GSE14003_Proteasome_inhibitor_bortezomib_treated_JEKO1_cells_10H_vs_untreated",
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
  setwd(path)
  dir.create(query.types[k])
  for(i in 1:length(cases)){
    setwd(paste0(path,query.types[k],"/"))
    dir.create(cases[i])
    setwd(cases[i])
    # need CMap output to replicate lines 41-58 (could not upload to GitHub without interfering with DMEA package)
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
    # write.csv(input.df, file="DMEA_input.csv", row.names=FALSE)
    input.df <- read.csv(file=paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Examples/Drug_rank_list/CMap/",query.types[k],"/",cases[i],"/DMEA_input.csv"))
    
    # run DMEA
    DMEA.results <- drugSEA(input.df, drug="pert_iname", rank.metric=rank.metric)
    saveRDS(DMEA.results$gmt, file="DMEA_gmt.rds")
    write.csv(DMEA.results$result, file="DMEA_results.csv", row.names=FALSE)
    ggsave(file="DMEA_volcano_plot.pdf", DMEA.results$volcano.plot)
    if(length(DMEA.results$mtn.plots)>0){
      for(j in 1:length(DMEA.results$mtn.plots)){
        # if moa name contains illegal filename characters, replace these characters with "-"
        moa.name <- selected.moa
        if(str_contains(moa.name, illegal.chars, logic="or")){
          for(m in 1:length(illegal.chars)){
            moa.name <- gsub(illegal.chars[m], "-", moa.name)
          }
        }else if(str_contains(moa.name, illegal.chars.need.brackets, logic="or")){
          for(m in 1:length(illegal.chars.need.brackets)){
            moa.name <- gsub(paste0("[",illegal.chars.need.brackets[m],"]"), "-", moa.name)
          }
        }
        ggsave(file=paste0("DMEA_mountain_plot_",moa.name,".pdf"),DMEA.results$mtn.plots[[j]])
      }
    }
  }
}
