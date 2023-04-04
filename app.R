#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Author: Belinda B. Garana; last edit: 2022-04-04

library(shiny);library(utils);library(GSA);library(DMEA);library(plyr);
library(dplyr);library(ggplot2);library(reshape2);library(gridExtra);library(sjmisc);
library(stats);library(foreach);library(cowplot);library(aplot);library(patchwork);library(ggrepel);
library(devtools);library(usethis);library(iterators);
library(tidyselect);library(parallel);library(testthat);library(snow);library(doSNOW);
library(BiocManager)
options(repos = BiocManager::repositories())
library(qvalue);

# set limit for upload file size
MB.limit <- 180

# set illegal filename characters
illegal.chars <- c("#","<","%",">","!","`","&","'","=","}","/",":","@") # source: https://www.mtu.edu/umc/services/websites/writing/characters-avoid/; would be nice to protect against " and \ too
illegal.chars.need.brackets <- c("$","+","*","|","{","?") # for these chars, gsub needs brackets but str_contains can't have brackets

# define functions for later use
as.filename <- function(moa.name){
  # if moa.name contains illegal file name characters, replace these characters with "-" or "_"
  if(sjmisc::str_contains(moa.name, illegal.chars, logic="or")){
    for(j in 1:length(illegal.chars)){
      moa.name <- gsub(illegal.chars[j], "-", moa.name)
    }
  }
  if(sjmisc::str_contains(moa.name, illegal.chars.need.brackets, logic="or")){
    for(j in 1:length(illegal.chars.need.brackets)){
      moa.name <- gsub(paste0("[",illegal.chars.need.brackets[j],"]"), "-", moa.name)
    }
  }
  if(sjmisc::str_contains(moa.name, " ")){
    moa.name <- gsub(" ", "_", moa.name)
  }
  return(moa.name)
}

load.CMap <- function(gene.symbol.type="19Q4"){
  # load gmt
  cat(file=stderr(), "About to get gmt", "\n")
  gmt <- GSA::GSA.read.gmt(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/MOA_gmt_file_n6_no_special_chars.gmt")

  # load PRISM drug AUC
  cat(file=stderr(), "About to get PRISM AUC data", "\n")
  PRISM.AUC <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv") #481 cell lines
  PRISM.AUC$X <- NULL

  # download RNAseq
  cat(file=stderr(), "About to get CCLE RNAseq data v19Q4", "\n")
  if(gene.symbol.type=="19Q4"){
    download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin",
                  destfile = "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
    load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
    download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin",
                  destfile = "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
    load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
  }else if(gene.symbol.type=="approved"){
    download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_1-200.Rbin",
                  destfile = "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_1-200.Rbin")
    load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_1-200.Rbin")
    download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_201-327.Rbin",
                  destfile = "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_201-327.Rbin")
    load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_201-327.Rbin")
  }
  RNA.df <- rbind(RNA.first200, RNA.rest)

  return(list(gmt=gmt, PRISM.AUC=PRISM.AUC, RNA.df=RNA.df))
}

# Define UI for application
ui <- fluidPage(
  # Application title
  titlePanel("Drug Mechanism Enrichment Analysis: Web Application"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      # bubble select: drug list or gene signature? default: drug list
      radioButtons(inputId = "type", label = "Would you like to input a drug rank list or gene signature?",
                   choices = c("Drug rank list from CMap L1000 query" = "L1000",
                               "Drug rank list from CMap PRISM query" = "PRISM",
                               "Drug rank list" = "drug",
                               "Gene signature" = "gene")),

      conditionalPanel(
        condition = "input.type == 'drug'",
        # if drug list, moa provided in third column? default: FALSE
        checkboxInput(inputId = "moa.included", label = "Optional: are mechanism-of-action (moa) annotations in the third column of your input?",
                      value = FALSE),
        # is score averaging needed? default: FALSE
        checkboxInput(inputId = "drug.avg", label = "Optional: are there multiple entries (rows) for each drug? If you select this option, the rank metric will be averaged for each drug.",
                      value = FALSE)
      ),

      conditionalPanel(
        condition = "input.type == 'gene'",
        # if gene signature, is score averaging needed? default: FALSE
        checkboxInput(inputId = "gene.avg", label = "Optional: are there multiple entries (rows) for each gene? If you select this option, the rank metric will be averaged for each gene.",
                      value = FALSE)
      ),

      # accept file input
      conditionalPanel(
        condition = "input.type == 'L1000' || input.type == 'PRISM'",
        fileInput(inputId = "gct", label = paste0("Upload a .gct file from your CMap output (",MB.limit," MB limit)"),
                  accept = ".gct")
      ),

      conditionalPanel(
        condition = "input.type == 'drug' || input.type == 'gene'",
        fileInput(inputId = "csv", label = paste0("Upload a .csv file with names in first column and ranks in second column (",MB.limit," MB limit)"),
                  accept = ".csv")
      ),

      # allow user to choose examples
      selectInput(inputId = "example", label = "Or choose an example below",
                  choices = c("No example selected" = "none",
                              "CMap L1000: JEKO1 treated with proteasome inhibitor bortezomib" = "Drug_rank_list/CMap/L1000/GSE14003_Proteasome_inhibitor_bortezomib_treated_JEKO1_cells_10H_vs_untreated",
                              "CMap L1000: CD34 treated with glucocorticoid receptor agonist dexamethasone" = "Drug_rank_list/CMap/L1000/GSE28896_Glucocorticoid_agonist_dexamethasone_treated_CD34_cells_24H_vs_untreated",
                              "CMap L1000: HUVEC treated with HMGCR inhibitor pitavastatin" = "Drug_rank_list/CMap/L1000/GSE32547_Pitavastatin_treated_HUVEC_cells_1_uM_at_4H_vs_DMSO_treated",
                              "CMap L1000: A2058 treated with PI3K/MTOR inhibitor BEZ235" = "Drug_rank_list/CMap/L1000/GSE33643_PI3K_MTOR_inhibitor_BEZ235_treated_A2058_cells_3_doses_at_24H_vs_DMSO_treated",
                              "CMap L1000: A375 treated with MEK inhibitor GSK212" = "Drug_rank_list/CMap/L1000/GSE35230_MEK_inhibitor_GSK212_treated_A375_clones_30nM_at_24H_vs_DMSO_treated",
                              "CMap L1000: Senescent HMEC" = "Drug_rank_list/CMap/L1000/Delfarah_et_al_HMEC_senescence",
                              "CMap PRISM: Sensitive to HMGCR inhibitor lovastatin" = "Drug_rank_list/CMap/PRISM/Cell_lines_sensitive_to_HMGCR_inhibitor_lovastatin",
                              "CMap PRISM: EGFR-activating mutation" = "Drug_rank_list/CMap/PRISM/Cell_lines_with_EGFR_activating_mutation_ELREA746DEL",
                              "CMap PRISM: High expression of PDGFRA" = "Drug_rank_list/CMap/PRISM/Cell_lines_with_high_gene_expression_of_PDGFRA",
                              "Gene signature: NSCLC sensitive to EGFR inhibitor gefitinib" = "Gene_signature/Coldren_et_al_NSCLC_sensitive_vs_resistant_to_gefitinib",
                              "Gene signature: BRCA sensitive to EGFR inhibitor erlotinib" = "Gene_signature/GSE12790_BRCA_sensitive_vs_resistant_to_erlotinib",
                              "Gene signature: NSCLC sensitive to EGFR inhibitor erlotinib" = "Gene_signature/GSE31625_NSCLC_sensitive_vs_resistant_to_erlotinib",
                              "Gene signature: SKCM sensitive to RAF inhibitor vemurafenib" = "Gene_signature/GSE66539_SKCM_sensitive_vs_resistant_to_vemurafenib",
                              "Gene signature: Senescent HMEC" = "Gene_signature/Delfarah_et_al_HMEC_senescence")),

      # get drug moa of interest (if any)
      textInput(inputId = "interest", label = "Optional: enter a moa of interest to view its mountain plot (case-sensitive; e.g., HMGCR inhibitor)"),

      # offer advanced settings
      checkboxInput(inputId = "advanced", label = "Optional: advanced settings", value = FALSE),
      conditionalPanel(
          condition = "input.advanced",
          sliderInput("n.min.per.set","Minimum drugs per set:",
                    min=1, max=40, value=6),
          sliderInput("p.cutoff","P-value cutoff:",
                      min=0, max=1, value=0.05),
          sliderInput("FDR.cutoff","False discovery rate cutoff:",
                        min=0, max=1, value=0.25),
          radioButtons("conv.syn","Convert drug synonyms if no moa annotations are provided with input drug rank list?",
                       choices = list("Yes" = TRUE, "No" = FALSE), selected = TRUE),
          radioButtons("gene.type","If inputting a gene signature, are you using currently approved gene symbols or those from the CCLE 19Q4 release?",
                       choices = list("Current HGNC-approved gene symbols" = "approved", "CCLE 19Q4 gene symbols" = "19Q4"), selected = "19Q4")
      ),

      # include "Run" button
      actionButton(inputId = "run", label = "Run"),

      # including loading message
      conditionalPanel(
        condition = "input.run && !output.msg",
        textOutput(outputId = "fyi")
      ),

      # indicate when run is completed
      textOutput("msg")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      # display plots in tabs
      tabsetPanel(
        tabPanel("Volcano Plot", plotOutput("volcPlot"),
                 conditionalPanel(
                   condition = "output.msg=='Run completed'",
                   downloadButton(outputId = "volcDownload", label = "Download volcano plot")
                 )),
        tabPanel("Mountain Plot", plotOutput("mtnPlot"),
                 conditionalPanel(
                   condition = "output.msg=='Run completed'",
                   downloadButton(outputId = "mtnDownload", label = "Download mountain plot")
                 )),
      ),

      # output result files (.zip)
      conditionalPanel(
        condition = "output.msg=='Run completed'",
        downloadButton(outputId = "results", label = "Download results")
      ),

      uiOutput("info"),
      textOutput("private"),
      textOutput("refresh")
    )
  )
)

# Define backend
server <- function(input, output) {
  url <- a("https://belindabgarana.github.io/DMEA", href = "https://belindabgarana.github.io/DMEA")
  output$info <- renderUI({tagList("For more information or to contact us, please visit: ", url)})
  output$private <- renderText({"No user data is stored on our secure server, so your data will remain private."})
  output$refresh <- renderText({"Please refresh your web browser after each analysis and format your inputs to match the examples on the 'How to Use' page at the url above to avoid errors. You will also need to refresh this webpage after 5 minutes of inactivity."})
  output$fyi <- renderText({"Running analysis... Please allow 1 to 5 minutes of run time"})
  options(shiny.maxRequestSize = MB.limit*1024^2)
  observeEvent(input$run, {
    selected.type <- input$type
    use.moa <- input$moa.included
    avg.drug <- input$drug.avg
    avg.gene <- input$gene.avg
    selected.ex <- input$example
    selected.moa <- input$interest
    min.per.set <- input$n.min.per.set
    p <- input$p.cutoff
    FDR <- input$FDR.cutoff
    convert.syn <- input$conv.syn
    gene.symbol.type <- input$gene.type
    if(selected.ex != "none"){
    ##### Using provided example #####
      # get last folder name for filename when outputting results
      selection.info <- strsplit(selected.ex, "/")[[1]]
      selected.case <- selection.info[length(selection.info)]

      if(sjmisc::str_contains(selected.ex, "Drug_rank_list/CMap")){
        #### Using drug rank list example ####
        # get inputs
        cat(file=stderr(), "About to read example input", "\n")
        rank.data <- read.csv(file=paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Examples/",selected.ex,"/DMEA_input.csv"))

        cat(file=stderr(), "About to read gmt", "\n")
        gmt <- readRDS(file=url(paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Examples/",selected.ex,"/DMEA_gmt.rds")))

        # run DMEA
        cat(file=stderr(), "About to run enrichment analysis", "\n")
        DMEA.output <- drugSEA(rank.data, gmt, drug="pert_iname", rank.metric=colnames(rank.data)[2],
                               FDR=FDR, min.per.set=min.per.set)

        # output results
        cat(file=stderr(), "About to output results", "\n")
        output$results <- downloadHandler(
          filename = function(){paste0("DMEA_results_",selected.case,"_",Sys.Date(),".zip")},
          content = function(file){
            mtn.plot.names <- c()
            for(i in 1:length(DMEA.output$mtn.plots)){
              temp.name <- paste0("DMEA_mountain_plot_",as.filename(names(DMEA.output$mtn.plots)[[i]]),".pdf")
              mtn.plot.names <- c(mtn.plot.names, temp.name)
              ggsave(temp.name, DMEA.output$mtn.plots[[i]], device="pdf")
            }
            all.files <- c("DMEA_input.csv",
                           "DMEA_gmt.Rbin",
                           "DMEA_results.csv",
                           "DMEA_removed_drug_sets.csv",
                           "DMEA_unannotated_drugs.csv",
                           "DMEA_volcano_plot.pdf",
                           mtn.plot.names)
            write.csv(rank.data, all.files[1], row.names = FALSE)
            saveRDS(gmt, all.files[2])
            write.csv(DMEA.output$result, all.files[3], row.names = FALSE)
            write.csv(DMEA.output$removed.sets, all.files[4], row.names = FALSE)
            write.csv(DMEA.output$unannotated.drugs, all.files[5], row.names = FALSE)
            ggsave(all.files[6], DMEA.output$volcano.plot, device="pdf")
            zip(zipfile=file, files=all.files)}
        )
      }else if(sjmisc::str_contains(selected.ex, "Gene_signature")){
        #### Using gene signature example ####
        # get inputs
        cat(file=stderr(), "About to read example input", "\n")
        rank.data <- read.csv(file=paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Examples/",selected.ex,"/Filtered_gene_signature_no_duplicates.csv"))
        inputs <- load.CMap()

        # run DMEA
        cat(file=stderr(), "About to run enrichment analysis", "\n")
        DMEA.output <- DMEA(inputs$PRISM.AUC, inputs$gmt, inputs$RNA.df, weights=rank.data, ylab="Drug AUC",
                            FDR=FDR, min.per.set=min.per.set)

        # output results
        cat(file=stderr(), "About to output results", "\n")
        output$results <- downloadHandler(
          filename = function(){paste0("DMEA_results_",selected.case,"_",Sys.Date(),".zip")},
          content = function(file){
            mtn.plot.names <- c()
            for(i in 1:length(DMEA.output$mtn.plots)){
              temp.name <- paste0("DMEA_mountain_plot_",as.filename(names(DMEA.output$mtn.plots)[[i]]),".pdf")
              mtn.plot.names <- c(mtn.plot.names, temp.name)
              ggsave(temp.name, DMEA.output$mtn.plots[[i]], device="pdf")
            }
            all.files <- c("DMEA_input.csv",
                           "DMEA_WGV_scores.csv",
                           "DMEA_unused_gene_weights.csv",
                           "DMEA_correlation_results.csv",
                           "DMEA_correlation_plots.pdf",
                           "DMEA_gmt.Rbin",
                           "DMEA_results.csv",
                           "DMEA_removed_drug_sets.csv",
                           "DMEA_unannotated_drugs.csv",
                           "DMEA_volcano_plot.pdf",
                           mtn.plot.names)
            write.csv(rank.data, all.files[1], row.names = FALSE)
            write.csv(DMEA.output$WV.scores, all.files[2], row.names = FALSE)
            write.csv(DMEA.output$unused.weights, all.files[3], row.names = FALSE)
            write.csv(DMEA.output$corr.result, all.files[4], row.names = FALSE)
            ggsave(all.files[5], DMEA.output$corr.scatter.plots, device="pdf")
            saveRDS(inputs$gmt, all.files[6])
            write.csv(DMEA.output$result, all.files[7], row.names = FALSE)
            write.csv(DMEA.output$removed.sets, all.files[8], row.names = FALSE)
            write.csv(DMEA.output$unannotated.drugs, all.files[9], row.names = FALSE)
            ggsave(all.files[10], DMEA.output$volcano.plot, device="pdf")
            zip(zipfile=file, files=all.files)}
        )
      }

      output$volcDownload <- downloadHandler(
        filename = function(){paste0("DMEA_volcano_plot_",selected.case,"_",Sys.Date(),".pdf")},
        content = function(file){ggsave(file, DMEA.output$volcano.plot, device="pdf")})
      output$volcPlot <- renderPlot({DMEA.output$volcano.plot})

      # if no moa is selected, display mountain plot for top moa passing p & FDR thresholds
      if(selected.moa==""){
        results.table <- DMEA.output$result
        sig.results <- results.table[results.table$p_value < p & results.table$FDR_q_value < FDR,]
        if(nrow(sig.results)>=1){
          top.result <- sig.results %>% slice_max(abs(sig.results$NES), n=1, with_ties = FALSE)
          selected.moa <- top.result$Drug_set
        }
      }

      output$mtnDownload <- downloadHandler(
        filename = function(){paste0("DMEA_mountain_plot_",as.filename(selected.moa),"_",selected.case,"_",Sys.Date(),".pdf")},
        content = function(file){ggsave(file, DMEA.output$mtn.plots[[selected.moa]], device="pdf")})
      output$mtnPlot <- renderPlot({DMEA.output$mtn.plots[[selected.moa]]})

      output$msg <- renderText({"Run completed"})
    }else if(selected.type == "L1000"){
      ##### Using CMap L1000 drug rank list input #####
      # get input
      cat(file=stderr(), "About to read input", "\n")
      rank.info <- input$gct
      rank.data <- read.delim(file=rank.info$datapath, skip = 2)
      rank.name <- paste0("CMap_L1000_",rank.data[1,ncol(rank.data)])

      # filter drug list
      cat(file=stderr(), "About to filter input", "\n")
      rank.data <- rank.data[rank.data$qc_pass==1,]
      rank.data <- rank.data[rank.data$pert_iname!="na",] # removes empty first row in CMap query output
      rank.data <- rank.data[rank.data$moa!="-666",] # remove drugs without known moa (moa="-666")
      rank.data[,ncol(rank.data)] <- as.numeric(rank.data[,ncol(rank.data)]) # assuming rank metric is last column
      rank.data <- rank.data[rank.data[,ncol(rank.data)]!=0,]
      rank.metric <- colnames(rank.data)[ncol(rank.data)] # store name of rank metric to restore later
      colnames(rank.data)[ncol(rank.data)] <- "rank_metric"

      # average rank metric
      cat(file=stderr(), "About to average rank metric", "\n")
      rank.data <- plyr::ddply(rank.data, .(pert_iname, moa), summarize, rank_metric = mean(rank_metric, na.rm=TRUE)) #average across cell lines
      rank.data <- dplyr::distinct(na.omit(rank.data[!duplicated(rank.data[,c("pert_iname")]), c("pert_iname", "rank_metric", "moa")]))
      colnames(rank.data)[2] <- rank.metric # restore name of rank metric

      # run DMEA
      cat(file=stderr(), "About to run enrichment analysis", "\n")
      DMEA.output <- drugSEA(rank.data, drug="pert_iname", rank.metric=rank.metric,
                             FDR=FDR, min.per.set=min.per.set)

      # output results
      cat(file=stderr(), "About to output results", "\n")
      output$results <- downloadHandler(
        filename = function(){paste0("DMEA_results_",rank.name,"_",Sys.Date(),".zip")},
        content = function(file){
          mtn.plot.names <- c()
          for(i in 1:length(DMEA.output$mtn.plots)){
            temp.name <- paste0("DMEA_mountain_plot_",as.filename(names(DMEA.output$mtn.plots)[[i]]),".pdf")
            mtn.plot.names <- c(mtn.plot.names, temp.name)
            ggsave(temp.name, DMEA.output$mtn.plots[[i]], device="pdf")
          }
          all.files <- c("DMEA_input.csv",
                         "DMEA_gmt.Rbin",
                         "DMEA_results.csv",
                         "DMEA_removed_drug_sets.csv",
                         "DMEA_unannotated_drugs.csv",
                         "DMEA_volcano_plot.pdf",
                         mtn.plot.names)
          write.csv(rank.data, all.files[1], row.names = FALSE)
          saveRDS(DMEA.output$gmt, all.files[2])
          write.csv(DMEA.output$result, all.files[3], row.names = FALSE)
          write.csv(DMEA.output$removed.sets, all.files[4], row.names = FALSE)
          write.csv(DMEA.output$unannotated.drugs, all.files[5], row.names = FALSE)
          ggsave(all.files[6], DMEA.output$volcano.plot, device="pdf")
          zip(zipfile=file, files=all.files)}
      )

      output$volcDownload <- downloadHandler(
        filename = function(){paste0("DMEA_volcano_plot_",rank.name,"_",Sys.Date(),".pdf")},
        content = function(file){ggsave(file, DMEA.output$volcano.plot, device="pdf")})
      output$volcPlot <- renderPlot({DMEA.output$volcano.plot})

      # if no moa is selected, display mountain plot for top moa passing p & FDR thresholds
      if(selected.moa==""){
        results.table <- DMEA.output$result
        sig.results <- results.table[results.table$p_value < p & results.table$FDR_q_value < FDR,]
        if(nrow(sig.results)>=1){
          top.result <- sig.results %>% slice_max(abs(sig.results$NES), n=1, with_ties = FALSE)
          selected.moa <- top.result$Drug_set
        }
      }

      output$mtnDownload <- downloadHandler(
        filename = function(){paste0("DMEA_mountain_plot_",as.filename(selected.moa),"_",rank.name,"_",Sys.Date(),".pdf")},
        content = function(file){ggsave(file, DMEA.output$mtn.plots[[selected.moa]], device="pdf")})
      output$mtnPlot <- renderPlot({DMEA.output$mtn.plots[[selected.moa]]})

      output$msg <- renderText({"Run completed"})
    }else if(selected.type == "PRISM"){
      ##### Using CMap PRISM drug rank list input #####
      # get input
      cat(file=stderr(), "About to read input", "\n")
      rank.info <- input$gct
      rank.data <- read.delim(file=rank.info$datapath, skip = 2)
      rank.name <- paste0("CMap_PRISM_",rank.data[1,ncol(rank.data)])

      # filter drug list
      cat(file=stderr(), "About to filter input", "\n")
      rank.data <- rank.data[rank.data$pert_iname!="na",] # removes empty first row in CMap query output
      rank.data <- rank.data[rank.data$moa!="-666",] # remove drugs without known moa (moa="-666")
      rank.data[,ncol(rank.data)] <- as.numeric(rank.data[,ncol(rank.data)]) # assuming rank metric is last column
      rank.data <- rank.data[rank.data[,ncol(rank.data)]!=0,]
      rank.metric <- colnames(rank.data)[ncol(rank.data)]
      rank.data <- distinct(na.omit(rank.data[!duplicated(rank.data[,c("pert_iname")]),c("pert_iname", rank.metric, "moa")]))

      # run DMEA
      cat(file=stderr(), "About to run enrichment analysis", "\n")
      DMEA.output <- drugSEA(rank.data, drug="pert_iname", rank.metric=rank.metric,
                             FDR=FDR, min.per.set=min.per.set)

      # output results
      cat(file=stderr(), "About to output results", "\n")
      output$results <- downloadHandler(
        filename = function(){paste0("DMEA_results_",rank.name,"_",Sys.Date(),".zip")},
        content = function(file){
          mtn.plot.names <- c()
          for(i in 1:length(DMEA.output$mtn.plots)){
            temp.name <- paste0("DMEA_mountain_plot_",as.filename(names(DMEA.output$mtn.plots)[[i]]),".pdf")
            mtn.plot.names <- c(mtn.plot.names, temp.name)
            ggsave(temp.name, DMEA.output$mtn.plots[[i]], device="pdf")
          }
          all.files <- c("DMEA_input.csv",
                         "DMEA_gmt.Rbin",
                         "DMEA_results.csv",
                         "DMEA_removed_drug_sets.csv",
                         "DMEA_unannotated_drugs.csv",
                         "DMEA_volcano_plot.pdf",
                         mtn.plot.names)
          write.csv(rank.data, all.files[1], row.names = FALSE)
          saveRDS(DMEA.output$gmt, all.files[2])
          write.csv(DMEA.output$result, all.files[3], row.names = FALSE)
          write.csv(DMEA.output$removed.sets, all.files[4], row.names = FALSE)
          write.csv(DMEA.output$unannotated.drugs, all.files[5], row.names = FALSE)
          ggsave(all.files[6], DMEA.output$volcano.plot, device="pdf")
          zip(zipfile=file, files=all.files)}
      )

      output$volcDownload <- downloadHandler(
        filename = function(){paste0("DMEA_volcano_plot_",rank.name,"_",Sys.Date(),".pdf")},
        content = function(file){ggsave(file, DMEA.output$volcano.plot, device="pdf")})
      output$volcPlot <- renderPlot({DMEA.output$volcano.plot})

      # if no moa is selected, display mountain plot for top moa passing p & FDR thresholds
      if(selected.moa==""){
        results.table <- DMEA.output$result
        sig.results <- results.table[results.table$p_value < p & results.table$FDR_q_value < FDR,]
        if(nrow(sig.results)>=1){
          top.result <- sig.results %>% slice_max(abs(sig.results$NES), n=1, with_ties = FALSE)
          selected.moa <- top.result$Drug_set
        }
      }

      output$mtnDownload <- downloadHandler(
        filename = function(){paste0("DMEA_mountain_plot_",as.filename(selected.moa),"_",rank.name,"_",Sys.Date(),".pdf")},
        content = function(file){ggsave(file, DMEA.output$mtn.plots[[selected.moa]], device="pdf")})
      output$mtnPlot <- renderPlot({DMEA.output$mtn.plots[[selected.moa]]})

      output$msg <- renderText({"Run completed"})
    }else if(selected.type == "drug" & use.moa == TRUE){
      ##### Using drug rank list input with moa annotations #####
      # get input
      cat(file=stderr(), "About to get input", "\n")
      rank.info <- input$csv
      rank.data <- read.csv(file=rank.info$datapath)
      rank.name <- substr(rank.info$name, 1, nchar(rank.info$name)-4)

      # average rank metric if needed
      if(avg.drug == TRUE){
        cat(file=stderr(), "About to average rank metric", "\n")
        rank.metric <- colnames(rank.data)[2] # store name of rank metric to restore later
        colnames(rank.data)[1:3] <- c("Drug", "rank_metric", "moa")
        rank.data$rank_metric <- as.numeric(rank.data$rank_metric)
        rank.data <- rank.data[rank.data$rank_metric!=0,]
        rank.data <- plyr::ddply(rank.data, .(Drug, moa), summarize, rank_metric = mean(rank_metric, na.rm=TRUE)) #average score for each drug
        rank.data <- distinct(na.omit(rank.data[!duplicated(rank.data[,c("Drug")]),c("Drug", "rank_metric", "moa")]))
        colnames(rank.data)[2] <- rank.metric # restore name of rank metric
      }

      # run DMEA
      cat(file=stderr(), "About to run enrichment analysis", "\n")
      DMEA.output <- drugSEA(rank.data, drug=colnames(rank.data)[1],
                             rank.metric=colnames(rank.data)[2], set.type=colnames(rank.data)[3],
                             FDR=FDR, min.per.set=min.per.set)

      # output results
      cat(file=stderr(), "About to output results", "\n")
      output$results <- downloadHandler(
        filename = function(){paste0("DMEA_results_",rank.name,"_",Sys.Date(),".zip")},
        content = function(file){
          mtn.plot.names <- c()
          for(i in 1:length(DMEA.output$mtn.plots)){
            temp.name <- paste0("DMEA_mountain_plot_",as.filename(names(DMEA.output$mtn.plots)[[i]]),".pdf")
            mtn.plot.names <- c(mtn.plot.names, temp.name)
            ggsave(temp.name, DMEA.output$mtn.plots[[i]], device="pdf")
          }
          all.files <- c("DMEA_input.csv",
                         "DMEA_gmt.Rbin",
                         "DMEA_results.csv",
                         "DMEA_removed_drug_sets.csv",
                         "DMEA_unannotated_drugs.csv",
                         "DMEA_volcano_plot.pdf",
                         mtn.plot.names)
          write.csv(rank.data, all.files[1], row.names = FALSE)
          saveRDS(DMEA.output$gmt, all.files[2])
          write.csv(DMEA.output$result, all.files[3], row.names = FALSE)
          write.csv(DMEA.output$removed.sets, all.files[4], row.names = FALSE)
          write.csv(DMEA.output$unannotated.drugs, all.files[5], row.names = FALSE)
          ggsave(all.files[6], DMEA.output$volcano.plot, device="pdf")
          zip(zipfile=file, files=all.files)}
      )

      output$volcDownload <- downloadHandler(
        filename = function(){paste0("DMEA_volcano_plot_",rank.name,"_",Sys.Date(),".pdf")},
        content = function(file){ggsave(file, DMEA.output$volcano.plot, device="pdf")})
      output$volcPlot <- renderPlot({DMEA.output$volcano.plot})

      # if no moa is selected, display mountain plot for top moa passing p & FDR thresholds
      if(selected.moa==""){
        results.table <- DMEA.output$result
        sig.results <- results.table[results.table$p_value < p & results.table$FDR_q_value < FDR,]
        if(nrow(sig.results)>=1){
          top.result <- sig.results %>% slice_max(abs(sig.results$NES), n=1, with_ties = FALSE)
          selected.moa <- top.result$Drug_set
        }
      }

      output$mtnDownload <- downloadHandler(
        filename = function(){paste0("DMEA_mountain_plot_",as.filename(selected.moa),"_",rank.name,"_",Sys.Date(),".pdf")},
        content = function(file){ggsave(file, DMEA.output$mtn.plots[[selected.moa]], device="pdf")})
      output$mtnPlot <- renderPlot({DMEA.output$mtn.plots[[selected.moa]]})

      output$msg <- renderText({"Run completed"})
    }else if(selected.type == "drug" & use.moa == FALSE){
      ##### Using drug rank list input without moa annotations #####
      # get input and remove special characters from drug names
      cat(file=stderr(), "About to get input", "\n")
      rank.info <- input$csv
      rank.data <- read.csv(file=rank.info$datapath)
      rank.name <- substr(rank.info$name, 1, nchar(rank.info$name)-4)
      special.chars <- c(" ", "-", "'", "/", "_", "[+]", "[(]", "[)]")
      for(j in 1:length(special.chars)){
        rank.data[,1] <- gsub(special.chars[j],".",rank.data[,1])
      }

      # average rank metric if needed
      if(avg.drug == TRUE){
        cat(file=stderr(), "About to average rank metric", "\n")
        rank.metric <- colnames(rank.data)[2] # store name of rank metric to restore later
        colnames(rank.data)[1:2] <- c("Drug", "rank_metric")
        rank.data$rank_metric <- as.numeric(rank.data$rank_metric)
        rank.data <- rank.data[rank.data$rank_metric!=0,]
        rank.data <- plyr::ddply(rank.data, .(Drug), summarize, rank_metric = mean(rank_metric, na.rm=TRUE)) #average score for each drug
        rank.data <- distinct(na.omit(rank.data[!duplicated(rank.data[,c("Drug")]),]))
        colnames(rank.data)[2] <- rank.metric # restore name of rank metric
      }

      # load gmt
      cat(file=stderr(), "About to get gmt", "\n")
      if(!require(GSA)){install.packages(GSA, repos = "http://cran.us.r-project.org")}
      library(GSA);
      gmt <- GSA.read.gmt(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/MOA_gmt_file_n6_no_special_chars.gmt")

      # run DMEA
      cat(file=stderr(), "About to run enrichment analysis", "\n")
      DMEA.output <- drugSEA(rank.data, gmt, drug=colnames(rank.data)[1], rank.metric=colnames(rank.data)[2],
                             FDR=FDR, min.per.set=min.per.set, convert.synonyms=convert.syn)

      # output results
      cat(file=stderr(), "About to output results", "\n")
      output$results <- downloadHandler(
        filename = function(){paste0("DMEA_results_",selected.ex,"_",Sys.Date(),".zip")},
        content = function(file){
          mtn.plot.names <- c()
          for(i in 1:length(DMEA.output$mtn.plots)){
            temp.name <- paste0("DMEA_mountain_plot_",as.filename(names(DMEA.output$mtn.plots)[[i]]),".pdf")
            mtn.plot.names <- c(mtn.plot.names, temp.name)
            ggsave(temp.name, DMEA.output$mtn.plots[[i]], device="pdf")
          }
          all.files <- c("DMEA_input.csv",
                         "DMEA_gmt.Rbin",
                         "DMEA_results.csv",
                         "DMEA_removed_drug_sets.csv",
                         "DMEA_replaced_drug_synonyms.csv",
                         "DMEA_unannotated_drugs.csv",
                         "DMEA_volcano_plot.pdf",
                         mtn.plot.names)
          write.csv(rank.data, all.files[1], row.names = FALSE)
          saveRDS(DMEA.output$gmt, all.files[2])
          write.csv(DMEA.output$result, all.files[3], row.names = FALSE)
          write.csv(DMEA.output$removed.sets, all.files[4], row.names = FALSE)
          write.csv(DMEA.output$replaced.drugs, all.files[5], row.names = FALSE)
          write.csv(DMEA.output$unannotated.drugs, all.files[6], row.names = FALSE)
          ggsave(all.files[7], DMEA.output$volcano.plot, device="pdf")
          zip(zipfile=file, files=all.files)}
      )

      output$volcDownload <- downloadHandler(
        filename = function(){paste0("DMEA_volcano_plot_",rank.name,"_",Sys.Date(),".pdf")},
        content = function(file){ggsave(file, DMEA.output$volcano.plot, device="pdf")})
      output$volcPlot <- renderPlot({DMEA.output$volcano.plot})

      # if no moa is selected, display mountain plot for top moa passing p & FDR thresholds
      if(selected.moa==""){
        results.table <- DMEA.output$result
        sig.results <- results.table[results.table$p_value < p & results.table$FDR_q_value < FDR,]
        if(nrow(sig.results)>=1){
          top.result <- sig.results %>% slice_max(abs(sig.results$NES), n=1, with_ties = FALSE)
          selected.moa <- top.result$Drug_set
        }
      }

      output$mtnDownload <- downloadHandler(
        filename = function(){paste0("DMEA_mountain_plot_",as.filename(selected.moa),"_",rank.name,"_",Sys.Date(),".pdf")},
        content = function(file){ggsave(file, DMEA.output$mtn.plots[[selected.moa]], device="pdf")})
      output$mtnPlot <- renderPlot({DMEA.output$mtn.plots[[selected.moa]]})

      output$msg <- renderText({"Run completed"})
    }else if(selected.type == "gene"){
      ##### Using gene signature input #####
      # get inputs
      cat(file=stderr(), "About to get input", "\n")
      rank.info <- input$csv
      rank.data <- read.csv(file=rank.info$datapath)
      rank.name <- substr(rank.info$name, 1, nchar(rank.info$name)-4)
      inputs <- load.CMap(gene.symbol.type)

      # average rank metric if needed
      if(avg.gene == TRUE){
        cat(file=stderr(), "About to average rank metric", "\n")
        rank.metric <- colnames(rank.data)[2] # store name of rank metric to restore later
        colnames(rank.data)[1:2] <- c("Gene", "rank_metric")
        rank.data$rank_metric <- as.numeric(rank.data$rank_metric)
        rank.data <- plyr::ddply(rank.data, .(Gene), summarize, rank_metric = mean(rank_metric, na.rm=TRUE)) #average score for each gene
        rank.data <- distinct(na.omit(rank.data[!duplicated(rank.data[,c("Gene")]),]))
        colnames(rank.data)[2] <- rank.metric # restore name of rank metric
      }

      # run DMEA
      cat(file=stderr(), "About to run enrichment analysis", "\n")
      if(nrow(rank.data)>0){
        DMEA.output <- DMEA(inputs$PRISM.AUC, inputs$gmt, inputs$RNA.df, weights=rank.data, ylab="Drug AUC",
                            FDR=FDR, min.per.set=min.per.set)

        # output results
        cat(file=stderr(), "About to output results", "\n")
        output$results <- downloadHandler(
          filename = function(){paste0("DMEA_results_",rank.name,"_",Sys.Date(),".zip")},
          content = function(file){
            mtn.plot.names <- c()
            for(i in 1:length(DMEA.output$mtn.plots)){
              temp.name <- paste0("DMEA_mountain_plot_",as.filename(names(DMEA.output$mtn.plots)[[i]]),".pdf")
              mtn.plot.names <- c(mtn.plot.names, temp.name)
              ggsave(temp.name, DMEA.output$mtn.plots[[i]], device="pdf")
            }
            all.files <- c("DMEA_input.csv",
                           "DMEA_WGV_scores.csv",
                           "DMEA_unused_gene_weights.csv",
                           "DMEA_correlation_results.csv",
                           "DMEA_correlation_plots.pdf",
                           "DMEA_gmt.Rbin",
                           "DMEA_results.csv",
                           "DMEA_removed_drug_sets.csv",
                           "DMEA_unannotated_drugs.csv",
                           "DMEA_volcano_plot.pdf",
                           mtn.plot.names)
            write.csv(rank.data, all.files[1], row.names = FALSE)
            write.csv(DMEA.output$WV.scores, all.files[2], row.names = FALSE)
            write.csv(DMEA.output$unused.weights, all.files[3], row.names = FALSE)
            write.csv(DMEA.output$corr.result, all.files[4], row.names = FALSE)
            ggsave(all.files[5], DMEA.output$corr.scatter.plots, device="pdf")
            saveRDS(inputs$gmt, all.files[6])
            write.csv(DMEA.output$result, all.files[7], row.names = FALSE)
            write.csv(DMEA.output$removed.sets, all.files[8], row.names = FALSE)
            write.csv(DMEA.output$unannotated.drugs, all.files[9], row.names = FALSE)
            ggsave(all.files[10], DMEA.output$volcano.plot, device="pdf")
            zip(zipfile=file, files=all.files)}
        )

        output$volcDownload <- downloadHandler(
          filename = function(){paste0("DMEA_volcano_plot_",rank.name,"_",Sys.Date(),".pdf")},
          content = function(file){ggsave(file, DMEA.output$volcano.plot, device="pdf")})
        output$volcPlot <- renderPlot({DMEA.output$volcano.plot})

        # if no moa is selected, display mountain plot for top moa passing p & FDR thresholds
        if(selected.moa==""){
          results.table <- DMEA.output$result
          sig.results <- results.table[results.table$p_value < p & results.table$FDR_q_value < FDR,]
          if(nrow(sig.results)>=1){
            top.result <- sig.results %>% slice_max(abs(sig.results$NES), n=1, with_ties = FALSE)
            selected.moa <- top.result$Drug_set
          }
        }

        output$mtnDownload <- downloadHandler(
          filename = function(){paste0("DMEA_mountain_plot_",as.filename(selected.moa),"_",rank.name,"_",Sys.Date(),".pdf")},
          content = function(file){ggsave(file, DMEA.output$mtn.plots[[selected.moa]], device="pdf")})
        output$mtnPlot <- renderPlot({DMEA.output$mtn.plots[[selected.moa]]})

        output$msg <- renderText({"Run completed"})
      }else{
        output$msg <- renderText({"No input genes were available in CCLE RNAseq v19Q4"})}
    }
  }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
