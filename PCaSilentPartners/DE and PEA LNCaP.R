#' ---
#' title: "Differential expression analysis and pathway enrichment analysis of LNCaP Co-Culture vs. Monoculture"
#' date: "02.08.2024"
#' ---

#'## Loading packages

library(BiocManager)
library(ggplot2)
library(tidyr)
library(dplyr)
library(openxlsx)

#'## Paths

# create shortcut paths
input <- "input_data/"
output_modified_data <- "output_data/LNCaP/0_modified_data/"
output_imputation <- "output_data/LNCaP/1_Imputation/"
output_DE_results <- "output_data/LNCaP/2_DE_results/"
output_Volcano_Plot <- "output_data/LNCaP/3_Volcano_Plots/"
output_pathfindR <- "output_data/LNCaP/4_pathfindR_results/"
output_GSEA_Reactome <- "output_data/LNCaP/5_GSEA_Reactome_results/"

#'## Import data

#' The Data contains values that are not a number (NaN) and blanks.
#' They are now assigned "NA" values.

raw_data <- read.delim(
  file = paste0(input, "input_data.txt"),
  comment.char = "#",
  na.strings = c("", "NaN")
)

data <- raw_data |>
  dplyr::rename(
    "1_CTRL_LNCaP_01" = "HOEMSG06_01_01",
    "1_CTRL_LNCaP_02" = "HOEMSG06_01_02",
    "2_CTRL_LNCaP" = "HOEMSG06_05_01",
    "3_CTRL_LNCaP" = "HOEMSG06_09_01",
    "1_Co_LNCaP" = "HOEMSG06_02_01",
    "2_Co_LNCaP" = "HOEMSG06_06_01",
    "3_Co_LNCaP" = "HOEMSG06_10_01"
  ) |>
  dplyr::select(
    "1_CTRL_LNCaP_01",
    "1_CTRL_LNCaP_02",
    "2_CTRL_LNCaP", 
    "3_CTRL_LNCaP", 
    "1_Co_LNCaP", 
    "2_Co_LNCaP", 
    "3_Co_LNCaP", 
    Genes)

#'## Modify data

#' Since single cells contain more than one gene symbol (isoforms,...) 
#' separated by semicolons. We use separate_rows from the tidyr package 
#' to solve this issue. Then we filter for NA vales within Gene symbols and
#' check for distinct Genes within the Genes column.

# separate rows and filter NA values in Gene symbols
data_filtered <- data |>
  tidyr::separate_rows(Genes, 
                       sep = ";") |>
  dplyr::filter(!is.na(Genes))

# filter for unique Gene names, keeps the first occurrence
data_filtered <- data_filtered |>
  dplyr::distinct(Genes, 
                  .keep_all = TRUE)

# save the modified data in an .xlsx file.
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb = wb, 
                       sheetName = "modified data LNCaP Co vs. Mono")

openxlsx::writeData(wb = wb,
                    sheet = "modified data LNCaP Co vs. Mono", 
                    x = data_filtered)

openxlsx::saveWorkbook(wb = wb, 
                       file = paste0(output_modified_data, "modified data LNCaP Co-Culture vs. Monoculture.xlsx"),
                       overwrite = TRUE)

#'## Imputation

#' Now we impute missing values using missForest.
#' Before the imputation takes place, we first filter for NA values:
#' max 1 NA in each condition.

# Function to count NA values in each group and the whole row
count_nas <- function(row) {
  row <- as.numeric(row)
  CTRL_NAs <- sum(is.na(row[1:4]))
  Co_NAs <- sum(is.na(row[5:7]))
  Total_NAs <- sum(is.na(row))
  
  return(list(CTRL_NAs = CTRL_NAs, 
              Co_NAs = Co_NAs, 
              Total_NAs = Total_NAs))
}

# Apply the function to each row and filter based on the conditions
data_for_imputation <- data_filtered |>
  dplyr::rowwise() |>
  dplyr::mutate(NA_counts = list(count_nas(c_across(1:7)))) |>
  dplyr::filter(
    NA_counts$CTRL_NAs <= 1 &
      NA_counts$Co_NAs <= 1 &
      NA_counts$Total_NAs <= 2
  ) |>
  dplyr::select(-NA_counts)

# select columns containing only numbers
data_numbers <- as.matrix(x = data_for_imputation[, 1:7])

rownames(data_numbers) <- data_for_imputation$Genes

# Imputation

# Open a file connection for writing messages
message_file <- file(paste0(output_imputation, "message imputation LNCaP Co-Culture vs. Monoculture.txt"), 
                     open = "wt")

sink(file = message_file, 
     append = TRUE, 
     type = "output")

sink(file = message_file, 
     append = TRUE, 
     type = "message")

library(missForest)

# starting time
current_timezone <- Sys.timezone()

starting_time <- format(x = Sys.time(), 
                        tz = current_timezone, 
                        format = "%Y-%m-%d %H:%M:%S")

print(starting_time)

# print begin message
cat("starting Imputation \n")

# Imputation
imputed_missForest <- missForest::missForest(xmis = data_numbers,
                                             verbose = TRUE)$ximp

# gene information to imputation

imputation_data <- as.data.frame(imputed_missForest)

imputation_data$Genes <- rownames(imputation_data)

# saving results in an .xlsx file
wb_Imputation <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb = wb_Imputation, 
                       sheetName = "Imputation LNCaP Co vs. Mono")

openxlsx::writeData(wb = wb_Imputation, 
                    sheet = "Imputation LNCaP Co vs. Mono", 
                    x = imputation_data)

openxlsx::saveWorkbook(wb = wb_Imputation, 
                       file = paste0(output_imputation, "Imputation LNCaP Co-Culture vs. Monoculture.xlsx"),
                       overwrite = TRUE)

# saving results .Rdata to easily load the file later
save(imputation_data, 
     file = paste0(output_imputation, "Imputation LNCaP Co-Culture vs. Monoculture.Rdata"))

cat("finished Imputation \n")

# finishing time
end_time <- format(x = Sys.time(), 
                   tz = current_timezone, 
                   format = "%Y-%m-%d %H:%M:%S")

print(end_time)

# stop redirecting messages
sink(type = "output")
sink(type = "message")

# Close the file connection
close(message_file)  

#'## Differential Expression Analysis with limma

library(limma)

# create a design matrix
group_data <- factor(c("CTRL", "CTRL", "CTRL", "CTRL",
                     "Co", "Co", "Co"),
                         levels = c("CTRL", "Co"))

biological_replicate <- c(1, 1, 2, 3, 1, 2, 3)

design <- model.matrix(~0 + group_data)

colnames(design) <- levels(group_data)

# Handle technical replicates
dupcor <- limma::duplicateCorrelation(object = imputed_missForest, 
                                      design  = design, 
                                      block = biological_replicate)

# define weights
aw <- limma::arrayWeights(object = imputed_missForest, 
                          design = design)

# fit linear model
fit <- limma::lmFit(object = imputed_missForest, 
                    design = design, 
                    block = biological_replicate, 
                    correlation = dupcor$consensus.correlation,
                    weights = aw)

# design contrasts
contrast_matrix <- limma::makeContrasts("LNCaP Co-Culture vs. Monoculture" = Co - CTRL,
                                        levels = design)

fit <- limma::contrasts.fit(fit = fit, 
                            contrasts = contrast_matrix)

# eBayes
fit <- limma::eBayes(fit = fit)

# saving results in a list containing the contrasts
DE <- limma::topTable(fit = fit, 
                      number=Inf, 
                      coef = "LNCaP Co-Culture vs. Monoculture", 
                      adjust.method = "BH")

DE$Genes <- rownames(DE)

rownames(DE) <- NULL

# saving results in an excel file
wb_DE <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb = wb_DE, 
                       sheetName = "LNCaP Co vs. Mono")

openxlsx::writeData(wb = wb_DE, 
                    sheet = "LNCaP Co vs. Mono", 
                    x = DE)

openxlsx::saveWorkbook(wb = wb_DE, 
                       file = paste0(output_DE_results, "DE LNCaP Co-Culture vs. Monoculture.xlsx"),
                       overwrite = TRUE)

#'## Volcano Plots

# create column with gene type: up, down, not significant
DE <- DE |>
  dplyr::mutate(gene_type = dplyr::case_when(logFC > 0 & adj.P.Val <= 0.05 ~ "up",
                                             logFC < 0 & adj.P.Val <= 0.05 ~ "down",
                                             TRUE ~ "not significant"))

# set colors
cols <- c("up" = "#ffad73", 
          "down" = "#26b3ff", 
          "not significant" = "grey") 

# Volcano plot
volcanoplot <- ggplot2::ggplot(data = DE,
                               mapping = aes(x = logFC, 
                                             y = -log10(adj.P.Val))) + 
  ggplot2::geom_point(mapping = aes(colour = gene_type), 
                      shape = 16,
                      size = 1.5) +
  ggplot2::geom_hline(yintercept = -log10(0.05),
                      linetype = "dashed") +
  ggplot2::geom_vline(xintercept = log2(1),
                      linetype = "dashed") +
  ggplot2::scale_x_continuous(breaks = c(seq(from = -6, 
                                             to = 6, 
                                             by = 1)),
                              limits = c(-6, 6)) +
  ggplot2::scale_y_continuous(breaks = c(seq(from = 0, 
                                             to = 6, 
                                             by = 1)),
                              limits = c(0, 6)) +
  ggplot2::scale_colour_manual(values = cols) +
  ggplot2::labs(title = paste("Volcano Plot: LNCaP Co-Culture vs. Monoculture"),
                x = "log2(fold change)",
                y = "-log10(adjusted p-value)",
                colour = "Expression \nchange") +
  ggplot2::theme(panel.border = element_rect(colour = "black", 
                                             fill = NA, 
                                             linewidth = 0.5),
                 legend.title = element_text(size = 12, 
                                             hjust = 0.5),
                 legend.text = element_text(size = 10),
                 plot.title = element_text(size = 13,
                                           hjust = 0.5),
                 axis.title.y = element_text(size = 10),
                 axis.title.x = element_text(size = 10),
                 text = element_text(family = "sans"),
                 panel.background = element_blank())

ggplot2::ggsave(filename = paste0(output_Volcano_Plot, "VP LNCaP Co-Culture vs. Monoculture.tiff"), 
                plot = volcanoplot,
                width = 5.5,
                height = 6,
                dpi = 800)

#+ fig.dim = c(5.5, 6)
volcanoplot

#'## Pathway enrichment analysis using pathfindR

# Open a file connection for writing messages
message_file <- file(paste0(output_pathfindR, "message pathfindR Reactome LNCaP Co-Culture vs. Monoculture.txt"), 
                     open = "wt")

sink(file = message_file, 
     append = TRUE, 
     type = "output")

sink(file = message_file, 
     append = TRUE, 
     type = "message")

# starting time
current_timezone <- Sys.timezone()

starting_time <- format(x = Sys.time(), 
                        tz = current_timezone, 
                        format = "%Y-%m-%d %H:%M:%S")

print(starting_time)

# pathway enrichment analysis

library(pathfindR)
library(org.Hs.eg.db)

# create workbooks for saving data in excel
wb <- openxlsx::createWorkbook()

# print begin message
cat("begin pathfindR Reactome LNCaP Co-Culture vs. Monoculture \n")

# beginning  time for one gene set within a DE
time <- format(x = Sys.time(), 
               tz = current_timezone, 
               format = "%Y-%m-%d %H:%M:%S")

print(time)

# pathfindR
diffdataexp <- DE |>
  dplyr::select(Genes, 
                logFC, 
                adj.P.Val)

pathfindR <- pathfindR::run_pathfindR(
  input = diffdataexp,
  gene_sets = "Reactome",
  min_gset_size = 10,
  max_gset_size = Inf,
  pin_name_path = "Biogrid",
  p_val_threshold = 0.05, # the p value threshold when filtering the input data
  adj_method = "bonferroni",
  enrichment_threshold = 0.05,
  convert2alias = TRUE,
  plot_enrichment_chart = FALSE,
  iterations = 10,
  list_active_snw_genes = TRUE,
  silent_option = TRUE
  )

# add worksheets and write data to worksheets
openxlsx::addWorksheet(wb = wb, 
                       sheetName = "LNCaP Co vs. Mono Reactome")

openxlsx::writeData(wb = wb,
                    sheet = "LNCaP Co vs. Mono Reactome", 
                    x = pathfindR)

# saving the workbook
openxlsx::saveWorkbook(wb = wb, 
                       file = paste0(output_pathfindR, "/output pathfindR Reactome LNCaP Co-Culture vs. Monoculture.xlsx"),
                       overwrite = TRUE)

cat("finished pathfindR Reactome LNCaP Co-Culture vs. Monoculture \n")

# finishing time 
end_time <- format(x = Sys.time(), 
                   tz = current_timezone, 
                   format = "%Y-%m-%d %H:%M:%S")

print(end_time)

# print the starting time and finishing time as well as the difference to see how long the process took
time_difference <- as.POSIXct(x = end_time, tz = current_timezone) - as.POSIXct(x = starting_time, tz = current_timezone)

print(paste("started", starting_time, "and finished", end_time))

print(time_difference)

# stop redirecting messages
sink(type = "output")
sink(type = "message")

# Close the file connection
close(message_file)  

# enrichment chart
# was later added so that the image is shorter:
# pathfindR[10, "Term_Description"] <- "Respiratory electron transport, ATP synthesis \nby chemiosmotic coupling, and heat production by uncoupling proteins."

enrichment_chart <- pathfindR::enrichment_chart(result_df = pathfindR,
                                                top_terms = 20) +
  ggplot2::scale_x_continuous(breaks = c(seq(from = 0, 
                                             to = 25, 
                                             by = 5)),
                              limits = c(0, 25)) +
  ggplot2::scale_size_continuous(breaks = seq(from = 0, 
                                              to = 300, 
                                              by = 50),
                                 limits = c(0, 300)) +
  ggplot2::scale_color_gradient(low = "#00CCFF", 
                                high = "#0000CC") +
  ggplot2::labs(title = paste0("LNCaP Co-Culture vs. Monoculture \nEnrichment bubble chart pathfindR Reactome")) +
  ggplot2::theme(legend.text = element_text(size = 10),
                 legend.title = element_text(size = 12),
                 plot.title = element_text(size = 16),
                 text = element_text(family = "sans"),
                 axis.text.y = element_text(size = 15)) 

ggplot2::ggsave(filename = paste0(output_pathfindR, "enrichment_chart/enrichment chart pathfindR Reactome LNCaP Co-Culture vs. Monoculture.tiff"), 
                plot = enrichment_chart,
                width = 14,
                height = 8,
                dpi = 800)

#+ fig.dim = c(14, 8)
enrichment_chart

#'##  Gene Set Enrichment Analysis Reactome

# Open a file connection for writing messages
message_file <- file(paste0(output_GSEA_Reactome, "message GSEA Reactome LNCaP Co-Culture vs. Monoculture.txt"), 
                     open = "wt")

sink(file = message_file, 
     append = TRUE, 
     type = "output")

sink(file = message_file, 
     append = TRUE, 
     type = "message")

library(clusterProfiler)
library(ReactomePA)

# starting time
current_timezone <- Sys.timezone()

starting_time <- format(x = Sys.time(), 
                        tz = current_timezone, 
                        format = "%Y-%m-%d %H:%M:%S")

print(starting_time)

# create workbooks for saving data in excel
wb <- openxlsx::createWorkbook(title = paste0("GSEA Reactome LNCaP Co-Culture vs. Monoculture"))

# print begin message
cat("begin GSEA Reactome LNCaP Co-Culture vs. Monoculture \n")

# prepare input
res <- DE

# transform Symbols into Entrez-IDs, which is required for gsePathway
ids <- clusterProfiler::bitr(res$Genes, 
                             fromType ="SYMBOL", 
                             toType = "ENTREZID", 
                             OrgDb = "org.Hs.eg.db", 
                             drop = FALSE)

# exclude duplicate
ids <- ids[ids$ENTREZID != "100187828",]

# attach Entrez-IDs to res
res$Entrez <- ids[,"ENTREZID"]

# remove NA, which were created by bitr
res <- res[complete.cases(res), ]

# ranking
rankings <- sign(res$logFC)*(-log10(res$adj.P.Val))

# genes as names
names(rankings) <- res$Entrez

# sort in decreasing order
rankings <- sort(rankings, 
                 decreasing = TRUE)

# Gene Set Enrichment Analysis
GSEA_Reactome <- ReactomePA::gsePathway(geneList = rankings,
                                        organism = "human",
                                        exponent = 1,
                                        minGSSize = 0,
                                        maxGSSize = Inf,
                                        eps = 1e-300,
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        verbose = TRUE,
                                        nPermSimple = 10000,
                                        by = "fgsea")

print(warnings())

# save as R.data to easily load the file later
save(GSEA_Reactome, 
     file = paste0(output_GSEA_Reactome, "GSEA Reactome LNCaP Co-Culture vs. Monoculture.Rdata"))

# add worksheets
openxlsx::addWorksheet(wb = wb, 
                       sheetName = "LNCaP Co vs. Mono")

openxlsx::writeData(wb = wb,
                    sheet = "LNCaP Co vs. Mono", 
                    x = GSEA_Reactome)

# saving the workbook
openxlsx::saveWorkbook(wb = wb, 
                       file = paste0(output_GSEA_Reactome, "GSEA Reactome LNCaP Co-Culture vs. Monoculture.xlsx"),
                       overwrite = TRUE)

cat("finished GSEA Reactome LNCaP Co-Culture vs. Monoculture \n")

# finishing time
end_time <- format(x = Sys.time(), 
                   tz = current_timezone, 
                   format = "%Y-%m-%d %H:%M:%S")
print(end_time)

# print the starting time and finishing time as well as the difference to see how long the process took
time_difference <- as.POSIXct(x = end_time, tz = current_timezone) - as.POSIXct(x = starting_time, tz = current_timezone)

print(paste("started", starting_time, "and finished", end_time))

print(time_difference)

# stop redirecting messages => see in directory messages.txt if errors occured and when
sink(type = "output")
sink(type = "message")

# Close the file connection
close(message_file)  

#'## Visualization of GSEA
#'### GSEA Plot Reactome

library(enrichplot)

# GSEA plot of selected pathways 
ids <- c("R-HSA-71291") # Metabolism of amino acids and derivatives

# GSEA plot of selected pathways saved in ids
gsea_plot_Reactome <- enrichplot::gseaplot2(x = GSEA_Reactome,
                                            geneSetID = ids,
                                            title = paste("GSEA Plot Reactome LNCaP Co-Culture vs. Monoculture"),
                                            pvalue_table = TRUE,
                                            color = "darkolivegreen3")

# the upper plot is saved within the gglist on position 1, 
# set y scale for graph comparison
gsea_plot_Reactome[[1]] <- gsea_plot_Reactome[[1]] +
  ggplot2::scale_y_continuous(breaks = c(seq(from = -0.3,
                                             to = 0.5, 
                                             by = 0.2)),
                              limits = c(-0.3, 0.5)) +
  ggplot2::theme(text = element_text(family = "sans"),
                 plot.title = element_text(size = 13,
                                           hjust = 0.5),
                 axis.title = element_text(size = 10))

# set alpha = 0 to not visualize colour
gsea_plot_Reactome[[2]]$layers[[2]]$aes_params$alpha <- 0 

# the lower plot is saved within gglist on position 3, 
# set x scale for graph comparison
gsea_plot_Reactome[[3]] <- gsea_plot_Reactome[[3]] +
  ggplot2::scale_y_continuous(breaks = c(seq(from = -15, 
                                             to = 12, 
                                             by = 5)),
                              limits = c(-15, 12))  +
  ggplot2::theme(text = element_text(family = "sans"),
                 axis.title = element_text(size = 10))

ggplot2::ggsave(filename = paste0(output_GSEA_Reactome, "GSEA Plot Reactome LNCaP Co-Culture vs. Monoculture.tiff"), 
                plot = gsea_plot_Reactome,
                width = 10,
                height = 7,
                dpi = 800)

#+ fig.dim = c(10, 7)
gsea_plot_Reactome

#'### Barplot NES Scores of selected pathways
# 13 / 15 pathways of the Metabolism tree of Reactome
Description_selected_pathways <- c("Metabolism of carbohydrates",
                                   "Inositol phosphate metabolism",
                                   "Metabolism of lipids",
                                   "Integration of energy metabolism",
                                   "Metabolism of nitric oxide: NOS3 activation and regulation",
                                   "Aerobic respiration and respiratory electron transport",
                                   "Metabolism of nucleotides",
                                   "Metabolism of vitamins and cofactors",
                                   "Metabolism of amino acids and derivatives",
                                   "Metabolism of porphyrins",
                                   "Biological oxidations",
                                   "Mitochondrial iron-sulfur cluster biogenesis",
                                   "Cytosolic iron-sulfur cluster assembly")

gsea_nes_score_Reactome <- GSEA_Reactome@result |>
  dplyr::filter(Description %in% Description_selected_pathways)

gsea_plot_nes <- ggplot2::ggplot(data = gsea_nes_score_Reactome,
                                 mapping = aes(x = NES, 
                                               y = Description, 
                                               fill = p.adjust < 0.05)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = c("TRUE" = "#ffad73", 
                                        "FALSE" = "#26b3ff")) +
  ggplot2::labs(x = "Normalized Enrichment Score",
                y = NULL,
                title = "LNCaP Co-Culture vs. Monoculture \nNES values Reactome") +
  ggplot2::geom_text(mapping = aes(label = sprintf("%.2f", NES)), 
                     position = position_stack(vjust = 0.5), 
                     color = "black") +
  ggplot2::scale_x_continuous(breaks = c(seq(from = -3, 
                                             to = 2, 
                                             by = 0.5)),
                              limits = c(-3, 2)) +
  ggplot2::theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
                 legend.title = element_text(size = 12, 
                                             hjust = 0.5),
                 legend.text = element_text(size = 10),
                 plot.title = element_text(size = 13,
                                           hjust = 0.5),
                 axis.text.y = element_text(size = 10),
                 text = element_text(family = "sans"),
                 panel.background = element_blank())

ggplot2::ggsave(filename = paste0(output_GSEA_Reactome, "NES Plot Reactome LNCaP Co-Culture vs. Monoculture.tiff"), 
                plot = gsea_plot_nes,
                width = 10,
                height = 7,
                dpi = 800)

#+ fig.dim = c(10, 7)
gsea_plot_nes

#'## session and packages
#' Here is important information about packages and session

# Open a file connection for writing messages
message_file <- file("packages LNCaP Co-Culture vs. Monoculture.txt", 
                     open = "wt")

sink(file = message_file, 
     append = TRUE, 
     type = "output")

sink(file = message_file, 
     append = TRUE, 
     type = "message")

# starting time
current_timezone <- Sys.timezone()
starting_time <- format(x = Sys.time(), 
                        tz = current_timezone, 
                        format = "%Y-%m-%d %H:%M:%S")
print(starting_time)

sessionInfo()

library(purrr)
c("BiocManager", 
  "ggplot2", 
  "tidyr",
  "dplyr",
  "openxlsx",
  "missForest", 
  "limma", 
  "pathfindR", 
  "org.Hs.eg.db", 
  "clusterProfiler",
  "ReactomePA",
  "enrichplot",
  "DOSE",
  "purrr") |>
  purrr::map(citation) |>
  print(style = "text")

java_version <- system("java -version", 
                       intern = TRUE)

print(java_version, 
      sep = "\n")

# stop redirecting messages
sink(type = "output")
sink(type = "message")

# Close the file connection
close(message_file)  
