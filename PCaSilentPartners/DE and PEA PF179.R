#' ---
#' title: "Differential expression analysis and pathway enrichment analysis of CAF Co-Culture vs. Monoculture"
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
output_modified_data <- "output_data/CAF/0_modified_data/"
output_imputation <- "output_data/CAF/1_Imputation/"
output_DE_results <- "output_data/CAF/2_DE_results/"
output_Volcano_Plot <- "output_data/CAF/3_Volcano_Plots/"
output_pathfindR <- "output_data/CAF/4_pathfindR_results/"

#'## Import data

#' The Data contains values that are not a number (NaN) and blanks.
#' They are now assigned "NA" values.

raw_data <- read.delim(
  file = paste0(input, "input_data.txt"),
  comment.char = "#",
  na.strings = c("", "NaN")
)

#' We need to rename the data and select relevant column for further analysis.
#' Here, we analyze CAF cells with Co-culutre of CAFs

data <- raw_data |>
  dplyr::rename(
    "1_Co_CAF" = "HOEMSG06_03_01",
    "2_Co_CAF" = "HOEMSG06_07_01",
    "3_Co_CAF" = "HOEMSG06_11_01",
    "1_CTRL_CAF" = "HOEMSG06_04_01",
    "2_CTRL_CAF" = "HOEMSG06_08_01",
    "3_CTRL_CAF" = "HOEMSG06_12_01",
  ) |>
  dplyr::select(
    "1_CTRL_CAF",
    "2_CTRL_CAF",
    "3_CTRL_CAF",
    "1_Co_CAF", 
    "2_Co_CAF", 
    "3_Co_CAF",
    Genes)

#'## Modify data

#' Since single cells contain more than one gene symbol (isoforms,...) 
#' separated by semicolons, we use separate_rows from the tidyr package 
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
                       sheetName = "modified data CAF Co vs. Mono")

openxlsx::writeData(wb = wb,
                    sheet = "modified data CAF Co vs. Mono", 
                    x = data_filtered)

openxlsx::saveWorkbook(wb = wb, 
                       file = paste0(output_modified_data, "modified data CAF Co-Culture vs. Monoculture.xlsx"),
                       overwrite = TRUE)

#'## Imputation

#' Now we impute missing values using missForest.
#' Before the imputation takes place, we first filter for NA values:
#' max 1 NA in each condition.

# Function to count NAs in each group and the whole row
count_nas <- function(row) {
  row <- as.numeric(row)
  CTRL_NAs <- sum(is.na(row[1:3]))
  Co_NAs <- sum(is.na(row[4:6]))
  Total_NAs <- sum(is.na(row))
  
  return(list(CTRL_NAs = CTRL_NAs, 
              Co_NAs = Co_NAs, 
              Total_NAs = Total_NAs))
}

# Apply the function to each row and filter based on the conditions
data_for_imputation <- data_filtered |>
  dplyr::rowwise() |>
  dplyr::mutate(NA_counts = list(count_nas(c_across(1:6)))) |>
  dplyr::filter(
    NA_counts$CTRL_NAs <= 1 &
      NA_counts$Co_NAs <= 1 &
      NA_counts$Total_NAs <= 2
  ) |>
  dplyr::select(-NA_counts)

# select columns containing only numbers
data_numbers <- as.matrix(x = data_for_imputation[, 1:6])

rownames(data_numbers) <- data_for_imputation$Genes

# Imputation

# Open a file connection for writing messages
message_file <- file(paste0(output_imputation, "message imputation CAF Co-Culture vs. Monoculture.txt"), 
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
                       sheetName = "Imputation CAF Co vs. Mono")

openxlsx::writeData(wb = wb_Imputation, 
                    sheet = "Imputation CAF Co vs. Mono", 
                    x = imputation_data)

openxlsx::saveWorkbook(wb = wb_Imputation, 
                       file = paste0(output_imputation, "Imputation CAF Co-Culture vs. Monoculture.xlsx"),
                       overwrite = TRUE)

# saving results .Rdata to easily load the file later
save(imputation_data, 
     file = paste0(output_imputation, "Imputation CAF Co-Culture vs. Monoculture.Rdata"))

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
group_data <- factor(rep(c("CTRL", "Co"), each = 3,
                          levels = c("CTRL", "Co")))

design <- model.matrix(~0 + group_data)

colnames(design) <- levels(group_data)

# define weights
aw <- limma::arrayWeights(object = imputed_missForest, 
                          design = design)

# fit linear model
fit <- limma::lmFit(imputed_missForest, 
                    design = design, 
                    weights = aw)

# design contrasts
contrast_matrix <- limma::makeContrasts("CAF Co-Culture vs. Monoculture" = Co - CTRL,
                                        levels = design)

fit <- limma::contrasts.fit(fit = fit, 
                            contrasts = contrast_matrix)

# eBayes
fit <- limma::eBayes(fit = fit)

# saving results in a list containing the contrasts
DE <- limma::topTable(fit = fit, 
                      number=Inf, 
                      coef = "CAF Co-Culture vs. Monoculture", 
                      adjust.method = "BH")

DE$Genes <- rownames(DE)

rownames(DE) <- NULL

# saving results in an excel file
wb_DE <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb = wb_DE, 
                       sheetName = "CAF Co vs. Mono")

openxlsx::writeData(wb = wb_DE, 
                    sheet = "CAF Co vs. Mono", 
                    x = DE)

openxlsx::saveWorkbook(wb = wb_DE, 
                       file = paste0(output_DE_results, "DE CAF Co-Culture vs. Monoculture.xlsx"),
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
  ggplot2::labs(title = paste("Volcano Plot: CAF Co-Culture vs. Monoculture"),
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

ggplot2::ggsave(filename = paste0(output_Volcano_Plot, "VP CAF Co-Culture vs. Monoculture.tiff"), 
                plot = volcanoplot,
                width = 5.5,
                height = 6,
                dpi = 800)

#+ fig.dim = c(5.5, 6)
volcanoplot

#'## Pathway enrichment analysis using pathfindR

# Open a file connection for writing messages
message_file <- file(paste0(output_pathfindR, "message pathfindR Reactome CAF Co-Culture vs. Monoculture.txt"), 
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
cat("begin pathfindR Reactome CAF Co-Culture vs. Monoculture \n")

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
                       sheetName = "CAF Co vs. Mono Reactome")

openxlsx::writeData(wb = wb,
                    sheet = "CAF Co vs. Mono Reactome", 
                    x = pathfindR)

# saving the workbook
openxlsx::saveWorkbook(wb = wb, 
                       file = paste0(output_pathfindR, "/output pathfindR Reactome CAF Co-Culture vs. Monoculture.xlsx"),
                       overwrite = TRUE)

cat("finished pathfindR Reactome CAF Co-Culture vs. Monoculture \n")

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
  ggplot2::scale_color_gradient(low = "thistle3", 
                                high = "firebrick") +
  ggplot2::labs(title = paste0("CAF Co-Culture vs. Monoculture \nEnrichment bubble chart pathfindR Reactome")) +
  ggplot2::theme(legend.text = element_text(size = 10),
                 legend.title = element_text(size = 12),
                 plot.title = element_text(size = 16),
                 text = element_text(family = "sans"),
                 axis.text.y = element_text(size = 15)) 

ggplot2::ggsave(filename = paste0(output_pathfindR, "enrichment_chart/enrichment chart pathfindR Reactome CAF Co-Culture vs. Monoculture.tiff"), 
                plot = enrichment_chart,
                width = 14,
                height = 8,
                dpi = 800)

#+ fig.dim = c(14, 8)
enrichment_chart

#'## session and packages
#' Here is important information about packages and session

# Open a file connection for writing messages
message_file <- file("packages CAF Co-Culture vs. Monoculture.txt", 
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
  "purrr") |>
  purrr::map(citation) |>
  print(style = "text")

java_version <- system("java -version", 
                       intern = TRUE)

print(java_version, 
      sep = "\n", 
      style = "text")

# stop redirecting messages
sink(type = "output")
sink(type = "message")

# Close the file connection
close(message_file)  
