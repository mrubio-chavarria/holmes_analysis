
# Libraries
library(oligo)
library(biomaRt)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(qusage)
library(limma)


# Functions
create_pheno_description = function(metadata){
  num_samples = nrow(metadata)
  num_classes = length(unique(metadata$Experimental_group))
  
  first_line = as.character(c(num_samples, num_classes, 1))
  second_line = paste('#', paste(unique(metadata$Experimental_group), collapse=' '), collapse='')
  third_line = metadata$Experimental_group
  
  return( list(first_line, second_line, third_line))
}


# Set data directory
setwd('/home/mario/Projects/holmes_analysis/data')

# Read data
cel_files <- list.files(path = getwd(), pattern = '*.CEL.gz', full.names = TRUE)

# Set working directory
setwd('/home/mario/Projects/holmes_analysis')

# Parse CEL files
parsed_cels <- oligo::read.celfiles(cel_files, verbose = TRUE)  

# Background correction of the microarrays
parsed_cels_rma <- oligo::rma(parsed_cels, normalize = TRUE, background = TRUE)  

# Translate from probes to genes
mart <- useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl')
mouse_probes <- row.names(parsed_cels_rma@assayData$exprs)
id_translation_table <- getBM(attributes = c('affy_mouse430_2', 'ensembl_gene_id', 'mgi_symbol'),
                             filters = 'affy_mouse430_2',
                             values = mouse_probes,
                             mart=mart)
expression_data <- parsed_cels_rma@assayData$exprs
expression_data <- as.data.frame(expression_data)
expression_data$affy_mouse430_2 <- rownames(expression_data)
expression_data <- merge(x = id_translation_table, y = expression_data, by = 'affy_mouse430_2')
expression_data$affy_mouse430_2 <- NULL
expression_data <- expression_data %>% dplyr::rename(Ensembl = ensembl_gene_id)
reduced_id_translation_table <- id_translation_table
reduced_id_translation_table$affy_mouse430_2 <- NULL
reduced_id_translation_table <- reduced_id_translation_table[!duplicated(reduced_id_translation_table), ]
reduced_id_translation_table$mgi_symbol <- toupper(reduced_id_translation_table$mgi_symbol)

# Remove probes that do not match to Ensembl IDs
expression_data  <-  data.table(expression_data)
expression_data <- expression_data[!is.na(expression_data$Ensembl),]


# Compute the mean of all the probes that match to the same gene
sample_columns <- row.names(parsed_cels_rma@phenoData@data)
expression_data <- expression_data[,lapply(.SD, mean), by=Ensembl, .SDcols=sample_columns]

# Format the data columns
colnames(expression_data) <- gsub('_Mouse430v2.CEL.gz', '', colnames(expression_data))
colnames(expression_data) <- substr(colnames(expression_data), 12, 35)
colnames(expression_data)[1] <- 'Ensembl'
expression_data <- as.data.frame(expression_data)
expression_data <- merge(x = expression_data, y = reduced_id_translation_table, by.x = 'Ensembl', by.y = 'ensembl_gene_id')


################################################################################
# SINGLE-GENE ANALYSIS (log2 fold-change intensities)
################################################################################
# Compare by effect

# KO VS WILD TYPE
wt <- expression_data[, grepl('WT', colnames(expression_data))]
ko <- expression_data[, grepl('KO', colnames(expression_data))]
first <- TRUE
for (i in 1:dim(expression_data)[1]) {
 test <- t.test(wt[i, ], ko[i, ])
 row <- data.frame(Ensembl = expression_data[i, 'Ensembl'],
            p.value = test$p.value,
            mean_diff = test$estimate[[2]] - test$estimate[[1]])
 if (first) {
   first <- FALSE
   ko_vs_wt <- row
   next
 }
 ko_vs_wt <- rbind(ko_vs_wt, row)
}
# Multiple testing correction
ko_vs_wt$adjusted.p.value <- p.adjust(ko_vs_wt$p.value, method = 'BH')
# Change IDs by names
ko_vs_wt <- merge(x = ko_vs_wt, y = reduced_id_translation_table,
                  by.x = 'Ensembl', by.y = 'ensembl_gene_id')
# Print the results (volcano plot)
ko_vs_wt <- ko_vs_wt %>% dplyr::mutate(mlog10PValue = -log10(p.value))
relevants <- ko_vs_wt[ko_vs_wt$adjusted.p.value <= 0.05, ]
relevants <- relevants[order(-abs(relevants$mean_diff)), ]
relevants <- relevants[1:25, ]
ko_vs_wt %>%
  ggplot + 
  geom_point(aes(x = mean_diff, y = mlog10PValue, colour = adjusted.p.value <= 0.05), size = 4) + 
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color = "red", size=1.5) +
  scale_color_brewer(palette="Set2") + 
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  geom_text_repel(data = relevants,
                  aes(x = mean_diff, y = mlog10PValue, label = mgi_symbol)) +
  labs(x = '-log2 fold-change', y = '-log10 p-value', colour = 'Adjusted p-value <= 0.05') +
  labs(title = 'Knockout vs Wild type')

# TREATED VS UNTREATED
ut <- expression_data[, grepl('Untreated', colnames(expression_data))]
t <- expression_data[, grepl('IL-33', colnames(expression_data))]
first <- TRUE
for (i in 1:dim(expression_data)[1]) {
  test <- t.test(ut[i, ], t[i, ])
  row <- data.frame(Ensembl = expression_data[i, 'Ensembl'],
                    p.value = test$p.value,
                    mean_diff = test$estimate[[2]] - test$estimate[[1]])
  if (first) {
    first <- FALSE
    t_vs_ut <- row
    next
  }
  t_vs_ut <- rbind(t_vs_ut, row)
}
# Multiple testing correction
t_vs_ut$adjusted.p.value <- p.adjust(t_vs_ut$p.value, method = 'BH')
# Change IDs by names
t_vs_ut <- merge(x = t_vs_ut, y = reduced_id_translation_table,
                  by.x = 'Ensembl', by.y = 'ensembl_gene_id')
t_vs_ut <- t_vs_ut[!duplicated(t_vs_ut), ]
# Print the results (volcano plot)
t_vs_ut <- t_vs_ut %>% dplyr::mutate(mlog10PValue = -log10(p.value))
relevants <- t_vs_ut[t_vs_ut$adjusted.p.value <= 0.05, ]
relevants <- relevants[order(-abs(relevants$mean_diff)), ]
relevants <- relevants[1:25, ]
t_vs_ut %>%
  ggplot + 
  geom_point(aes(x = mean_diff, y = mlog10PValue, colour = adjusted.p.value <= 0.05), size = 4) + 
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color = "red", size=1.5) +
  scale_color_brewer(palette="Set2") + 
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  geom_text_repel(data = relevants,
                  aes(x = mean_diff, y = mlog10PValue, label = mgi_symbol)) +
  labs(x = '-log2 fold-change', y = '-log10 p-value', colour = 'Adjusted p-value <= 0.05') +
  labs(title = 'Treated vs Untreated')


################################################################################
# HYPERGEOMETRIC TEST (TREATED VS UNTREATED)
################################################################################
# Compare by phenotype (effect)

# Read gene sets
gene_sets <- read.gmt('gene_sets/biocarta_gene_sets.gmt')
expression_data$mgi_symbol <- toupper(expression_data$mgi_symbol)

# Translate genes to entrez ids
id_translation_table_2 <- getBM(attributes = c('mgi_symbol', 'entrezgene_id'),
                                filters = 'mgi_symbol',
                                values = expression_data$mgi_symbol,
                                mart=mart)
id_translation_table_2$mgi_symbol <- toupper(id_translation_table_2$mgi_symbol)

# Perform the hypergeometric test per gene set
altered_t_vs_ut <- t_vs_ut
differentially_expressed_genes <- altered_t_vs_ut[altered_t_vs_ut$adjusted.p.value <= 0.05, ]$mgi_symbol

# Perform the hypergeometric test
# Genes in the arrays
N <- length(expression_data$Ensembl)
# Number of differentiated genes
n <- length(differentially_expressed_genes)  
# Test p-values
hyper.p.values <- c()  
# Number of genes in the set
n_genes_set <- c()
# Number of differentially expressed genes in the set
n_genes_in_the_set <- c()
for (gene_set in gene_sets) {
  # Number of differentially expressed genes in the set
  x <- sum(differentially_expressed_genes %in% gene_set)
  # Number of genes in the set
  k <- length(unlist(gene_set))
  # Compute the test
  p.value <- phyper(x, k, N - k, n, lower.tail = FALSE)
  hyper.p.values <- c(hyper.p.values, p.value)
  n_genes_set <- c(n_genes_set, k)
  n_genes_in_the_set <- c(n_genes_in_the_set, x)
}
hyper.p.values <- p.adjust(hyper.p.values, 'BH')
hyper_results <- data.frame(gene.set = names(gene_sets),
                            adjusted.p.value = hyper.p.values,
                            n.set = n_genes_set,
                            n.in.set = n_genes_in_the_set)
relevant_hyper_results <- hyper_results[hyper_results$adjusted.p.value <= 0.05, ]


################################################################################
# PREPARE FILE FOR GSEA
################################################################################

# RANKED (TREATED VS UNTREATED)
# Prepare the input
t_vs_unt_gsea_input <- expression_data[, c(22, 2:21)]
t_vs_unt_gsea_input <- cbind(t_vs_unt_gsea_input[, 1],
                             rep("https://www.google.com", each=dim(expression_data)[1]),
                             t_vs_unt_gsea_input[, 2:21])
colnames(t_vs_unt_gsea_input)[1:2] <- c("NAME", "DESCRIPTION")
t_vs_unt_gsea_input <- t_vs_unt_gsea_input[t_vs_unt_gsea_input$NAME != '', ]
write.table(t_vs_unt_gsea_input, 'gsea_input.txt', sep='\t', quote = F, row.names = F)
# Prepare the labels
experimental_group <- substr(colnames(t_vs_unt_gsea_input)[3:22], 15, 23)
num_samples = length(experimental_group)
num_classes = 2  # Treated and untreated
first_line <- as.character(c(num_samples, num_classes, 1))
second_line <- paste(c('#', unique(experimental_group)), collapse=' ')
third_line <- experimental_group
labels_file <- file('gsea_labels.cls')
open(labels_file, 'w')
for(line in list(first_line, second_line, third_line)){
  line = paste(line, collapse=' ')
  writeLines(line, con = labels_file)
}
close(labels_file)




