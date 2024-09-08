

# load the required libraries
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(curatedTCGAData)
library(pheatmap)
library(SummarizedExperiment)
library(TCGAbiolinks)
# get a list of projects
gdcprojects <- getGDCprojects()
gdcprojects$disease_type
summarypro <- getProjectSummary('TCGA-PRAD')

#tcga_ids <- read.csv("/home/admin/ids.csv",stringsAsFactors = FALSE)

######## 1 make an TCGA object #######
# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-PRAD',data.category = 'Transcriptome Profiling',experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts', access = 'open')
GDCquery()

# download data - GDCdownload
####### 2 download the object ########
GDCdownload(query_TCGA, directory = "/home/admin/Sung-prostate-cancer/")

####### 3 prepare data #########
tcga_prost_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE,directory = "/home/admin/Sung-prostate-cancer/")

######## 4 get the assay ########
prostate_matrix <- assay(tcga_prost_data, 'stranded_first')

prostate_matrix <- as.data.frame(prostate_matrix)

######## 5 convert Ensemble Ids to gene symbols ########

# Install and load org.Hs.eg.db package
library(org.Hs.eg.db)

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(prostate_matrix)

# Remove version numbers for compatibility
ensembl_ids <- sub("\\..*", "", ensembl_ids)

# Map IDs
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Remove NA values and keep only mapped genes
mapped_genes <- !is.na(gene_symbols)
prostate_matrix <- prostate_matrix[mapped_genes, ]
gene_symbols <- gene_symbols[mapped_genes]

unique_gene_symbols <- make.unique(gene_symbols, sep = "_")
length(unique_gene_symbols)

# Assign gene symbols to prostate_matrix
rownames(prostate_matrix) <- unique_gene_symbols
# Clean up gene symbols
valid_symbols <- !is.na(rownames(prostate_matrix)) & 
  !grepl("^[0-9.]+$", rownames(prostate_matrix)) &
  !grepl("^LOC", rownames(prostate_matrix)) &
  grepl("^[A-Za-z0-9-]+$", rownames(prostate_matrix))

prostate_matrix <- prostate_matrix[valid_symbols, ]
prostate_matrix <-   prostate_matrix[, grep("01A|11A", colnames(prostate_matrix)), drop = FALSE]
dim(prostate_matrix)
# Get column names
column_names <- colnames(prostate_matrix)

# Count occurrences of '01A'
count_01A <- sum(str_count(column_names, "01A")) # 483 primary tumor samples  

# Count occurrences of '11A'
count_11A <- sum(str_count(column_names, "11A")) # 51 normal samples

# Identify columns with '01A' and '11A'
tumor_columns <- column_names[str_detect(column_names, "01A")]
normal_columns <- column_names[str_detect(column_names, "11A")]

# Create a new data frame with labeled samples
tumor_samples <- prostate_matrix[, tumor_columns]
normal_samples <- prostate_matrix[, normal_columns]

# Label the columns
colnames(tumor_samples) <- paste0(colnames(tumor_samples), "_tumor")
colnames(normal_samples) <- paste0(colnames(normal_samples), "_normal")

# Combine the labeled samples into one data frame
labeled_prostate_matrix <- cbind(tumor_samples, normal_samples)


write.csv(labeled_prostate_matrix, "labeled_prostate_matrix.csv")
write.csv(prostate_matrix, "prostate_matrix.csv")
############### Deseq2 Normalizations ####################
tail(labeled_prostate_matrix)
# my phenodata 
phdata <- colnames(labeled_prostate_matrix)
phdata <- as.data.frame(phdata)

# Assuming your data frame is called 'df'
phdata$disease_state <- ifelse(grepl("01A", phdata$phdata), "tumor", 
                               ifelse(grepl("11A",  phdata$phdata), "normal", NA))

names(phdata)[1] <- "ids"
rownames(phdata) <- phdata$ids

# making the rownames and column names identical
all(rownames(phdata) %in% colnames(labeled_prostate_matrix))
all(rownames(phdata) == colnames(labeled_prostate_matrix))

labeled_prostate_matrix <- round(labeled_prostate_matrix)
dds <- DESeqDataSetFromMatrix(countData = labeled_prostate_matrix,
                              colData = phdata,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds) #  2158 genes


# perform variance stabilization
dds_norm <- varianceStabilizingTransformation(dds)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

# norm.counts <- norm.counts[tumor$ids,]
norm.counts <- t(norm.counts) %>% as.data.frame()

write.csv(norm.counts, "normalized_labeled_prostate_data.csv")





























