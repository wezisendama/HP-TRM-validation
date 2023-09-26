library(dplyr)
library(DESeq2)
library(gplots)
library(R.utils)

# Download and unzip matrix metadata from GEO, then strip out the header comments and keep potentially useful fields
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150910/matrix/GSE150910_series_matrix.txt.gz",
              "series_matrix.txt.gz")

gunzip("series_matrix.txt.gz",
       "series_matrix.txt",
       overwrite = TRUE)

series_matrix <- read.csv("series_matrix.txt",
                          skip = 30,
                          sep = "\t")

series_matrix <- series_matrix[,-1]
series_matrix <- t(series_matrix)
series_matrix <- cbind(rownames(series_matrix), data.frame(series_matrix, row.names=NULL))


series.metadata <- series_matrix[,c(1,11,12,13,14,15,16,20,24)]
metadata_colnames <- c("sample",
                       "sex",
                       "age",
                       "diagnosis",
                       "ever_smoked",
                       "race",
                       "genotype",
                       "antigen_identified",
                       "institution")
colnames(series.metadata) = metadata_colnames

# Return names of samples where the institution is not recorded as "LTRC"
nonLTRCsamples <- series.metadata[series.metadata$institution != 'institution: LTRC',1]

# Download gene count matrix from GEO and strip out all the non-LTRC samples
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150910/suppl/GSE150910_gene-level_count_file.csv.gz", "gene-level_count_file.csv.gz")
gunzip("gene-level_count_file.csv.gz", "gene-level_count_file.csv", overwrite = TRUE)
main.seqdata <- read.csv("gene-level_count_file.csv", row.names = 1)
nonLTRC.seqdata <- main.seqdata[,nonLTRCsamples]

# Arrange the metadata in alphabetical order by sample name so metadata can match up with gene count matrix
series.metadata <- arrange(series.metadata, sample)

# Remove the useless bits of the strings in the metadata
series.metadata$sex <- sub(".*: ", "", series.metadata$sex)
series.metadata$age <- sub(".*: ", "", series.metadata$age)
series.metadata$diagnosis <- sub(".*: ", "", series.metadata$diagnosis)
series.metadata$ever_smoked <- sub(".*: ", "", series.metadata$ever_smoked)
series.metadata$race <- sub(".*: ", "", series.metadata$race)
series.metadata$genotype <- sub(".*: ", "", series.metadata$genotype)
series.metadata$antigen_identified <- sub(".*: ", "", series.metadata$antigen_identified)
series.metadata$institution <- sub(".*: ", "", series.metadata$institution)

# Format columns in the metadata
series.metadata$sex <- as.factor(series.metadata$sex)
series.metadata$age <- as.integer(series.metadata$age)
series.metadata$diagnosis <- as.factor(series.metadata$diagnosis)
series.metadata$diagnosis <- factor(series.metadata$diagnosis, levels=c('control', 'ipf', 'chp'))
series.metadata$ever_smoked <- as.factor(series.metadata$ever_smoked)
series.metadata$race <- as.factor(series.metadata$race)
series.metadata$genotype <- as.factor(series.metadata$genotype)
series.metadata$antigen_identified <- as.factor(series.metadata$antigen_identified)
series.metadata$institution <- as.factor(series.metadata$institution)

nonLTRC.metadata <- series.metadata[series.metadata$institution != 'LTRC',]

# Sort the count matrix by column name, then check the count matrix column names
# and metadata row names line up, answer should be TRUE
nonLTRC.seqdata <- nonLTRC.seqdata[,order(colnames(nonLTRC.seqdata))]
all(nonLTRC.metadata$sample == colnames(nonLTRC.seqdata))

# Differential gene expression analysis between diagnoses with DESeq
dds <- DESeqDataSetFromMatrix(countData = nonLTRC.seqdata,
                              colData = nonLTRC.metadata,
                              design = ~ diagnosis)
dds <- DESeq(dds)

chpres <- results(dds, contrast=c("diagnosis", "chp", "control"))
ipfres <- results(dds, contrast=c("diagnosis", "ipf", "control"))
diffres <- results(dds, contrast = c("diagnosis", "chp", "ipf"))

# Worth noting that some of the genes/probes in the cluster (CXCL10, IFNG, TRAC) identified from the discovery dataset
# do not appear (?no reads mapped) in this validation dataset, hence not in this list
upreg.genecluster <- c("CXCL9",
                       "ZNF683",
                       "CALHM6",
                       "GBP5",
                       "STAT1",
                       "JAKMIP1",
                       "FASLG",
                       "CD2",
                       "CXCR3",
                       "CD3D",
                       "LAG3",
                       "ITGB7",
                       "LGALS2")
chpres.genecluster <- chpres[upreg.genecluster,]
sig.upreg.genecluster <- rownames(chpres.genecluster[which(chpres.genecluster$log2FoldChange > 0.58 & chpres.genecluster$padj < 0.05),])

# Construct a matrix of differential expression values to print as heatmap
heatmapmatrix <- rbind(t(chpres[sig.upreg.genecluster,]$log2FoldChange), t(ipfres[sig.upreg.genecluster,]$log2FoldChange), t(diffres[sig.upreg.genecluster,]$log2FoldChange))
rownames(heatmapmatrix) = c("HP vs healthy control", "IPF vs healthy control", "HP vs IPF")
colnames(heatmapmatrix) = sig.upreg.genecluster

# Output supplementary data files showing differential expression
write.csv(chpres[sig.upreg.genecluster,], "Differential expression of cluster in HP in validation data.csv")
write.csv(ipfres[sig.upreg.genecluster,], "Differential expression of cluster in IPF in validation data.csv")

# Draw heatmap. lhei, lwid and margin values affect sizes of cells for visualisation. Good luck
heatmap.2(heatmapmatrix,
          col = greenred(10),
          density.info = "none",
          dendrogram = "column",
          Rowv = FALSE,
          trace = "none",
          keysize = 1.1,
          colsep = 1:10,
          rowsep = 1:2,
          sepwidth = c(0.001, 0.001),
          lhei = c(1.1,4), lwid = c(1,2),
          margins = c(23,17),
          cexRow = 1, cexCol = 1,
          key.title = "Log2-ratio",
          key.xlab = "Down-regulated              Up-regulated")
