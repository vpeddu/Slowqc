library("tximport")
library("tidyverse")
library("EnsDb.Hsapiens.v86")

args = commandArgs(trailingOnly=TRUE)

f_path<-args[1]
ref_file<-read.csv(args[2], sep = '\t', header = FALSE, stringsAsFactors = FALSE, col.names = c('ens_transcript_id','ens_gene_id','gene_symbol'))
ref_file$ens_gene_id<-sapply(strsplit(as.character(ref_file$ens_gene_id), "[.]"), "[[", 1)

filenames<-list.files(f_path, pattern = '*.tsv', full.names = TRUE)

txdb <- EnsDb.Hsapiens.v86
tx2gene <- transcripts(txdb, return.type="DataFrame")
tx2gene <- as.data.frame(tx2gene[,c("tx_id", "gene_id")])

txi <- tximport(as.character(filenames), type = "kallisto", tx2gene = tx2gene, ignoreTxVersion=TRUE)
df<-txi$counts

filenames_short<-list.files(f_path, pattern = "*.tsv")

colnames(df)<-filenames_short
df<-as.data.frame(df)

# Annotates by gene symbol


df$gene_id<-rownames(df)

merged<-merge(df, ref_file, by.x="gene_id", by.y="ens_gene_id") 

merged$gene_id<-NULL
merged$ens_transcript_id<-NULL

# Merges counts by gene symbol
merged_aggregated<-aggregate(merged[,c(1:ncol(merged)-1)], by=list(Category=merged$gene_symbol), FUN=sum)

colnames(merged_aggregated)[1]<-'gene'

# Writes DEBrowser ready table
write.table(merged_aggregated, file = "DEBrowser_input.txt", sep = "\t", row.names = FALSE)
