library("tximport")
library("tidyr")


args = commandArgs(trailingOnly=TRUE)

f_path<-args[1]
ref_file<-read.csv(args[2], sep = '\t', header = FALSE, stringsAsFactors = FALSE, col.names = c('ens_transcript_id','ens_gene_id','gene_symbol'))
ref_file$ens_gene_id<-sapply(strsplit(as.character(ref_file$ens_gene_id), "[.]"), "[[", 1)

filenames<-list.files(f_path, pattern = '*.tsv', full.names = TRUE)

#txdb <- EnsDb.Hsapiens.v86
#tx2gene <- transcripts(txdb, return.type="DataFrame")
#tx2gene <- as.data.frame(tx2gene[,c("tx_id", "gene_id")])

#txi <- tximport(as.character(filenames), type = "kallisto", tx2gene = tx2gene, ignoreTxVersion=TRUE)

read_reports<-function(csv_list){ 
  for(i in 1:length(csv_list)){ 
    if(i == 1){ 
        temp<-read.csv(csv_list[i], header = TRUE, stringsAsFactors = FALSE, sep = '\t')
        df<-temp[,c(1,4)]
        colnames(df)[ncol(df)]<-basename(strsplit(filenames[i], '_abundance')[[1]][1])
        df[,ncol(df)]<-as.numeric(as.character(df[,ncol(df)]))
    }
    else{ 
      temp<-read.csv(csv_list[i], header = TRUE, stringsAsFactors = FALSE, sep = '\t')
      df<-cbind(df,temp[,4])
      colnames(df)[ncol(df)]<-basename(strsplit(filenames[i], '_abundance')[[1]][1])
      df[,ncol(df)]<-as.numeric(as.character(df[,ncol(df)]))
      
      }
  }
  return(df)
  }

#df<-txi$counts

#filenames_short<-list.files(f_path, pattern = "*.tsv")

#colnames(df)<-filenames_short
df<-read_reports(filenames)

#Annotates by gene symbol





merged<-merge(df, ref_file, by.x="target_id", by.y="ens_transcript_id") 

merged$target_id<-NULL
merged$ens_transcript_id<-NULL
merged$ens_gene_id<-NULL

# Merges counts by gene symbol
merged_aggregated<-aggregate(merged[,c(2:ncol(merged)-1)], by=list(Category=merged$gene_symbol), FUN=sum)

colnames(merged_aggregated)[1]<-'gene'

# Writes DEBrowser ready table
write.table(merged_aggregated, file = "DEBrowser_input.txt", sep = "\t", row.names = FALSE)
