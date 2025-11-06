########------infercnv
library(tidyverse)
library(Seurat)
library(infercnv)
###-----1. Data loading
epi_so <- readRDS('Epithelial_cells_sample.rds')
ref_so <- readRDS('Ref_cells_sample_for_inferCNV.rds')

so <- merge(epi_so,ref_so)
so$group <- ifelse(so$cluster_standard %in% c('Cancer cells'),
                   so$sample,
                   so$cluster_standard
)

# ----2. Expression matrix preparation
raw_counts_matrix <- so@assays[["RNA"]]@counts

# ----3. Cell annotation file preparation
df_meta <- so@meta.data
cell.annotation <- data.frame(group = so$group,
                              row.names = rownames(so@meta.data))
write.table(cell.annotation, paste(getwd(), "cell.annotation.txt", sep = "/"),
            col.names = FALSE, sep = "\t",quote = F)

# ----4. Gene order file (chromosomal order)
all.genes <- c(rownames(so))

df_gene <- as.data.frame(table(all.genes), row.names = NULL)
colnames(df_gene) <- c("hgnc_symbol", "freq")
head(df_gene,10)

gene_ref <- read.table('gencode_v19_gene_pos.txt',stringsAsFactors = F,header = F)
names(gene_ref) <- c("hgnc_symbol",'chromosome_name',
                     'start_position', 'end_position')
gene_po <- inner_join(df_gene,gene_ref,by = 'hgnc_symbol')
gene_po$freq <- NULL
df_gene$hgnc_symbol %>% unique() %>% length()
write.table(gene_po, paste(getwd(), "gene_chromopos.txt", sep = "/"),
            col.names = FALSE, row.names = FALSE,  sep = "\t", quote = FALSE)
# filter the counts matrix according to results of chromosome positions
counts_matrix <- raw_counts_matrix[c(gene_po$hgnc_symbol), ]


# -----5. infercnv
out_dir <- paste(getwd(), "/InferCNV_sample/", sep = "")
if (dir.exists(out_dir)){
  out_dir <- out_dir
} else {dir.create(out_dir)}

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=paste(getwd(), "cell.annotation.txt", sep = "/"),
                                    delim="\t",
                                    gene_order_file= paste(getwd(), "gene_chromopos.txt", sep = "/"),
                                    ref_group_names=c("Club cells", 
                                                      "Alveolar epithelial cells", "Ciliated cells", "Endothelial cells", 
                                                      "Fibroblasts", "Mural cells")) 


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)

print(paste("files are written in:", out_dir, sep = " "))
