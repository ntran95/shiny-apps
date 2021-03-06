files <- list.files("./", pattern = "TRIMMED", full.names = TRUE)

Idents(file_list[[1]]) <- file_list[[1]]@meta.data$tree.ident
seurat_obj <- file_list[[1]]

gene_df <- read.table("./Danio_Features_unique_Ens98_v1.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

gene_df <- gene_df[gene_df$Gene.name.uniq %in% rownames(seurat_obj),]

colnames(gene_df)
gene_df <- gene_df[,c(1:3,6,4:5)]
gene_df$ZFIN.ID <- paste0("=HYPERLINK(", '"', gene_df$ZFIN.ID, '"',")")

write.table(gene_df, "./Danio_Features_unique_Ens98_v1.1.tsv", sep = "\t")