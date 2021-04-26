# ==== concatenate original lovain clustering n= 10 with n=19 to allow users to ocillate between two options ====

seurat_obj_cl10 <- 
  readRDS(paste0("/Volumes/projects/nt2473/Analysis/Scripts/",
  "Piotrowski-Lab/shiny-apps/pseudotime_regen_traj_Sungmin/data/",
"TRIMMED_SeurObj_ptime_subset_experiments_seurat3_10hr_central_v1.0_dim_6_.RDS"))

seurat_obj_c19 <- 
  readRDS(paste0("/Volumes/projects/nt2473/Analysis/Data/sb2191-regen/",
"SeurObj_ptime_subset_experiments_seurat3_10hr_central_v1.1_dim_6_res_2_.RDS"))

#rename tree.ident column
meta_cl10 <- seurat_obj_cl10@meta.data
#colnames(meta_cl10)[colnames(meta_cl10) =="tree.ident"] <- "tree.ident.cl10"
meta_cl10$tree.ident.cl10 <- meta_cl10$tree.ident

#concatenate tree.ident to original shiny app seurat obj metadata
meta_cl19 <- seurat_obj_c19@meta.data
meta_cl10$tree.ident.cl19 <- meta_cl19$tree.ident

#remove unused metadata columns
meta_cl10 <- meta_cl10[,!(grepl("snn_res|seurat_clusters", colnames(meta_cl10)))]


seurat_obj_cl10@meta.data <- meta_cl10
seurat_obj_cl10$tree.ident.cl10 <- factor(seurat_obj_cl10$tree.ident.cl10)
seurat_obj_cl10$tree.ident.cl19 <- factor(seurat_obj_cl10$tree.ident.cl19)
#check
Idents(seurat_obj_cl10) <- "tree.ident.cl19"
DimPlot(seurat_obj_cl10)

Idents(seurat_obj_cl10) <- "tree.ident.cl10"
DimPlot(seurat_obj_cl10)

#save
setwd(paste0("/Volumes/projects/nt2473/Analysis/Scripts/",
  "Piotrowski-Lab/shiny-apps/pseudotime_regen_traj_Sungmin/data/"))
saveRDS(seurat_obj_cl10, 
        file = "TRIMMED_SeurObj_ptime_subset_experiments_seurat3_10hr_central_v1.0_dim_6_.RDS")
