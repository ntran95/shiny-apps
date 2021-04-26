library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
server <- FALSE
if(server){
devtools::load_all("/n/projects/nt2473/Analysis/Scripts/SeuratExtensions")
devtools::load_all("/n/projects/nt2473/Analysis/Scripts/CellTrajectoryExtensions")


dataPath <- function(filename){paste0("/n/projects/",
                                      "nt2473/Analysis/Data/sb2191-regen/", 
                                      filename)}

  }else{
    devtools::load_all("/Volumes/projects/nt2473/Analysis/Scripts/SeuratExtensions")
    devtools::load_all("/Volumes/projects/nt2473/Analysis/Scripts/CellTrajectoryExtensions")
    
    
    dataPath <- function(filename){paste0("/Volumes/projects/",
                                          "nt2473/Analysis/Data/sb2191-regen/", 
                                          filename)}
}
seurat_obj <- readRDS(dataPath(paste0(
  "SeurObj_ptime_subset_experiments_seurat3_10hr_central_v1.0_dim_6_.RDS")))

cds <- readRDS(dataPath(paste0("saved_obj_", 
              "HC_lineage_homeo_10hr_monocle3_10_central_v1.0", "_.RDS")))

cell_type_trt <- levels(cds$cell.type.and.trt)
type_trt_cols <- gg_color_hue(length(cell_type_trt))

traject_plot <- plot_cells(cds, color_cells_by = "cell.type.and.trt",
                           label_leaves = TRUE, 
                           show_trajectory_graph = TRUE, 
                           cell_size = 0.70,
                           trajectory_graph_segment_size = 0.55, 
                           label_cell_groups = FALSE,
                           group_label_size = 8) + 
  theme(legend.position="bottom")

traject_plot <- cleanUMAP(traject_plot)
traject_plot <- traject_plot + scale_color_manual(values = type_trt_cols)

traject_plot

ptime_plot <- plot_cells(cds, color_cells_by = "pseudotime",
                         label_leaves = TRUE, 
                         show_trajectory_graph = TRUE, 
                         cell_size = 0.80,
                         trajectory_graph_segment_size = 0.60, 
                         label_cell_groups = FALSE)

ptime_plot <- ptime_plot +
  guides(fill = guide_legend(ncol = 1))
ptime_plot <- cleanUMAP(ptime_plot)
ptime_plot
# ============================================ overlap traj line with UMAP coordinate
x <- 1
y <- 2
ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
  as.data.frame() %>%
  dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))

dp_mst <- cds@principal_graph[["UMAP"]]

edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select_(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       source="sample_name",
                       source_prin_graph_dim_1="prin_graph_dim_1",
                       source_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       target="sample_name",
                       target_prin_graph_dim_1="prin_graph_dim_1",
                       target_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "target")

mst_branch_nodes <- branch_nodes(cds, reduction_method ="UMAP")
branch_point_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
  dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))

mst_leaf_nodes <- leaf_nodes(cds, reduction_method = "UMAP")
leaf_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
  dplyr::mutate(leaf_idx = seq_len(dplyr::n()))

mst_root_nodes <- root_nodes(cds, reduction_method = "UMAP")
root_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
  dplyr::mutate(root_idx = seq_len(dplyr::n()))



features <- c("atoh1a", "her4.1")

trajectory_graph_segment_size=0.75
graph_label_size=2

g <- FeaturePlot(seurat_obj, features = features, slot = "data")
g <- cleanUMAP(g)
g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
                                 y="source_prin_graph_dim_2",
                                 xend="target_prin_graph_dim_1",
                                 yend="target_prin_graph_dim_2"),
                      linetype="solid",
                      na.rm=TRUE,
                      data=edge_df) 
g <- g + #plot branching
  geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
             shape = 21, stroke=I(0.75),
             color="white",
             fill="black",
             size=I(2 * 1.5),
             na.rm=TRUE, branch_point_df) +
  geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                       label="branch_point_idx"),
            size=I(2), color="white", na.rm=TRUE,
            branch_point_df)

g <- g + #plot leaves
  geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
             shape = 21, stroke=I(0.75),
             color="black",
             fill="lightgray",
             size=I(2 * 1.5),
             na.rm=TRUE,
             leaf_df) +
  geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                       label="leaf_idx"),
            size=I(2), color="black", na.rm=TRUE, leaf_df)

g <- g + #plot root node
  geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
             shape = 21, stroke=I(0.75),
             color="black",
             fill="white",
             size=I(2 * 1.5),
             na.rm=TRUE,
             root_df) +
  geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                       label="root_idx"),
            size=I(2), color="black", na.rm=TRUE, root_df)
g

# ======================================================== line dynamic individual genes
gene <- c("atoh1a")

plot_cells(ptime_central_traj)
plot_cells(ptime_main_traj)

make_plot_df <- function(cds_sub, gene){
  cds_sub <- cds_sub[rownames(cds_sub) %in% gene,]
  cds_exprs <- as.matrix(SingleCellExperiment::counts(cds_sub))
  
  count_mtx_sub <- cds_exprs
  #count_mtx_sub <- count_mtx_sub[!apply(count_mtx_sub,1,sum)==0,]
  count_mtx_sub <- t(count_mtx_sub)
  count_mtx_sub <- scale(log1p(count_mtx_sub))
  count_mtx_sub <- as.data.frame(Seurat::MinMax(count_mtx_sub , min = -2.5, 
                                                max = 2.5))
  mod1_df <- as.data.frame(count_mtx_sub)
  mod1_df$Cell <- rownames(mod1_df)
  
  mod1_df <- reshape2::melt(mod1_df)
  colnames(mod1_df)[2:3] <- c("Gene.name.uniq","expression")
  
  
  ptime_df <- data.frame(pseudotime = t(pseudotime(cds))[1,],
                         Cell = names(pseudotime(cds)),
                         cell_group = colData(cds)$cell.type.and.trt)
  
  #order cells in correspondence to ptime
  mod1_df$Cell <- factor(mod1_df$Cell, levels = ptime_df$Cell)
  
  plot_dt <- inner_join(mod1_df, ptime_df)
  
  
  return(plot_dt)
}

central_traj_df <- make_plot_df(cds_sub = ptime_central_traj,gene = gene)
nrow(central_traj_df)
main_HC_traj_df <- make_plot_df(cds_sub = ptime_main_traj, gene = gene)
nrow(main_HC_traj_df)

central_traj_geom_smooth_col <- "#00BE67"

main_traj_geom_smooth_col <- "#F8766D"

# color_indx <- which(cell_type_trt %in% levels(cds@colData$cell.type.and.trt))
# select_colors <- type_trt_cols[color_indx]
cell_type_trt <- levels(seurat_obj$cell.type.and.trt)
type_trt_cols <- gg_color_hue(length(cell_type_trt))

test_lg <- ggplot(central_traj_df) +
  geom_smooth(
    colour=central_traj_geom_smooth_col,
    span=0.2,
    method='loess',
    fill = central_traj_geom_smooth_col,
    aes( x=pseudotime, y=expression,fill="Central Cell Lineage"),
    fullrange = TRUE)

test_lg <- test_lg +
  geom_smooth(data = main_HC_traj_df,
              colour=main_traj_geom_smooth_col,
              span=0.2,
              method='loess',
              fill = main_traj_geom_smooth_col,
              aes( x=pseudotime, y=expression, fill="HC Lineage"))


test_lg <- test_lg + geom_rug(data=main_HC_traj_df, sides='b', 
                              alpha=.10, aes(x=pseudotime,
                                             color = cell_group) )

test_lg <- test_lg + geom_rug(data=central_traj_df, sides='t', alpha=.10, 
                              aes(x=pseudotime,
                                 color = cell_group) )

test_lg <- test_lg + scale_color_manual(values = type_trt_cols)
test_lg <-   test_lg + theme_bw() +
  theme(
  ) +
  labs(
    x     = 'pseudotime'
    ,y    = 'z-scored gene expression'
  )
test_lg <- test_lg +theme(legend.position="bottom", legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

#add facet title to 
test_lg$data$title <- unique(test_lg$data$Gene.name.uniq)
test_lg <- test_lg + facet_wrap(~title)
test_lg <- test_lg+ 
  theme(strip.text.x = element_text(size = 18))
test_lg

#========================== multiple gene line plot
library(ggnewscale)
smpl_genes_lg <- paste0("atoh1a her4.1 hes2.2 dld sox4a*1")
selected <- unique(unlist(strsplit(smpl_genes_lg, " ")))

make_plot_df <- function(cds_sub, gene){
  cds_sub <- cds_sub[rownames(cds_sub) %in% gene,]
  cds_exprs <- as.matrix(SingleCellExperiment::counts(cds_sub))
  
  count_mtx_sub <- cds_exprs
  #count_mtx_sub <- count_mtx_sub[!apply(count_mtx_sub,1,sum)==0,]
  count_mtx_sub <- t(count_mtx_sub)
  count_mtx_sub <- scale(log1p(count_mtx_sub))
  count_mtx_sub <- as.data.frame(Seurat::MinMax(count_mtx_sub , min = -2.5, 
                                                max = 2.5))
  mod1_df <- as.data.frame(count_mtx_sub)
  mod1_df$Cell <- rownames(mod1_df)
  
  mod1_df <- reshape2::melt(mod1_df)
  colnames(mod1_df)[2:3] <- c("Gene.name.uniq","expression")
  
  
  ptime_df <- data.frame(pseudotime = t(pseudotime(cds))[1,],
                         Cell = names(pseudotime(cds)),
                         cell_group = colData(cds)$cell.type.and.trt)
  
  #order cells in correspondence to ptime
  mod1_df$Cell <- factor(mod1_df$Cell, levels = ptime_df$Cell)
  
  plot_dt <- inner_join(mod1_df, ptime_df)
  
  
  return(plot_dt)
}

central_traj_df <- make_plot_df(cds_sub = ptime_central_traj,gene = selected)
#order genes in order of input sequence 
central_traj_df$Gene.name.uniq <- factor(central_traj_df$Gene.name.uniq, 
                                         levels = selected)
#specify trajectory trail
central_traj_df $trajectory <- "Central Cell Lineage"
nrow(central_traj_df)
main_HC_traj_df <- make_plot_df(cds_sub = ptime_main_traj, gene = selected)
#specify trajectory trail
main_HC_traj_df $trajectory <- "HC Lineage"
main_HC_traj_df$Gene.name.uniq <- factor(main_HC_traj_df$Gene.name.uniq, 
                                         levels = selected)
nrow(main_HC_traj_df)

#cancatenate two df 
plot_dt <- rbind(central_traj_df,main_HC_traj_df)
p <- ggplot(data = plot_dt,
       mapping = aes(x = pseudotime, y = expression, color = Gene.name.uniq,
                     fill = Gene.name.uniq)) +
  geom_smooth(method ="loess", span = 0.2, fullrange = TRUE) +
  new_scale_color() + # add new scale color for geom_rug
  new_scale_fill() +# add new scale color for geom_rug
  geom_rug(data=plot_dt[plot_dt$trajectory == "Central Cell Lineage",], sides='b', 
           alpha=.10, aes(x=pseudotime, color = cell_group) ) +
  scale_color_manual(values = type_trt_cols, guide = "none") +
  new_scale_color() + # add new scale color for geom_rug
  new_scale_fill() +# add new scale color for geom_rug
  geom_rug(data=plot_dt[plot_dt$trajectory == "HC Lineage",], sides='b', 
           alpha=.10, aes(x=pseudotime, color = cell_group) ) +
  scale_color_manual(values = type_trt_cols) + theme_bw() +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme(
  ) +
  labs(
    x     = 'pseudotime'
    ,y    = 'z-scored gene expression'
  ) + 
  theme(legend.position="bottom", 
            legend.box = "vertical",
            legend.title=element_blank()) +
  facet_wrap(~ trajectory) + #split plot by trajectory path
  theme(strip.text.x = element_text(size = 18, face = "bold")) +#specify facet title size 
  coord_cartesian(ylim=c(min(ggplot_build(p)$data[[1]]$y,na.rm = TRUE), #dynamically plot min  and max y lim
                         max(ggplot_build(p)$data[[1]]$y,na.rm = TRUE)))

p

#================== manually apply facet panel colors 
g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))

pal <- c(central_traj_geom_smooth_col, main_traj_geom_smooth_col)

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i] #change facet  background label colors
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- "white" #change facet letter colors
}

plot(g)

# ============== add vertical line at branching point

#get embedding info
embeddings <- as.data.frame(seurat_obj[["umap"]]@cell.embeddings)
embeddings$Cell <- rownames(embeddings)
#round
embeddings[,1:2] <- round(embeddings[,1:2],digits = 2)

#get branching info from monocle princle graph
branch_point_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
  dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))

#identify branching point bt central + HC Lineage, #1
branch_point_df <- branch_point_df %>% filter(branch_point_idx == "1")
#round
branch_point_df[,1:2] <- round(branch_point_df[,1:2], digits = 2)

#retrieve pseudotime, cell, cell group info
ptime_df <- data.frame(pseudotime = t(pseudotime(cds))[1,],
                       Cell = names(pseudotime(cds)),
                       cell_group = colData(cds)$cell.type.and.trt)
#join embedding info with pseudotime info
ptime_df <- inner_join(ptime_df, embeddings)

branching_embedding <- ptime_df %>% filter(abs(UMAP_1 - as.numeric(branch_point_df$prin_graph_dim_1)) == min(abs(UMAP_1 - as.numeric(branch_point_df$prin_graph_dim_1)))) %>%
  filter(abs(UMAP_2 - as.numeric(branch_point_df$prin_graph_dim_2)) == min(abs(UMAP_2 - as.numeric(branch_point_df$prin_graph_dim_2))))

#extract first option, retrieve only the pseudotime val
branching_ptime <- branching_embedding[1,]

get_branching_point <- function(seurat_obj, cds){
  x <- 1
  y <- 2
  ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
    as.data.frame() %>%
    dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    dplyr::mutate(sample_name = rownames(.),
                  sample_state = rownames(.))
  
  #get embedding info
  embeddings <- as.data.frame(seurat_obj[["umap"]]@cell.embeddings)
  embeddings$Cell <- rownames(embeddings)
  #round
  embeddings[,1:2] <- round(embeddings[,1:2],digits = 2)
  
  #get branching info from monocle princle graph
  branch_point_df <- ica_space_df %>%
    dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
    dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
  
  #identify branching point bt central + HC Lineage, #1
  branch_point_df <- branch_point_df %>% filter(branch_point_idx == "1")
  #round
  branch_point_df[,1:2] <- round(branch_point_df[,1:2], digits = 2)
  
  #retrieve pseudotime, cell, cell group info
  ptime_df <- data.frame(pseudotime = t(pseudotime(cds))[1,],
                         Cell = names(pseudotime(cds)),
                         cell_group = colData(cds)$cell.type.and.trt)
  #join embedding info with pseudotime info
  ptime_df <- inner_join(ptime_df, embeddings)
  
  branching_embedding <- ptime_df %>% filter(abs(UMAP_1 - as.numeric(branch_point_df$prin_graph_dim_1)) == min(abs(UMAP_1 - as.numeric(branch_point_df$prin_graph_dim_1)))) %>%
    filter(abs(UMAP_2 - as.numeric(branch_point_df$prin_graph_dim_2)) == min(abs(UMAP_2 - as.numeric(branch_point_df$prin_graph_dim_2))))
  
  #extract first option, retrieve only the pseudotime val
  branching_ptime <- branching_embedding[1,]
  
  return(branching_ptime)
}

branching_ptime <- get_branching_point(seurat_obj = seurat_obj, cds = cds)

b <- DimPlot(seurat_obj, cells.highlight = "3hr_TGCGGGTGTAAGAGAG", cols.highlight = "blue") +
  scale_color_manual(labels = c("unselected", "branching"), values = c("grey", "blue"))  + NoAxes()
b <- b+ #plot branching
  geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
             shape = 21, stroke=I(0.75),
             color="white",
             fill="black",
             size=I(2 * 1.5),
             na.rm=TRUE, branch_point_df) +
  geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                       label="branch_point_idx"),
            size=I(2), color="white", na.rm=TRUE,
            branch_point_df)

#add smoothing gene line
p <- ggplot(data = plot_dt,
            mapping = aes(x = pseudotime, y = expression, color = Gene.name.uniq,
                          fill = Gene.name.uniq)) +
  geom_smooth(method ="loess", span = 0.2, fullrange = TRUE) 

#add vertical line at branching point
p <- p +
  new_scale_color() + # add new scale color for geom_rug
  new_scale_fill() +# add new scale color for geom_rug
  geom_vline(data = branching_ptime, aes(xintercept = pseudotime, color = "cell_group"),
             linetype=2) +
  scale_color_manual(name = " ",  values = "red",
                     labels = c("branching point"))

#add geom_rug for dashed cell types
p <- p+
  new_scale_color() + # add new scale color for geom_rug
  new_scale_fill() +# add new scale color for geom_rug
  geom_rug(data=plot_dt[plot_dt$trajectory == "Central Cell Lineage",], sides='b', 
           alpha=.10, aes(x=pseudotime, color = cell_group) ) +
  scale_color_manual(values = type_trt_cols, guide = "none") +
  new_scale_color() + # add new scale color for geom_rug
  new_scale_fill() +# add new scale color for geom_rug
  geom_rug(data=plot_dt[plot_dt$trajectory == "HC Lineage",], sides='b', 
           alpha=.10, aes(x=pseudotime, color = cell_group) ) +
  scale_color_manual(values = type_trt_cols) 

#legend themes
p <- p + theme_bw() +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme(
  ) +
  labs(
    x     = 'pseudotime'
    ,y    = 'z-scored gene expression'
  ) + 
  theme(legend.position="bottom", 
        legend.box = "vertical",
        legend.title=element_blank()) +
  facet_wrap(~ trajectory) + #split plot by trajectory path
  theme(strip.text.x = element_text(size = 18, face = "bold")) +#specify facet title size 
  coord_cartesian(ylim=c(min(ggplot_build(p)$data[[1]]$y,na.rm = TRUE), #dynamically plot min  and max y lim
                         max(ggplot_build(p)$data[[1]]$y,na.rm = TRUE)))


p


lg <- ggplot(central_traj_df) +
  geom_smooth(
    #colour=central_traj_geom_smooth_col,
    span=0.2,
    method='loess',
    #fill = central_traj_geom_smooth_col,
    aes( x=pseudotime, y=expression,color="Central Cell Lineage", fill = "Central Cell Lineage"),
    fullrange = TRUE)

lg <- lg +
  geom_smooth(data = main_HC_traj_df,
              #colour=main_traj_geom_smooth_col,
              span=0.2,
              method='loess',
              #fill = main_traj_geom_smooth_col,
              aes( x=pseudotime, y=expression, color="HC Lineage", fill = "HC Lineage"))

lg <- lg + scale_color_manual(name="Branching Trajectories",guide = 'legend',
                                                values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
                                                           "HC Lineage" = main_traj_geom_smooth_col)) + guides(guide_legend(override.aes = list(linetype = c("black","black"))))
lg <- lg + scale_fill_manual(name="Branching Trajectories",guide = 'legend',
                                               values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
                                                         "HC Lineage" = main_traj_geom_smooth_col))
lg <- lg  +
  new_scale_color() + # add new scale color for branching point
  new_scale_fill() # add new scale color for branching point

# add branching vertical line
lg <- lg +
  new_scale_color() + # add new scale color for geom_rug
  new_scale_fill() +# add new scale color for geom_rug
  geom_vline(data = branching_ptime, 
             aes(xintercept = pseudotime, color = "cell_group"),
             linetype=2) +
  scale_color_manual(name = " ",  values = "red",
                     labels = c("branching point"))

lg <- lg  +
  new_scale_color() + # add new scale color for geom_rug
  new_scale_fill() # add new scale color for geom_rug
  
lg <- lg + geom_rug(data=main_HC_traj_df, sides='b',
                    alpha=.10, aes(x=pseudotime,
                                   color = cell_group) )

lg <- lg + geom_rug(data=central_traj_df, sides='t', alpha=.10,
                    aes(x=pseudotime,
                        color = cell_group) )

lg <- lg + scale_color_manual(values = type_trt_cols)
lg <-   lg + theme_bw() +
  theme(
  ) +
  labs(
    x     = 'pseudotime'
    ,y    = 'scaled gene expression'
  )
lg <- lg +theme(legend.position="bottom",
                legend.box = "vertical",
                legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

#add facet title to
lg$data$title <- unique(lg$data$Gene.name.uniq)
lg <- lg + facet_wrap(~title)
lg <- lg+
  theme(strip.text.x = element_text(size = 18))


lg

