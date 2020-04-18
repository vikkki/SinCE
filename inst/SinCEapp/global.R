
#### static file paths ###
path_to_demo_data_folder <- file.path(system.file("SinCEapp", package = "SinCE"))

#### -- data filetration parameter in loading matrix to seurat objects
min_cells <- 3 # genes expressed in >= 3 cells (~0.1% of the data, will be dynamic in the future)
min_features <- 200 # keep all cells with at least 200 detected genes

#### largest PC number ####
n_pca = 150

#### color set ####
my_color = c("#0B90AA","#7dce94","#B1B336","#04384A","#66638B","#D74B4B","#FF652D","#F6AE2D","#AE8D65","#D8DEAE","#70AB8F")
#
# -> show color in bar plot:
# df <- data.frame(x = c(1:11), y = rep(1,11))
# ggplot(df) +
#   geom_bar(aes(x = x, y = y), stat="identity",
#            fill = my_color)
par_templete <- list(
  "clean_feature" = c(200,3500),
  "clean_count" = c(500, 30000),
  "clean_mito" = c(0, 15),

  "seurat_nomalize_method" = "LogNormalize",
  "seurat_nomalize_scale_factor" = 10000,
  "seurat_nomalize_margin" = "Features",

  "seurat_cluster_pc" = 15,
  "seurat_cluster_resolution" = 0.205,
  "plot_cluster_feature" = "CD79A,S100A9",
  "cor_feature_x" = "CD79A",
  "cor_feature_y" = "S100A9",

  "included_cell_number" = 0,
  "cluster_list" = "1",

  "seurat_tsne_run_method" = "Rtsne",
  "seurat_tsne_pc" = 15,
  "seurat_tsne_seed" = 12,
  "seurat_tsne_max_iter" = 1000,

  "seurat_umap_run_method" = "umap-learn",
  "umap_max_pc" = 15,
  "umap_learning_rate" = 1,
  "umap_n_neighbor" = 10,
  "umap_min_dist" = 0.3,
  "umap_spread" = 1,
  "seurat_umap_seed" = 12
)

