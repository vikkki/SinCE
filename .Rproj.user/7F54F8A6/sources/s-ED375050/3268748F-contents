
#### static file paths ###
path_to_demo_data_folder <- file.path(system.file("SinCEapp", package = "SinCE"))

#### -- data filetration parameter in loading matrix to seurat objects
min_cells <- 3 # genes expressed in >= 3 cells (~0.1% of the data, will be dynamic in the future)
min_features <- 200 # keep all cells with at least 200 detected genes

#### PCA ####
n_pca = 150

#### color set ####
my_color = c("#0B90AA","#7dce94","#B1B336","#04384A","#66638B","#D74B4B","#FF652D","#F6AE2D","#AE8D65","#D8DEAE","#70AB8F")
#
# -> show color in bar plot:
# df <- data.frame(x = c(1:11), y = rep(1,11))
# ggplot(df) +
#   geom_bar(aes(x = x, y = y), stat="identity",
#            fill = my_color)
