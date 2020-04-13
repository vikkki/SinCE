library(Seurat)
library(Matrix)


########################################
# loading functions
########################################

load_10x_3_files <- function(barcode_path,feature_path,matrix_path){
  library(Matrix)
  mat <- readMM(file = matrix_path)
  barcode.names <- read.delim(barcode_path, 
                              header = FALSE,
                              stringsAsFactors = FALSE)
  feature.names <- read.delim(feature_path, 
                              header = FALSE,
                              stringsAsFactors = FALSE)
  
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V2
  
  return(mat*1)
}


load_demo_file <- function(demo_path){
  library(Matrix)
  df <- read.csv(demo_path,
                 header = TRUE,
                 sep = ";",
                 row.names = NULL)
  # -- remove duplicated feature(gene) names
  dup = duplicated(df$X) 
  df <- df[!dup,]
  rownames(df)<-df$X
  df = df[,-1]
  
  mat <- as.matrix(df,rownames = T)
  
  return(mat)
}

load_csv_file <- function(csv_path,sep){
  library(Matrix)
  df <- as.data.frame(read.csv2(csv_path,
                 header = TRUE,
                 sep = sep,
                 row.names = NULL))
  # -- remove duplicated feature(gene) names
  dup = duplicated(df$X)
  df <- df[!dup,]
  rownames(df)<-df$X
  df = df[,-1]
  
  as.matrix(df,rownames = T)
}

####################################
# -- update the cluster with identity table: row names are barcode, column ident is identity(factor)
ident_update = function(sc,ident_table){
  barcodes = row.names(ident_table)
  withProgress(message="Processing ...",{
  for(i in 1:length(barcodes)){
    Idents(sc, cells = barcodes[i]) <- ident_table[barcodes[i],"ident"]
  }
  })
  sc
}