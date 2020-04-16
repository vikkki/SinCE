
library(shiny, verbose=FALSE)
library(ggplot2, verbose=FALSE)
library(Seurat, verbose=FALSE)
library(RColorBrewer, verbose=FALSE)
library(plotly)

## parameters here ##


options(shiny.maxRequestSize=800*1024^2) # limit upload file to 700M
project_name <- "sc_project"
set.seed(4422)

## functions here ##
source("sc_functions.R")
getPalette = colorRampPalette(my_color) # use color palette my_color, in global.R


shinyServer(function(input, output, clientData, session) {

  ## dynamic user input ui ##
  output$load_alart_ui <- renderText("Beeeeee")
  ## load data ##

  # output$load_alart_ui <- reactive({
  #   input$load_sample_button
  #   r = "Ready to start."
  #   isolate({
  #     r =  paste0("Ready to start.",input$load_sample_button)})
  #   return(r)})

  demo_name <- reactive({switch(input$demo_select,
                                "demo_2.5k_pbmc" = "PBMC1_nova_2500_filtered",
                                "demo_1k_pbmc" = "10x_pbmc_1k_v3_filtered",
                                "demo_5k_pbmc" = "10x_pbmc_5k_v3_filtered",
                                "demo_10k_pbmc" = "10x_pbmc_10k_v3_filtered"
                                #,"mini_sample" = "sample"
  )
  })

  #### load data to matrix based on input####
  sc_matrix <- reactive({
    input$load_sample_button
    isolate({withProgress(message="Loading data ...",{
      output$load_alart_ui <- renderText("Ready to start.")
      if(input$load_sample_button == 0){output$load_alart <- renderText("No data loaded. Ready to start loading.")
      return(NULL)}
      else {
        if(is.null(input$csv_matrix_file)){
          if(is.null(input$barcods_file) || is.null(input$features_file) || is.null(input$matrix_file)){
            # no csv or tsv or 10x matrix uploaded
            # -- path parameter set before in the head
            demo_file_path = paste0(path_to_demo_data_folder, "/", demo_name(), ".csv")

            # -- function in sc_functions.R
            mat = load_demo_file(demo_file_path)
            output$load_alart_ui <- renderText("Demo data loaded.")
            mat
          }else{
            # only 10x matrix uploaded
            tryCatch({mat = load_10x_3_files(input$barcods_file$datapath, input$features_file$datapath, input$matrix_file$datapath)},
                     error = function(e){output$create_mat_alart <- renderText("Loading 10x matrix error: please double check the input file(s).")
                     output$load_alart_ui <- renderText("10x data loaded failed.")
                     mat = NULL}
            )
            if(!is.null(mat)) output$load_alart_ui <- renderText("10x data loaded.")
            mat
          }
        }else {
          if(is.null(input$barcods_file) || is.null(input$features_file) || is.null(input$matrix_file)){
            # only the csv uploaded
            sep <- input$csv_seperator
            mat = tryCatch({load_csv_file(input$csv_matrix_file$datapath,sep)},
                           error = function(e){output$create_mat_alart <- renderText("Loading csv matrix error: please double check the input file.")
                           output$load_alart_ui <- renderText("Csv data load failed.")
                           NULL})
            if(!is.null(mat)){

              if(dim(mat)[2] == 1){
                output$create_mat_alart <- renderText("Loading csv matrix error 1 please double check the input file(for example the seperator).")
                output$load_alart_ui <- renderText("Csv data load failed.")
                mat=NULL
              }else {
                output$create_mat_alart <- renderText("csv matrix loaded.")
                output$load_alart_ui <- renderText("csv data loaded.")
              }

            }else {
              output$create_mat_alart <- renderText("Loading csv matrix error 2: please double check the input file(for example the seperator).")
              output$load_alart_ui <- renderText("csv data load failed.")
              return(NULL)}

            mat
          }else{
            output$load_alart_ui <- renderText("Too many inputs, please hint reset and start over to have one of the data types uploaded.")
            return(NULL)

          }

        }
      }
    })})
  })

  #### download demo data ####

  output$download_demo_data <- downloadHandler(
    filename <- function(){paste0(demo_name(),".csv")},
    content <- function(file){
      download_demo_path = paste0(path_to_demo_data_folder,"/",demo_name(),".csv")
      file.copy(download_demo_path,file)},
    contentType = "text/csv"

  )

  #### prepare the seurat object####
  sc_seurat_base <- reactive({
    # -- the object will be (re)generate if submit button clicked or matrix chenged.
    if(input$load_sample_button == 0 || is.null(sc_matrix())) return(NULL)
    withProgress(message="Creating sc object ...",{
      # -- initialize the Seurat object
      tryCatch({sc_base = CreateSeuratObject(counts = sc_matrix(), project = project_name, min.cells = min_cells, min.features = min_features)},
               error = function(e){
                 output$create_seurat_alart <- renderText("Something wrong in loading scRNA data, please double check the input file(s)")
                 return(NULL)})
      # -- chech if 'nCount_RNA' exist in seurat object
      sc_valid <- tryCatch(sc_base$nCount_RNA[1],
                           error = function(e){NULL})
      if(is.null(sc_valid)) {
        output$create_seurat_alart <- renderText("Something wrong in loading scRNA data, please double check the input file(s)")
        return(NULL)
      }
      sc_base
    })
  })

  output$seurat_load_message <- renderText({
    if(is.null(sc_seurat_base())) return(NULL)
    sc = sc_seurat_base()
    dim = dim(sc)
    paste0(dim[1], " features and ", dim[2], " barcodes(cells) loaded.")})

  output$load_summary <- renderTable({
    mat = sc_matrix()
    mat[2:5,1:5]
  },rownames = TRUE, colnames = TRUE)
  output$path<- renderText(input$csv_matrix_file$datapath)

  #############################################
  # test zone
  #############################################


  #############################################
  # test zone end
  #############################################

  ## view data ##

  #### QC  ####
  sc_seurat_qc_sample <- reactive({
    # -- the object will be (re)generate if submit button clicked or sc_seurat_base chenged.
    if(input$load_sample_button == 0 || is.null(sc_seurat_base())) return(NULL)
    withProgress(message="Processing QC info ...",{
      sc_base <- sc_seurat_base()
      # --The [[ operator can add columns to object metadata. This is a great place to stash QC stats
      sc_base[["percent.mt"]] <- PercentageFeatureSet(sc_base, pattern = "^MT-")
      sc_base
    })
  })
  # -- visualize QC metrics as a violin plot
  output$sample_qc_violin <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_qc_sample())) return(NULL)
    VlnPlot(sc_seurat_qc_sample(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  })

  output$cluster_qc_violin <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL)
    sc = sc_seurat_cluster()
    sc = ident_update(sc,seurat_cluster_ident_table())
    VlnPlot(sc_seurat_cluster(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  })


  # -- visualize QC features
  output$sample_qc_feature <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_qc_sample())) return(NULL)
    plot1 <- FeatureScatter(sc_seurat_qc_sample(), feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(sc_seurat_qc_sample(), feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

    CombinePlots(plots = list(plot1, plot2))
  })


  #### cleaning data ####

  # -- slider UI update
  observe(label = "qc_sliders",{
    if(input$load_sample_button == 0 || is.null(sc_seurat_qc_sample())) return(NULL)
    updateSliderInput(session, "clean_feature",
                      max = max(sc_seurat_qc_sample()[["nFeature_RNA"]])
    )
    updateSliderInput(session, "clean_count",
                      max = max(sc_seurat_qc_sample()[["nCount_RNA"]])
    )
  })


  sc_seurat_filted <- reactive({
    tem = input$clean_button
    if(input$load_sample_button == 0 || is.null(sc_seurat_qc_sample())) return(NULL)
    withProgress(message="Cleaning data ...",{
      # -- now we have a hard coded threshold for cleaning data, and in the future there will be a dynamic way based either on data or user input
      if(input$clean_button == 0)
        subset(sc_seurat_qc_sample(), subset = nFeature_RNA >= 200 & nFeature_RNA <= 3500 & percent.mt <= 15) # if the clean data button isn't clicked, default filter will be apply
      else
        isolate({
          sc = sc_seurat_qc_sample()
          fmin = input$clean_feature[1]
          fmax = input$clean_feature[2]
          cmin = input$clean_count[1]
          cmax = input$clean_count[2]
          mmin = input$clean_mito[1]
          mmax = input$clean_mito[2]

          cells_use <- colnames(sc)[which(sc[[]]['nFeature_RNA'] >= fmin & sc[[]]['nFeature_RNA'] <= fmax &
                                            sc[[]]['nCount_RNA'] >= cmin & sc[[]]['nCount_RNA'] <= cmax &
                                            sc[[]]['percent.mt'] >= mmin & sc[[]]['percent.mt'] <= mmax)]
          subset(sc, cells = cells_use)
        })
    })
  })
  output$cleaned_cell_number <- renderText(paste0(dim(sc_seurat_filted())[2]," cell barcodes selected."))

  output$cleaned_qc_violin <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_filted())) return(NULL)
    VlnPlot(sc_seurat_filted(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  })

  output$cleaned_qc_feature <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_filted())) return(NULL)
    plot1 <- FeatureScatter(sc_seurat_filted(), feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(sc_seurat_filted(), feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

    CombinePlots(plots = list(plot1, plot2))
  })



  #### normalizing data ####

  sc_seurat_normalized <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_filted())) return(NULL)
    withProgress(message="Normalizing data ...",{
      # -- employ a global-scaling normalization method “LogNormalize” that normalizes the gene expression measurements for each cell by the total expression,
      # -- multiplies this by a scale factor (10,000 by default), and log-transforms the result.
      NormalizeData(sc_seurat_filted(), normalization.method = "LogNormalize",
                    scale.factor = 10000)
    })
  })

  #### identification of highly variable features ####
  sc_seurat_variable <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_normalized())) return(NULL)
    withProgress(message="Identifing highly variable features ...",{
      # -- defualt sets now
      FindVariableFeatures(sc_seurat_normalized(), selection.method = "vst", nfeatures = 2500)
    })
  })


  seurat_topn_variable_features <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_variable())) return(NULL)
    head(VariableFeatures(sc_seurat_variable()), input$topn)
  })

  output$seurat_topn_variable_features <- renderText({
    if(input$load_sample_button == 0 || is.null(sc_seurat_variable())) return(NULL)
    seurat_topn_variable_features()
  })

  output$seurat_variable_feature_plot <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_variable())) return(NULL)
    # -- plot variable features with and without labels
    plot1 <- VariableFeaturePlot(sc_seurat_variable())
    plot2 <- LabelPoints(plot = plot1, points = seurat_topn_variable_features(), repel = TRUE)
    return(plot2)
  })

  output$seurat_variable_density_plot <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_variable())) return(NULL)
    # -- plot variable features with and without labels
    hvf.info =HVFInfo(object = sc_seurat_variable()[["RNA"]], selection.method = 'vst',status = TRUE)
    var.status <- c('no', 'yes')[unlist(x = hvf.info[, ncol(x = hvf.info)]) + 1] # convert TRUE/FALSE to "yes" and "no"
    hvf.info$is_Variable <- factor(var.status)
    ggplot(hvf.info,aes(x = mean, fill = is_Variable))+
      geom_density(alpha = 0.3)+ scale_x_log10()+
      xlab("Average Expression")+
      scale_fill_manual(values=c("black","red"))+
      theme_set(theme_bw())+
      theme(panel.grid.major=element_line(colour=NA))+
      theme(panel.border = element_blank())
  })

  output$seurat_variable_count_plot <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_variable())) return(NULL)
    # -- plot variable features with and without labels
    hvf.info =HVFInfo(object = sc_seurat_variable()[["RNA"]], selection.method = 'vst',status = TRUE)
    var.status <- c('no', 'yes')[unlist(x = hvf.info[, ncol(x = hvf.info)]) + 1] # convert TRUE/FALSE to "yes" and "no"
    hvf.info$is_Variable <- factor(var.status)
    ggplot(hvf.info,aes(x = mean, fill = is_Variable))+
      geom_histogram(position = "identity",alpha = 0.4)+
      scale_x_log10()+
      xlab("Average Expression")+
      scale_fill_manual(values=c("black","red"))+
      theme_set(theme_bw())+
      theme(panel.grid.major=element_line(colour=NA))+
      theme(panel.border = element_blank())
  })

  #### scaling ####
  sc_seurat_scaled <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_variable())) return(NULL)
    withProgress(message="Scaling ...",{
      # -- defualt sets now,eurat constructs linear models to predict gene expression based on user-defined variables.
      sc <- sc_seurat_variable()
      all.genes <- rownames(sc)
      ScaleData(sc, features = all.genes)
    })
  })

  ####  run PCA ####
  sc_seurat_PCA <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_scaled())) return(NULL)
    withProgress(message="Running linear dimensional reduction...",{
      sc = sc_seurat_scaled()
      RunPCA(sc, npcs=n_pca, verbose = FALSE)
    })
  })

  output$seuret_pca_plot <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # depend on cluster
    sc = sc_seurat_cluster()
    sc = ident_update(sc,seurat_cluster_ident_table())
    DimPlot(sc,
            pt.size = input$pca_point_size,
            reduction = "pca",
            label = input$pca_label,
            dims = c(input$pcx,input$pcy),
            cols = getPalette(length(levels(sc_seurat_cluster()))))
  })

  sc_seurat_pca_inter_base <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL)
    sc = sc_seurat_cluster()
    sc = ident_update(sc,seurat_cluster_ident_table())
    ident = as.data.frame(sc@active.ident)
    colnames(ident) <- "ident"
    embeds = as.data.frame(Embeddings(sc[["pca"]]),col.names = T)
    embeds = cbind.data.frame(embeds, ident)
    embeds$nCount_RNA  = sc@meta.data[["nCount_RNA"]]
    embeds$nFeature_RNA = sc@meta.data[["nFeature_RNA"]]
    embeds$cluster <- sc@meta.data[["seurat_clusters"]]
    embeds$keys <- rownames(embeds)
    return(embeds)
  })

  pca_inter <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_pca_inter_base())) return(NULL)
    dim =paste0("PC_",c(input$pcx,input$pcy))
    ax <- list(
      zeroline = FALSE,
      title = dim[1]
      # ,gridcolor = #bdc3c7,
      # gridwidth = 0.6
    )
    ay <- list(
      zeroline = FALSE,
      title = dim[2]
      # ,gridcolor = #bdc3c7,
      # gridwidth = 0.6
    )
    sc = sc_seurat_pca_inter_base()
    plot_ly(sc,type="scatter", mode = "markers",
            x = sc[,input$pcx], y = sc[,input$pcy],
            color = ~ident,
            colors = getPalette(length(levels(sc$ident))),
            hoverinfo = "all",
            hovertext = paste0(sc$keys,"\n","nCount:",sc$nCount_RNA,"\n","nFeature:",sc$nFeature_RNA)) %>% layout(dragmode = "lasso", xaxis = ax, yaxis = ay)

  })

  output$seuret_pca_plot_inter <- renderPlotly({pca_inter()})

  pcn <- reactive(c(input$pcn[1]:input$pcn[2]))

  output$seurat_pca_genes <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_PCA())) return(NULL)
    VizDimLoadings(sc_seurat_PCA(), dims = pcn(), nfeatures = input$pc_genes_n, reduction = "pca", combine = TRUE)
  })

  output$seurat_pca_heatmap<- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_PCA())) return(NULL)
    DimHeatmap(sc_seurat_PCA(), dims = pcn(), cells = input$pca_heatmap_gene, nfeatures = input$pc_genes_n, reduction = "pca", balanced = TRUE)
  })

  output$seurat_pca_score <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_PCA())) return(NULL)
    temp = input$pca_score_button
    isolate({
      if(input$pca_score_button >0)
        withProgress(message="Counting for PCs ...",{
          sc <- JackStraw(sc_seurat_PCA(), num.replicate = 100)
          sc <- ScoreJackStraw(sc, dims = 1:15)
          JackStrawPlot(sc, dims = 1:15)
        })
    })

  })


  ####  cluster the cells ####

  cluster_modi_counter <- reactiveVal(0)

  # -- initail cluster name table
  seurat_cluster_ident_table <- reactiveVal(NULL)

  cluster_names <- reactiveVal(NULL)

  sc_seurat_cluster <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_PCA())) return(NULL)
    temp = input$seurat_cluster_resolution
    temp = input$seurat_cluster_pc
    isolate({
      withProgress(message="Looking for clusters ...",{
        cluster_modi_counter(0) #initialize counter when cluster regenerated
        output$name_message = NULL
        output$rename_selected_message = NULL
        sc <- sc_seurat_PCA()
        sc <- FindNeighbors(sc, dims = 1:input$seurat_cluster_pc)
        sc <- FindClusters(sc, resolution = input$seurat_cluster_resolution)
        # sc_seurat_cluster(sc)
        ident <- sc@active.ident
        ident <- as.data.frame(ident)
        names = levels(sc)
        cluster_names(names)
        seurat_cluster_ident_table(ident)
        return(sc)
      })
    })
  })


  cluster_names_raw <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on cluster
    levels(sc_seurat_cluster())
  })

  output$cluster_names_raw <- renderText(cluster_names_raw())

  # -- convert user input names
  user_cluster_names <- reactive({
    strsplit(input$input_cluster_names,",")[[1]]
  })

  # -- update ident table when name_button clicked
  observeEvent(input$name_button,
               {
                 if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) ident_table = NULL
                 else{
                   if(cluster_modi_counter() == 0){
                     # cluster haven't been synconized to current raw cluster result
                     ident <- sc_seurat_cluster()@active.ident
                     ident_table <- as.data.frame(ident)
                     cluster_modi_counter(1)
                   }
                   else{ # name have been changed
                     current_names <- levels(seurat_cluster_ident_table()$ident)
                     if(length(unique(current_names)) == length(user_cluster_names())){
                       new.cluster.ids <- user_cluster_names()
                       names(new.cluster.ids) <- current_names
                       ident_table <- seurat_cluster_ident_table()
                       ident_table$ident <- as.character(ident_table$ident)
                       for(i in 1: dim(ident_table)[1])
                         ident_table$ident[i] <- new.cluster.ids[ident_table$ident[i]]

                       ident_table$ident <- as.factor(ident_table$ident)
                       output$name_message = renderText("Cluster identification changed")
                       cluster_modi_counter(1)
                     }
                     else {
                       ident_table = seurat_cluster_ident_table()
                       output$name_message = renderText("Inputed names invalid")
                     }
                   }
                 }
                 names = levels(ident_table$ident)
                 cluster_names(names)
                 seurat_cluster_ident_table(ident_table)
               },ignoreNULL = TRUE,ignoreInit = TRUE)

  output$cluster_names <- renderText({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) NULL
    cluster_names()})

  #-- group2 UI update, only show clusters which haven't been choosen in group 1
  observeEvent(cluster_names(),{
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster()) || is.null(cluster_names())) return(NULL)
    updateCheckboxGroupInput(session, "cluster_group1",
                             choices = cluster_names(),
                             select = cluster_names()[1]
    )
  },label = "cluster_group1_observer")
  #
  # -- remain clusters after selecting group1
  remain = reactive({
    if(input$load_sample_button == 0 || is.null(input$cluster_group1)) return(NULL)
    setdiff(cluster_names(),input$cluster_group1)

  })

  observeEvent(remain(),{
    if(input$load_sample_button == 0 || is.null(remain()) || is.null(cluster_names())) return(NULL)
    updateCheckboxGroupInput(session, "cluster_group2",
                             choices = remain(),
                             select = remain()[1]
    )
  },label = "cluster_group2_observer")

  cluster_marker_filename <- reactiveVal("Cluster_makers.csv")

  seurat_cross_cluster_markers <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL)
    temp = input$find_markers_button
    isolate({
      if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL)
      if(input$find_markers_button == 0) return(NULL)
      group1 = input$cluster_group1
      group2 = input$cluster_group2
      g1 = paste(input$cluster_group1,collapse = ",")
      g2 = paste(input$cluster_group2,collapse = ",")
      filename = paste0("Cluster_", g1, "_vs_", g2, ".csv", sep ="")
      cluster_marker_filename(filename) # update marker download file name
      sc <- sc_seurat_cluster()
      sc = ident_update(sc,seurat_cluster_ident_table())
      withProgress(message="Looking for markers ...",{
        FindMarkers(
          sc,
          ident.1 = group1,
          ident.2 = group2,
          only.pos = FALSE,
          logfc.threshold = 0.25,
          min.pct = 0.15
        )
      })
    })
  })

  seurat_cross_cluster_markers_show <- reactive({
    temp = input$find_markers_button
    temp = input$topn_cluster
    isolate({
      if(input$load_sample_button == 0 || is.null(seurat_cross_cluster_markers())) return(NULL)
      if(input$find_markers_button == 0) return(NULL)
      seurat_cross_cluster_markers()

    })
  })

  output$seurat_cross_cluster_markers_show <- DT::renderDataTable({
    temp = input$find_markers_button
    temp = input$topn_cluster
    DT::datatable(head(seurat_cross_cluster_markers_show(), input$topn_cluster))
  })

  output$dl_cluster_markers <- downloadHandler(
    filename = function(){
      name = cluster_marker_filename()
      return(name)},
    content = function(file) {
      markers = seurat_cross_cluster_markers_show()
      write.table(markers,file, quote = TRUE, sep = ",", col.names = TRUE, row.names = TRUE)
    }
  )


  tsne_for_brush <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_tsne_inter_base())) return(NULL)
    # -- axits settings
    ax <- list(
      zeroline = FALSE
      # gridcolor = #bdc3c7,
      # gridwidth = 1
    )
    sc = sc_seurat_tsne_inter_base()

    plot_ly(sc,type="scatter", mode = "markers",
            x = ~tSNE_1, y = ~tSNE_2,
            key = ~keys,
            color = ~ident,
            colors = getPalette(length(levels(sc$ident))),
            hoverinfo = "all",
            hovertext = paste0(sc$keys,"\n","nCount:",sc$nCount_RNA,"\n","nFeature:",sc$nFeature_RNA)) %>% layout(dragmode = "lasso", xaxis = ax, yaxis = ax)

  })

  output$tsne_for_brush <- renderPlotly(tsne_for_brush())

  cluster_brush_cells <- reactiveVal(NULL)

  output$cluster_brush_info <- renderText({
    d <- event_data("plotly_selected")
    if(is.null(d)) "Click and drag to select cells (double-click to clear)."
    else {
      cluster_brush_cells(d$key)
      paste0(nrow(d)," barcodes have been selected.")
    }
  })

  # -- rename selected cells in ident table
  observeEvent(input$rename_selected_cell_button,
               {
                 if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) ident_table = NULL
                 else{
                   if(is.null(cluster_brush_cells())) ident_table = seurat_cluster_ident_table()
                   else {
                     ident_table = seurat_cluster_ident_table()
                     ident_table$ident <- as.character(ident_table$ident)
                     ident_table$ident[row.names(ident_table) %in% cluster_brush_cells()] <- input$selected_cell_name
                     ident_table$ident <- as.factor(ident_table$ident)
                     output$rename_selected_message <- renderText("Cells identity renamed")
                     cluster_modi_counter(1) # switch counter to edited value
                   }
                 }
                 names = levels(ident_table$ident)
                 cluster_names(names)
                 seurat_cluster_ident_table(ident_table)
               },ignoreNULL = TRUE,ignoreInit = TRUE)

  # -- rename selected cells ident

  output$dl_select_cells <- downloadHandler(
    filename = function() {
      paste0("Barcode_",input$selected_cell_name, ".csv")
    },
    content = function(file) {
      barcodes = cluster_brush_cells()
      write.table(barcodes,file, quote = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
    }
  )

  #### gene expression####

  # -- convert user input genes
  plot_features <- reactive({
    # if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL)
    strsplit(input$plot_cluster_feature,",")[[1]]
  })

  output$cluster_feature <- renderText(plot_features())

  # build plot bases to enable the seperation of identity and input change
  seurat_feature_plot_pca_base <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on umap, pca or t-sne
    withProgress(message="Generating feature plots ...",{
      sc = sc_seurat_PCA()
      ident_update(sc,seurat_cluster_ident_table())
    })
  })

  seurat_feature_plot_tsne_base <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on umap, pca or t-sne
    withProgress(message="Generating feature plots ...",{
      sc = sc_seurat_tsne()
      ident_update(sc,seurat_cluster_ident_table())
    })
  })

  seurat_feature_plot_umap_base <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on umap, pca or t-sne
    withProgress(message="Generating feature plots ...",{
      sc = sc_seurat_umap()
      ident_update(sc,seurat_cluster_ident_table())
    })
  })

  output$seurat_cluster_feature_plot <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on umap, pca or t-sne
    withProgress(message="Generating feature plots ...",{
      switch (input$featureplot_reduction,
              "umap" = {
                FeaturePlot(seurat_feature_plot_umap_base(), features = plot_features(), reduction = "umap")},
              "pca" = {
                FeaturePlot(seurat_feature_plot_pca_base(), features = plot_features(), reduction = "pca")},
              "tsne" = {
                FeaturePlot(seurat_feature_plot_tsne_base(), features = plot_features(), reduction = "tsne")}
      )
    })
  })
  output$seurat_feature_inter <- renderPlotly({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on umap, pca or t-sne
    withProgress(message="Generating interactive feature plot ...",{
      switch (input$featureplot_reduction,
              "umap" = {
                FeaturePlot(seurat_feature_plot_umap_base(), features = plot_features(), reduction = "umap")},
              "pca" = {
                FeaturePlot(seurat_feature_plot_pca_base(), features = plot_features(), reduction = "pca")},
              "tsne" = {
                FeaturePlot(seurat_feature_plot_tsne_base(), features = plot_features(), reduction = "tsne")}
      )
    })
  })


  output$seurat_cluster_feature_vln <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on cluster
    sc = sc_seurat_cluster()
    sc = ident_update(sc,seurat_cluster_ident_table())
    VlnPlot(sc, features = plot_features())
  })

  output$seurat_cluster_ridge <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on cluster
    sc = sc_seurat_cluster()
    sc = ident_update(sc,seurat_cluster_ident_table())
    RidgePlot(sc, features = plot_features(), ncol = 3)
  })

  output$seurat_cluster_dot <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on cluster
    sc = sc_seurat_cluster()
    sc = ident_update(sc,seurat_cluster_ident_table())
    DotPlot(sc, features = plot_features(),cols = c("lightgrey", "darkgreen")) + RotatedAxis()
  })

  output$seurat_feature_scatter <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on cluster
    sc = sc_seurat_cluster()
    sc = ident_update(sc,seurat_cluster_ident_table())
    FeatureScatter(sc,
                   feature1 = input$featurex,
                   feature2 = input$featurey,
                   cols = getPalette(length(levels(sc_seurat_cluster()))))
  })

  output$seurat_feature_coex <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on cluster
    sc = sc_seurat_cluster()
    sc = ident_update(sc,seurat_cluster_ident_table())
    FeaturePlot(sc, features = c(input$featurex, input$featurey), blend = TRUE)
  })

  #### heatmap ####

  # -- selector UI update based on cluster_names()
  observe(label = "heatmap_cluster_select",{
    if(input$load_sample_button == 0 || is.null(sc_seurat_qc_sample())) return(NULL)
    updateCheckboxGroupInput(session, "heatmap_clusters",
                      choices = cluster_names(),
                      selected = cluster_names())
  })

  seurat_all_markers <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL) # based on cluster
    sc = sc_seurat_cluster()
    sc = ident_update(sc,seurat_cluster_ident_table())
    withProgress(message="Looking for markers (could take minutes) ...",{
      FindAllMarkers(sc, only.pos = FALSE, min.pct = 0.15, logfc.threshold = 0.25)
    })
  })

  seurat_heatmap_base <- reactive({
    if(input$load_sample_button == 0 || is.null(seurat_all_markers())) return(NULL) # based on all makers
    sc = sc_seurat_cluster()
    sc = ident_update(sc,seurat_cluster_ident_table())
    cells_use = WhichCells(sc, idents = input$heatmap_clusters)
    sc_subset <- subset(sc, cells = cells_use)
    topn <- seurat_all_markers() %>% group_by(cluster) %>% dplyr::top_n(n = input$topn_heatmap, wt = avg_logFC)
    withProgress(message="Generating heatmap ...",{
    DoHeatmap(sc_subset, features = topn$gene, slot = input$heatmap_slot) + NoLegend() + scale_fill_gradientn(colors = c("#062170", "white", "#880045"))
    })
  })

  output$seurat_marksers_heatmap <- renderPlot({
    if(input$load_sample_button == 0 || is.null(seurat_all_markers())) return(NULL) # based on all makers
    seurat_heatmap_base()
  })



  #### t-SNE ####

  sc_seurat_tsne <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL)
    # -- load in acv to make this function reactive to these variables
    acv = input$suerat_tsne_run_method
    acv = input$seurat_tsne_max_iter
    acv = input$seurat_tsne_max_pc
    # -- end acv part
    isolate({
      method = input$suerat_tsne_run_mathod
      maxiter = input$seurat_tsne_max_iter
      max_pc = input$seurat_tsne_max_pc
      withProgress(message="Running t-SNE ...",{
        # -- defualt method is rtsne
        RunTSNE(sc_seurat_cluster(), dims = 1:20, method = method, nthreads = 4, max_iter = maxiter)
      })
    })
  })

  output$seuret_tsne_plot <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_tsne())) return(NULL)
    sc = sc_seurat_tsne()
    sc = ident_update(sc,seurat_cluster_ident_table())
    DimPlot(sc,
            reduction = "tsne",
            pt.size = input$tsne_point_size,
            label = input$tsne_label,
            cols = getPalette(length(levels(sc)))) + ggtitle(label = "t-SNE")
  })

  output$seuret_tsne_plot_2 <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_tsne())) return(NULL)
    sc = sc_seurat_tsne()
    sc = ident_update(sc,seurat_cluster_ident_table())
    DimPlot(sc,
            reduction = "tsne",
            pt.size = 1,
            label = TRUE,
            cols = getPalette(length(levels(sc)))) + ggtitle(label = "t-SNE")
  })


  # -- prepare inter active tsne plot, extraction data from S4 object and convert it to dataframe
  sc_seurat_tsne_inter_base <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_tsne())) return(NULL)
    sc = sc_seurat_tsne()
    sc = ident_update(sc,seurat_cluster_ident_table())
    ident = as.data.frame(sc@active.ident)
    colnames(ident) <- "ident"
    embeds = as.data.frame(Embeddings(sc[["tsne"]]),col.names = T)
    embeds = cbind.data.frame(embeds, ident)
    embeds$nCount_RNA  = sc@meta.data[["nCount_RNA"]]
    embeds$nFeature_RNA = sc@meta.data[["nFeature_RNA"]]
    embeds$cluster <- sc@meta.data[["seurat_clusters"]]
    embeds$keys <- rownames(embeds)

    return(embeds)
  })

  output$seuret_tsne_plot_inter <- renderPlotly(tsne_for_brush())

  #### UMAP ####

  sc_seurat_umap <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster()) || is.null(sc_seurat_cluster())) return(NULL)
    withProgress(message="Running UMAP ...",{
      # -- default
      RunUMAP(sc_seurat_cluster(),
              dims = 1:20,
              learning.rate = input$umap_learning_rate,
              spread = input$umap_spread,
              min.dist = input$umap_min_dist,
              verbose = FALSE)
    })
  })

  output$seurat_umap_plot <- renderPlot({
    if(input$load_sample_button == 0 || is.null(sc_seurat_umap())) return(NULL)
    sc = sc_seurat_umap()
    sc = ident_update(sc,seurat_cluster_ident_table())
    DimPlot(sc_seurat_umap(), reduction = "umap",
            pt.size = input$umap_point_size,
            label = input$umap_label,
            cols = getPalette(length(levels(sc_seurat_umap()))))+ ggtitle(label = "UMAP")
  })

  # output$seurat_umap_plot_inter_1 <- renderPlotly({
  #   if(input$load_sample_button == 0 || is.null(sc_seurat_umap())) return(NULL)
  #   sc = sc_seurat_umap()
  #   plot = DimPlot(sc, reduction = "umap",
  #           pt.size = input$umap_point_size,
  #           label = input$umap_label,
  #           cols = getPalette(length(levels(sc_seurat_umap()))))+ ggtitle(label = "UMAP")
  #   HoverLocator(plot = plot, information = FetchData(sc, vars = c("ident","UMAP_1","UMAP_2","nFeature_RNA","nCount_RNA")))
  # })

  # -- prepare inter active tsne plot, extraction data from S4 object and convert it to dataframe
  sc_seurat_umap_inter_base <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_umap())) return(NULL)
    sc = sc_seurat_umap()
    sc = ident_update(sc,seurat_cluster_ident_table())
    ident = as.data.frame(sc@active.ident)
    colnames(ident) <- "ident"
    embeds = as.data.frame(Embeddings(sc[["umap"]]),col.names = T)
    embeds = cbind.data.frame(embeds, ident)
    embeds$nCount_RNA  = sc@meta.data[["nCount_RNA"]]
    embeds$nFeature_RNA = sc@meta.data[["nFeature_RNA"]]
    embeds$cluster <- sc@meta.data[["seurat_clusters"]]
    embeds$keys <- rownames(embeds)

    return(embeds)
  })

  umap_inter <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_umap_inter_base())) return(NULL)
    # -- axits settings
    ax <- list(
      zeroline = FALSE
      # gridcolor = #bdc3c7,
      # gridwidth = 1
    )
    sc = sc_seurat_umap_inter_base()
    plot_ly(sc,type="scatter", mode = "markers",
            x = ~UMAP_1, y = ~UMAP_2,
            color = ~ident,
            colors = getPalette(length(levels(sc$ident))),
            hoverinfo = "all",
            hovertext = paste0(sc$keys,"\n","nCount:",sc$nCount_RNA,"\n","nFeature:",sc$nFeature_RNA)) %>% layout(dragmode = "lasso", xaxis = ax, yaxis = ax)

  })

  output$seurat_umap_plot_inter <- renderPlotly(umap_inter())

  #### summary page ####

  dl_seurat_ob_base <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return(NULL)

    acv = input$dl_analysis
    # -- end acv part
    isolate({
      sc = sc_seurat_cluster()
      if("t-SNE" %in% input$dl_analysis){
        sc = sc_seurat_tsne()
        if("UMAP" %in% input$dl_analysis){
          withProgress(message="Running UMAP ...",{
            # -- default
            sc = RunUMAP(sc,
                         dims = 1:20,
                         learning.rate = input$umap_learning_rate,
                         spread = input$umap_spread,
                         min.dist = input$umap_min_dist,
                         verbose = FALSE)
          })
        }
      }
      else if("UMAP" %in% input$dl_analysis){
        sc = sc_seurat_umap()
      }
      return(sc)
    })
  })

  output$dl_seurat_ob <- downloadHandler(

    filename = "cells_object.rds",
    content = function(file) {
        saveRDS(object = dl_seurat_ob_base(),file)
      }
    )

  analysis_info <- reactive({
    if(input$load_sample_button == 0 || is.null(sc_seurat_cluster())) return("NULL")
    return(paste0("default info","__",input$dl_analysis,"__"))
  })

  output$analysis_info_text <- renderText(analysis_info())

})
