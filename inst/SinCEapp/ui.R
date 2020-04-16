
library(shiny, verbose=FALSE)
library(shinythemes, verbose=FALSE)
library(plotly)
#library(heatmaply, verbose=FALSE) # for makeing interactive heatmap

# Define UI for application that draws a histogram
shinyUI(
  navbarPage("SinCE pipe",
             theme = shinytheme("simplex"),

             #### upload data ####
             tabPanel("Load data",
                      sidebarPanel(

                        ###### upload demo ######
                        selectInput("demo_select", label = h4("Select a demo dataset:"),
                                    choices = list(
                                                   #"2.5k pbmc, MGC" = "demo_2.5k_pbmc",
                                                   "1k pbmc, 10x chemistry v3" = "demo_1k_pbmc"
                                                   #"5k pbmc, 10x chemistry v3, cellranger 3" = "demo_5k_pbmc",
                                                   #"10k pbmc, 10x chemistry v3" = "demo_10k_pbmc"
                                                   #,"mini sample(download test, not for analyysis)" = "mini_sample"
                                                   ),
                                    selected = "demo_1k_pbmc"),



                        p(downloadLink("download_demo_data",label = "Download the selected demo data")),
                        hr(),

                        ###### upload csv/tsv #######
                        h4("Or, upload your expression matrix:"),
                        radioButtons("csv_seperator", "select the upload file seperator:",
                                     choices = c("comma" = ",",
                                                 "tab" = "\t",
                                                 "semicolon" = ";"),
                                     selected = ","
                        ),
                        fileInput('csv_matrix_file',
                                               h5('Upload a barcode ~ expression matrix file(CSV or TSV)'),
                                               accept = c('text/csv',
                                                          'text/comma-separated-values',
                                                          '.csv')),
                        hr(),

                        ###### upload matrix ######
                        h4("Or, upload the output files of cellranger count:"),

                        fileInput('barcods_file', "barcodes.tsv.gz:", accept = c('.tsv.gz','.tsv')),
                        fileInput('features_file', "features.tsv.gz:", accept = c('.tsv.gz','.tsv')),
                        fileInput('matrix_file', "matrix.mtx.gz", accept = c('.mtx.gz','.mtx')),

                        br(),
                        textOutput("load_alart_ui"),
                        div(actionButton("load_sample_button", "Load data"),align = "right"),
                        p(HTML("<div align=\"right\"> <A HREF=\"javascript:history.go(0)\">Reset</A></div>" )),
                        p(HTML("<div align=\"right\"> <A HREF=\"#\">Need Help?</A></div>" ))
                      ),
                      mainPanel(
                        h4(textOutput("seurat_load_message")),
                        h4("A cornor of loaded matrix:"),
                        tableOutput("load_summary"),
                        #verbatimTextOutput("path"),
                        br(),
                        h4(textOutput("create_mat_alart"), style = "color:red"),
                        h4(textOutput("create_seurat_alart"),style = "color:red"),

                        h5("Welcome! To start, first you need to load a matrix for analyse, we now provide three ways to load data to our app:"),
                        h5("1. Select a demo data;"),
                        h5("2. Upload matrix data in csv/tsv;"),
                        h5("3. Upload matrix output of 10x cellranger count."),
                        br(),

                        p(HTML("<h5 align=\"left\">(A sample be can found at <A HREF=\"https://github.com/vikkki/SinCE\">Github page</A>).</h5>" )),
                        h5("After you select a demo dataset or complete uploading your file(s), click the \"Load data\" button and it would take minutes to prepare your data."),

                        h5("Then click on each tab for different fuctions."),
                        p("Note: If you need to reload the data, please click the reset button first.",style = "color:#7D463B"),
                        br(),
                        h5(strong("WorkFlow:")),
                        HTML('<img src="workflow.png", height="500px" />'),


                        h4("Building progress:"),
                        p("QC --------  QC plot implemented, user input threshold finished;"),
                        p("Variable features --------  function completed;"),
                        p("PCA -------- completed, allowed visualization PCs:"),
                        p("t-SNE -------- function completed, allowed user input max.iter and number of PCs used"),
                        p("UMAP -------- function completed:"),
                        p("Heatmap -------- basic functions completed"),
                        p("Clusters -------- functions completed, cluster vs cluster analysis finished; user difined cluster finished"),
                        p("Expression -------- gene plots (violin, distribution, ridge, scattar, co-expression) finished"),
                        p("Help documents -------- not start yet(host on a seperated static website or github wiki)"),
                        p("Seurat object download -------- basic function completed"),
                        p("Parameter records and reproducibility -------- not start"),

                      )# end main panel(load data)
             ), # end load data tab

             #### data prep ####
             tabPanel("QC",
                      sidebarPanel(
                        h4("Input the criteria of cells to include in following analysis (intersection would be applied):"),
                        sliderInput("clean_feature", "Number of features range:",
                                    min = 1, max = 5000, value = par_templete[["clean_feature"]]), # -- max of this slider depends on sample (by observe())
                        sliderInput("clean_count", "RNA count range:",
                                    min = 1, max = 30000, value = par_templete[["clean_count"]]), # -- max of this slider depends on sample (by observe())
                        sliderInput("clean_mito", "Mitochondrial gene percentige range:",
                                    min = 0, max = 100, value = par_templete[["clean_mito"]]),
                        div(actionButton("clean_button", "Clean data"),align = "right"),
                        h4(textOutput("cleaned_cell_number")),
                        helpText("Note: if the clean data button isn't clicked, default filter will be apply."),
                        p(HTML("<div align=\"right\"> <A HREF=\"#\">?</A></div>" ))),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Sample QC",
                                   #checkboxInput("show_range", label = "Show the criteria range of cleaning on raw data violin plots", value = FALSE),
                                   p("Visualize QC metrics as a violin plot"),
                                   plotOutput("sample_qc_violin"),
                                   p("Scatter plots of cells, show RNA count vs the percentage of reads that map to the mitochondrial genome,
                                     and RNA count vs the number of features of each cell"),
                                   plotOutput("sample_qc_feature")
                                  ),
                          tabPanel("Cleaned data",
                                   p("Visualize QC metrics as a violin plot"),
                                   plotOutput("cleaned_qc_violin"),
                                   p("Scatter plots of cells, show RNA count vs the percentage of reads that map to the mitochondrial genome,
                                     and RNA count vs the number of features of each cell"),
                                   plotOutput("cleaned_qc_feature")
                                   ),
                          tabPanel("Cluster QC",
                                   p("Visualize QC metrics as a violin plot based on clusters, data scaled"),
                                   plotOutput("cluster_qc_violin")
                                   )
                          )
                        ) # end main panel
                      ), # end QC tab

             #### variation ####
             tabPanel("Variation",
                      sidebarPanel(
                        p("calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). "),
                        numericInput('topn', 'Display most highly variable genes:', 10, min = 1, max = 150),
                        p(HTML("<div align=\"right\"> <A HREF=\"#\">?</A></div>" ))
                      ),
                      mainPanel(
                        p("Most highly variable genes:"),
                        verbatimTextOutput("seurat_topn_variable_features"),
                        plotOutput("seurat_variable_feature_plot"),
                        fluidRow(
                          column(6,plotOutput("seurat_variable_count_plot")),
                          column(6,plotOutput("seurat_variable_density_plot")))
                      )
             ),


             #### clusters and genes ####
             tabPanel("Cluster",
                      sidebarPanel(
                        h4("Clustering options:"),
                        sliderInput("seurat_cluster_pc", "Number of PCs used in finding neighbors:",
                                    min = 2, max = n_pca,
                                    value = par_templete[["seurat_cluster_pc"]]),

                        numericInput("seurat_cluster_resolution", label = "Resolution of finding clusters:",
                                    value = par_templete[["seurat_cluster_resolution"]]),
                        helpText("For datasets of around 3K cells, reselution of 0.4-1.2 typically returns good results."),

                        hr(),

                        numericInput("topn_cluster", label = "Show top n markers of each cluster:", value = 3),
                        p(HTML("<div align=\"right\"> <A HREF=\"#\">?</A></div>" ))
                      ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Define clusters",
                                   plotOutput("seuret_tsne_plot_2", height = "400px"),
                                   p("Raw cluster names:"),
                                   verbatimTextOutput("cluster_names_raw"),
                                   p("Current cluster names:"),
                                   verbatimTextOutput("cluster_names"),
                                   hr(),
                                   textInput("input_cluster_names", label = "Input cluster names:",width = "100%"),
                                   p("Note: type in new names coresponding to current clusters and seperate them by ',' . Cluster would combine if same new name are assigned to them.
                                     And if you need to re-cluster the sample, adjust the cluster option sliderbars in side panel."),
                                   verbatimTextOutput("input_cluster_names_input"),
                                   h4(textOutput("name_message", inline = TRUE), style = "color:red"),
                                   div(actionButton("name_button", "Rename clusters"),align = "right", inline = TRUE)
                                   ),
                          # tabPanel("Cluster makers",
                          #          DT::dataTableOutput("seurat_cluster_marker_table")
                          # ),
                          tabPanel("Clusters vs Clusters",

                                   br(),
                                   fluidRow(
                                     column(6,wellPanel(
                                       checkboxGroupInput("cluster_group1",
                                                          label = h4("Select cluster(s) as group 1:"),
                                                          choices = c(0,1),
                                                          selected = 0)
                                     )),
                                     column(6,wellPanel(
                                       checkboxGroupInput("cluster_group2",
                                                          label = h4("Select cluster(s) as group 2:"),
                                                          choices = c(1),
                                                          selected = 1)
                                     )#group UI will be update by observer()
                                     ),
                                     div(actionButton("find_markers_button", "Find markers!"),align = "right", inline = TRUE)
                                   ),
                                   # end fluiDrow
                                   hr(),
                                   p("Markers between two groups:"),
                                   DT::dataTableOutput("seurat_cross_cluster_markers_show"),
                                   div(downloadButton("dl_cluster_markers", "Download marker table"), align = "right")
                          ),
                          tabPanel("Custormize clusters",
                                   plotlyOutput("tsne_for_brush"),
                                   hr(),
                                   h4(textOutput("cluster_brush_info")),
                                   div(downloadButton("dl_select_cells", "Download selected cell barcodes"), align = "right"),
                                   hr(),
                                   textInput("selected_cell_name","Cluster name for selected cells:",value = "MyCluster"),
                                   h4(textOutput("rename_selected_message", inline = TRUE), style = "color:red"),
                                   div(actionButton("rename_selected_cell_button", "Rename"),align = "left")

                          )
                        )
                      )
             ), # cluster tab end


             #### data PCA ####
             tabPanel("PCA",
                      sidebarPanel(
                        numericInput('pcx', 'PC on x-axis:', 1, min = 1, max = 150),
                        numericInput('pcy', 'PC on y-axis:', 2, min = 1, max = 150),
                        sliderInput('pca_point_size','Point size:', 1, min = 0.1, max = 2),
                        checkboxInput("pca_label",label = "Show cluster identities in plot.", value = TRUE),
                        hr(),
                        numericInput('pc_genes_n',"Number of genes:", 30, min = 1, max = 100),
                        sliderInput("pcn", "Visualize top genes associated and heatmap of principal components:",
                                    min = 1, max = 50,
                                    value = c(1,3)
                                    ),
                        numericInput('pca_heatmap_gene',"Number od cells in PC heatmap:",
                                     min = 1, max = 1000,
                                     value = 500), ################### should be a reactive ui, responding to input cell number
                        verbatimTextOutput("pcn"),
                        hr(),

                        p(HTML("<div align=\"right\"> <A HREF=\"#\">?</A></div>" ))),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("PCA plot",
                                   checkboxInput("pca_inter",label = "Show interactive scatter.", value = TRUE),
                                   conditionalPanel("input.pca_inter > 0",
                                                    plotlyOutput("seuret_pca_plot_inter")),
                                   conditionalPanel("input.pca_inter == 0",
                                                    plotOutput("seuret_pca_plot"))
                          ),
                          tabPanel("Genes on PCs",
                                   plotOutput("seurat_pca_genes",height="600px")),
                          tabPanel("PCs heatmap",
                                   plotOutput("seurat_pca_heatmap",height="600px")),
                          tabPanel("Dimensionality",
                                   p("Calculation of dimention scores could take more than 10 min. Click the button to start and have some coffee."),
                                   div(actionButton("pca_score_button", "Calculate scores"),align = "right", inline = TRUE),
                                   plotOutput("seurat_pca_score"))
                          )
                        )


             ),

             #### makers  ####
             tabPanel("Expression",
                      sidebarPanel(
                        h4("Ploting options:"),
                        textInput("plot_cluster_feature", label = "Input features to visualize:", value = par_templete[["plot_cluster_feature"]]),
                        helpText("Note: seperate multiple features by comma."),
                        #hr(),
                        #verbatimTextOutput("cluster_feature"),

                        hr(),
                        selectInput("featureplot_reduction", label = "Select a reduction method the feature expression plot based on:",
                                    choices = list(
                                      "UMAP" = "umap",
                                      "PCA" = "pca",
                                      "t-SNE" = "tsne"
                                    ),
                                    selected = "pca"),

                        p(HTML("<div align=\"right\"> <A HREF=\"#\">?</A></div>" ))
                      ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Violin plot",
                                   plotOutput("seurat_cluster_feature_vln")
                          ),

                          tabPanel("Features distribution",
                                   plotOutput("seurat_cluster_feature_plot", height = "800px"),
                                   hr(),
                                   textInput("feature_inter","View feature:", value = "nCount_RNA"),
                                   plotlyOutput("seurat_feature_inter")
                          ),

                          tabPanel("Ridge plot",
                                   p("Visualize single cell expression distributions in each cluster."),
                                   plotOutput("seurat_cluster_ridge")

                          ),
                          tabPanel("Dot plot",
                                   p("In the dot plots, the size of the dot represents the percentage of cells expressing the feature in each cluster. The color represents the average expression level"),
                                   plotOutput("seurat_cluster_dot")
                                   ),
                          tabPanel("Features correlation",
                                   textInput("featurex","Feature on x_axis:", value = par_templete[["cor_feature_x"]]),
                                   textInput("featurey","Feature on y_axis:", value = par_templete[["cor_feature_y"]]),
                                   hr(),
                                   p("Scatter plot of two features feature expression."),
                                   plotOutput("seurat_feature_scatter"),
                                   p("Visualize co-expression of two features."),
                                   plotOutput("seurat_feature_coex")

                          )

                        )# tabset end
                      )# main of marker panel end
             ), # markers tab end

             # #### Heatmap ####
             # tabPanel("Heatmap",
             #          sidebarPanel(
             #            numericInput("topn_heatmap", label = "Input a top n number:" , value = 10),
             #            helpText("Generate an expression heatmap with top n markers for each cluster"),
             #            br(),
             #            checkboxGroupInput("heatmap_clusters","Select clusters to include in heatmap",
             #                               choices = c(0,1),
             #                               selected = 0), # will be update by observer
             #            br(),
             #            radioButtons("heatmap_slot","Slot use in heapmap:",
             #                         choices = c("Raw counts" = "counts",
             #                                     "Normalized data" = "data",
             #                                     "Scaled.data" = "scale.data"),
             #                         selected = "scale.data")
             #          ),
             #          mainPanel(
             #            br(),
             #            plotOutput("seurat_marksers_heatmap",height = "500px")
             #          )
             # ),



             #### t-SNE ####
             tabPanel("t-SNE",
                      sidebarPanel(
                        radioButtons("seurat_tsne_run_method", "Select a t-SNE running method:",
                                     choices = c("Rtsne" ,
                                                 "FIt-SNE"),
                                     selected = par_templete[["seurat_tsne_run_method"]]
                        ),

                        p("Rtsne: Use the Rtsne package Barnes-Hut implementation of tSNE."),
                        p("FIt-SNE: Use the FFT-accelerated Interpolation-based t-SNE."),
                        hr(),

                        sliderInput("seurat_tsne_max_pc", "Number of PCs used in t-SNE:",
                                    min = 2, max = n_pca,
                                    value = par_templete[["seurat_tsne_pc"]]),

                        sliderInput("seurat_tsne_max_iter", "Max iteration to run t-SNE:",
                                    min = 100, max = 6000,
                                    value = par_templete[["seurat_tsne_max_iter"]]),
                        hr(),
                        sliderInput('tsne_point_size','Point size:', 1, min = 0.1, max = 2),
                        checkboxInput("tsne_label",label = "Show cluster identities in plot.", value = TRUE),
                        p(HTML("<div align=\"right\"> <A HREF=\"#\">?</A></div>" ))
                        ),
                      mainPanel(
                        checkboxInput("tsne_inter",label = "Show interactive scatter.", value = TRUE),
                        conditionalPanel("input.tsne_inter > 0",
                                         plotlyOutput("seuret_tsne_plot_inter", height = "500px")),
                        conditionalPanel("input.tsne_inter == 0",
                                         plotOutput("seuret_tsne_plot", height = "500px"))

                      )
             ),

             #### UMAP ####
             tabPanel("UMAP",
                      sidebarPanel(
                        sliderInput('umap_learning_rate',"Initial learning rate for the embedding optimization:",
                                    par_templete[["umap_learning_rate"]], min = 0.1, max = 5),
                        sliderInput('umap_min_dist',"min.dist:",
                                    par_templete[["umap_min_dist"]], min = 0.001, max = 0.5),
                        helpText("This controls how tightly the embedding is allowed compress points together.
                                 Larger values ensure embedded points are moreevenly distributed, while smaller
                                 values allow the algorithm to optimise more accurately with regard to local structure."),
                        sliderInput('umap_spread',"Spread:",
                                    par_templete[["umap_spread"]], min = 0.01, max = 5),
                        helpText("The effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are."),
                        hr(),
                        sliderInput('umap_point_size','Point size:', 1,  min = 0.1, max = 2),
                        checkboxInput("umap_label",label = "Show cluster identities in plot.", value = TRUE),
                        p(HTML("<div align=\"right\"> <A HREF=\"#\">?</A></div>" ))),
                      mainPanel(
                        checkboxInput("umap_inter",label = "Show interactive scatter.", value = TRUE),
                        conditionalPanel("input.umap_inter > 0",
                                         plotlyOutput("seurat_umap_plot_inter", height = "500px")),
                        conditionalPanel("input.umap_inter == 0",
                                         plotOutput("seurat_umap_plot", height = "500px"))
                      )
             ),

             tabPanel("Download",
                      sidebarPanel(
                        p("Here you can Download the Suerat object of current single cell expression data for further study:) The defualt file for download contains cleaned expression matrix itself and basic analysis info includes clustering and PCA. Aditional t-SNE and UMAP result are optional."),
                        br(),
                        checkboxGroupInput("dl_analysis","Select extra analysis to include:",
                                           choices = c("UMAP","t-SNE")),
                        br(),
                        helpText("There could be a few seconds before download box popup."),
                        div(downloadButton("dl_seurat_ob", "Download Seurat object"), align = "right")
                      ),
                      mainPanel(
                        h5("Analysis parameter summary(not work):."),
                        verbatimTextOutput("analysis_info_text")
                        #div(downloadButton("dl_sc_pars", "Download parameter JSON file"), align = "right")
                        #plotlyOutput("inter_test")
                        )
             )
             #### future functions ####

))
