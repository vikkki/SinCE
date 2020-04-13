![](https://img.shields.io/github/license/vikkki/SinCE.svg)
[![Build
Status](https://travis-ci.org/vikkki/SinCE.svg?branch=master)](https://travis-ci.org/vikkki/SinCE)
[![Coverage
status](https://codecov.io/gh/vikkki/SinCE/branch/master/graph/badge.svg)](https://codecov.io/github/vikkki/SinCE?branch=master)

# SinCE

SinCE (single cell cluster explorer) is shinyapp for the analysis and visualization of single cell expression data. Based on R packages including [Seurat](https://github.com/satijalab/seurat), and [plotly](https://plotly.com/r/), SinCE working on providing a graph user interface for researchers who want to investigate thier single cell expression data derectly, and providing plots and data for further anaylisis and publication.

To interact, you can adjust parameters of analysis including PCA, t-SNE, and UMAP, as well as compare and costumize clusters. For further investment, selected cell barcodes can be downloaded.

While developing SinCE, I'm more than glad to have your feedback;)
Hope we can have a good time play with data.

### Installation and lauch

```R
# install.packages("devtools")
devtools::install_github("vikkki/SinCE")

# lauch
library(SinCE)
SinCE()

```

### Input data type
Now there are three ways to load data to our app:
1. Demo data;
2. Matrix data in csv/tsv;
3. Matrix output of 10X cellranger count.

##### For .csv file,
A typical expression matrix could look like these:
```
"","AAACATACAACCAC.1","AAACATTGAGCTAC.1","AAACATTGATCAGC.1","AAACCGTGCTTCCG.1","AAACCGTGTATGCG.1","AAACGCACTGGTAC.1"
"MIR1302.10",0,0,0,0,0,0
"FAM138A",0,0,0,0,0,0
"RP4.669L17.2",0,0,0,0,0,0
"TNFRSF18",0,0,0,0,0,0
```
Seperator could change, and when you upload CSV file, it's a good chioce ot make sure that you have the right seperator selected.

##### For 10X output file,
The matrices output may looks like:
```shell
$ cd /home/jdoe/runs/sample345/outs
$ tree filtered_feature_bc_matrix
filtered_feature_bc_matrix
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
0 directories, 3 files
```
These three file in matrix output folder could be upload indivadually into the app. For the further information of 10X output, check out [here on 10X's website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices).

---
### Workflow
Once uploaded the data, by click on different tabs, different functions would be excuted. You can change and adjust your analysis by adjusting costumized input.

<img src = "https://raw.githubusercontent.com/vikkki/SinCE/master/inst/SinCEapp/www/workflow.png" height ="500" align = "center" />

---
### Acknowledgement
This on-going shinyapp is part of the master's capstone project has a great support by [the Department of Translational Genomics](https://dtg.usc.edu/site/index.php/bioinformatics/). 
