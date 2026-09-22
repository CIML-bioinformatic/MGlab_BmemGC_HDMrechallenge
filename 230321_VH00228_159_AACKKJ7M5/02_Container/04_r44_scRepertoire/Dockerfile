# Based on rocker - https://github.com/rocker-org/rocker-versioned
FROM rocker/tidyverse:4.4


# ##############################
# R packages system dependencies
# ##############################

RUN apt-get update && apt-get install -y optipng \
                                         pandoc \
                                         libpng-dev \
                                         zlib1g-dev \
                                         libssl-dev \
                                         libcurl4-openssl-dev


# Install pip3 and set config to ignore error linked to not using virtual environment
RUN apt update && apt install -y python3-pip
RUN pip3 config set global.break-system-packages true


# Add Leiden Algorithm for clustering
RUN pip3 install leidenalg

# ##########
# R PACKAGES 
# ##########

#### BiocManager
RUN Rscript -e 'install.packages( "BiocManager")'


#### RMarkdown & Reporting
RUN Rscript -e 'install.packages( "knitr")'
RUN Rscript -e 'install.packages( "rmarkdown")'
RUN Rscript -e 'install.packages( "bookdown")'
RUN Rscript -e 'install.packages( "htmltools")'      # browsable


#### Figures & layout management
# ggplot2
RUN Rscript -e 'install.packages( "ggplot2")'
RUN Rscript -e 'install.packages( "cowplot")'        # plot_grid, themes, ...
RUN Rscript -e 'install.packages( "ggpubr")'         # add_summary, geom_signif, ...
RUN Rscript -e 'install.packages( "ggrepel")'        # geom_text_repel, geom_label_repel
RUN Rscript -e 'install.packages( "gridExtra")'      # grid.arrange, ...
RUN Rscript -e 'install.packages( "hexbin")'         # stat_binhex
RUN Rscript -e 'BiocManager::install( "patchwork")'  # +/ operators for ggplots
RUN Rscript -e 'BiocManager::install("parcats")'

# plotly
RUN Rscript -e 'install.packages( "plotly")'

# general
RUN Rscript -e 'install.packages( "gplots")'         # heatmap.2
RUN Rscript -e 'install.packages( "heatmaply")'      # heatmaply (interactive)
RUN Rscript -e 'BiocManager::install( "iheatmapr")'  # iheatmap (interactive, uses plotly), dependencies OK with BiocManager
RUN Rscript -e 'install.packages( "pheatmap")'       # pheatmap
RUN apt-get install -y libxt-dev libproj-dev
RUN Rscript -e 'BiocManager::install( "EnhancedVolcano")'
RUN Rscript -e 'BiocManager::install( "ComplexHeatmap")'

#### Reporting
RUN Rscript -e 'install.packages( "DT")'             # datatable
RUN Rscript -e 'install.packages( "htmlwidgets")'    # JS (for datatable)
RUN Rscript -e 'install.packages( "kableExtra")'     # kable_styling
RUN Rscript -e 'install.packages( "pander")'         # pander


#### General
RUN Rscript -e 'install.packages( "digest")'         # digest (hashing)
RUN Rscript -e 'install.packages( "dplyr")'
RUN Rscript -e 'install.packages( "forcats")'        # fct_inorder (factors)
RUN Rscript -e 'install.packages( "funr")'           # get_script_path
RUN Rscript -e 'install.packages( "future")'         # plan (multicore)
RUN Rscript -e 'install.packages( "fs")'             # path_sanitize
RUN Rscript -e 'install.packages( "magrittr")'
RUN Rscript -e 'install.packages( "reshape")'        # melt
RUN Rscript -e 'install.packages( "scales")'         # hue_pal


#### Technology specific
#RUN Rscript -e 'BiocManager::install( "Seurat")'
RUN Rscript -e 'install.packages( "Matrix")'         # Force install Matrix 1.5.3 instead of 1.5.1 that cause a bug in FindNeighbors function


#### Functional analysis
RUN Rscript -e 'BiocManager::install( "org.Mm.eg.db" )' # Mouse annotation reference
RUN Rscript -e 'BiocManager::install( "org.Hs.eg.db" )' # Human annotation reference

RUN apt-get update && apt-get install -y libxml2-dev libigraph*
RUN Rscript -e 'BiocManager::install( "clusterProfiler")'
RUN Rscript -e 'BiocManager::install( "rrvgo")'

#### Custom dependencies

RUN apt-get update && apt-get install -y curl \
                                         libbz2-dev \
                                         libjpeg-dev \
                                         liblzma-dev \
                                         libxml2-dev

RUN Rscript -e 'BiocManager::install( "GenomicRanges")'
RUN Rscript -e 'BiocManager::install( "rtracklayer")'
RUN Rscript -e 'BiocManager::install( "topGO")'
RUN Rscript -e 'BiocManager::install( "Gviz")'

# Additional manual investigation tools
RUN Rscript -e 'BiocManager::install( "pcaExplorer")' # Also depends on libigraph installed earlier for clusterProfiler  
RUN Rscript -e 'devtools::install_github("JSB-UCLA/scDEED")'
RUN Rscript -e 'BiocManager::install( "babelgene")'

RUN apt install -y libhdf5-dev
RUN apt install -y libgsl-dev


RUN Rscript -e 'devtools::install_github("satijalab/seurat", "seurat5")'
RUN Rscript -e 'devtools::install_github("satijalab/seurat-data",ref= "seurat5")'
RUN Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'
RUN Rscript -e 'BiocManager::install( "Signac")'
RUN Rscript -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")'



# Additional graphic tools for reporting
RUN Rscript -e 'BiocManager::install( c("pathview", "ggh4x", "ggnewscale"))'

# Install scReprtoire which requires R>=4.3
RUN Rscript -e 'remotes::install_github(c("BorchLab/immApex", "BorchLab/scRepertoire"))'


