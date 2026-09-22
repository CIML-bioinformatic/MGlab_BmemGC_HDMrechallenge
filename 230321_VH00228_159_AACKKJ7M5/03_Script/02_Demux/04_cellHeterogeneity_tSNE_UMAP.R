# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# DIMENSIONAL REDUCTION (TSNE/UMAP)
###################################

# ..........................................................................................................
## @knitr heterogeneity_dimReduc
# ..........................................................................................................

nbPC_dimreduc=DIMREDUC_USE_PCA_NBDIMS
if(DIMREDUC_USE_PCA_NBDIMS>nbPC)
{
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'dimreduc' (", DIMREDUC_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_dimreduc = nbPC
}

sc10x = RunUMAP( sc10x, dims = 1:nbPC_dimreduc);
sc10x = RunTSNE( sc10x, dims = 1:nbPC_dimreduc);

# Save resulting coordinates for all cells as 'tsv' files
write.table( Embeddings(sc10x, reduction = "umap"), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_umap.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

write.table( Embeddings(sc10x, reduction = "tsne"), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_tsne.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");




# ..........................................................................................................
## @knitr dimreduc_ggplotly_colorBatch
# ..........................................................................................................

# Interactive plot using plotly
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)


# Gather data to be visualized together (cell name + numID + metrics)
cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                   sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = if(length( mito.genes)) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                   "percent.ribo" = if(length( ribo.genes)) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                   "Batch" = sc10x[["orig.ident", drop = TRUE]],
                   "HTO" = sc10x[["HTO_classification", drop = TRUE]]);

# Create text to show under cursor for each cell
hoverText = do.call(paste, c(Map( paste,
                                  c( "",
                                     "Cell ID: ",
                                     "# UMIs: ",
                                     "# Genes: ",
                                     if(length( mito.genes)) "% Mito: ",
                                     if(length( ribo.genes)) "% Ribo: "),
                                  cellsData[-ncol( cellsData)],
                                  sep = ""),
                             sep = "\n"));

## Do it in pure plotly
dimReducData = cbind(as.data.frame(Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"))),
                     Cluster = Idents( sc10x),
                     Batch   = sc10x[["orig.ident", drop = TRUE]],
                     HTO     = sc10x[["HTO_classification", drop = TRUE]]);

# Reusing 'hoverText' summarizing cell info (generated for split violin plots)
plotDimReduc = plot_ly(dimReducData,
                       x = as.formula( paste( "~", colnames( dimReducData)[1])),
                       y = as.formula( paste( "~", colnames( dimReducData)[2])),
                       #                       name = ~paste("Cluster", Cluster),
                       color = ~Batch,
                       height = 800) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = 5,
                            opacity = 0.6),
             text = hoverText,
             hoverinfo = "text+name") %>%
  #  hide_colorbar() %>%
  #  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list(format='svg'),
          modeBarButtons = list(list('toImage'),
                                list('zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d')));

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
     div(plotDimReduc, style = paste("flex : 0 0 800px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;")))
#)




# ..........................................................................................................
## @knitr dimreduc_ggplotly_colorHTOs
# ..........................................................................................................

# Interactive plot using plotly
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)


# Gather data to be visualized together (cell name + numID + metrics)
cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                   sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = if(length( mito.genes)) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                   "percent.ribo" = if(length( ribo.genes)) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                   "Batch" = sc10x[["orig.ident", drop = TRUE]],
                   "HTO" = sc10x[["HTO_classification", drop = TRUE]]);

# Create text to show under cursor for each cell
hoverText = do.call(paste, c(Map( paste,
                                  c( "",
                                     "Cell ID: ",
                                     "# UMIs: ",
                                     "# Genes: ",
                                     if(length( mito.genes)) "% Mito: ",
                                     if(length( ribo.genes)) "% Ribo: "),
                                  cellsData[-ncol( cellsData)],
                                  sep = ""),
                             sep = "\n"));

## Do it in pure plotly
dimReducData = cbind(as.data.frame(Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"))),
                     Cluster = Idents( sc10x),
                     Batch   = sc10x[["orig.ident", drop = TRUE]],
                     HTO     = sc10x[["HTO_classification", drop = TRUE]]);

# Reusing 'hoverText' summarizing cell info (generated for split violin plots)
plotDimReduc = plot_ly(dimReducData,
                       x = as.formula( paste( "~", colnames( dimReducData)[1])),
                       y = as.formula( paste( "~", colnames( dimReducData)[2])),
                       #                       name = ~paste("Cluster", Cluster),
                       color = ~HTO,
                       height = 800) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = 5,
                            opacity = 0.6),
             text = hoverText,
             hoverinfo = "text+name") %>%
  #  hide_colorbar() %>%
  #  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list(format='svg'),
          modeBarButtons = list(list('toImage'),
                                list('zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d')));

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
     div(plotDimReduc, style = paste("flex : 0 0 800px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;")))
#)




# ..........................................................................................................
## @knitr compareBatchHTO_table
# ..........................................................................................................

kable(table(sc10x[[c("HTOMICE_classification", "HTOTISSUE_classification")]]))




# ..........................................................................................................
## @knitr dimreduc_ggplot_filtering_selection
# ..........................................................................................................

# if QC EXPLORATION MODE is active, plot if the cells are outliers or not
print( DimPlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                group.by = "outlier") + 
         ggtitle( "Map of cells by \nQC filtering status"))



# ..........................................................................................................
## @knitr dimreduc_ggplot_filtering_variables
# ..........................................................................................................

# Plot the cells percentage of ribosomal genes in UMAP
FeaturePlot( sc10x, features = "percent.ribo") +
              scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
              ggtitle( "Map of cells with level of\n percentage of ribosomal genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by %of ribosomal genes or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                                  group.by = "outlier.percent.ribo") + 
                 ggtitle( "Map of cells by filter status\n on %ribosomal gene")
      )
}

# Plot the cells percentage of mitochondrial genes in UMAP  
FeaturePlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "percent.mito") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with level of percentage of\n mitochondrial genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by %of mitochondrial genes or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                  group.by = "outlier.percent.mito") + 
           ggtitle( "Map of cells by filter status\n on %mitochondrial gene")
  )
}

# Plot the cells RNA counts in UMAP
FeaturePlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nCount_RNA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with RNA counts")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by number of UMI count or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                  group.by = "outlier.nCount_RNA") + 
           ggtitle( "Map of cells by filter status \non nCount_RNA value")
  )
}

# Plot the cells genes counts in UMAP
FeaturePlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nFeature_RNA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with \nnumber of detected genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by number of feature count or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                  group.by = "outlier.nFeature_RNA") + 
           ggtitle( "Map of cells by filter status \non nFeature_RNA value")
  )
}



