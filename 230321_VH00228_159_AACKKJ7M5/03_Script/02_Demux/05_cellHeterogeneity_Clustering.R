# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# IDENTIFY CLUSTERS
###################

# ..........................................................................................................
## @knitr heterogeneity_identifyClusters
# ..........................................................................................................

# Identify clusters of cells by graph approach
nbPC_findclusters=FINDCLUSTERS_USE_PCA_NBDIMS
if(FINDCLUSTERS_USE_PCA_NBDIMS>nbPC)
{
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'findclusters' (", FINDCLUSTERS_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_findclusters = nbPC
}  
sc10x <- FindNeighbors(object    = sc10x,
                       k.param   = FINDNEIGHBORS_K, 
                       reduction = "pca",
                       dims      = 1:nbPC_findclusters,
                       verbose   = .VERBOSE);

sc10x <- FindClusters(object             = sc10x,
                      resolution         = FINDCLUSTERS_RESOLUTION,
                      algorithm          = FINDCLUSTERS_ALGORITHM,
                      temp.file.location = "/tmp/",
                      verbose            = .VERBOSE);


# Rename the clusters if required
if( exists( "RENAME_CLUSTERS") && SAMPLE_NAME %in% names( RENAME_CLUSTERS)){
  Idents( sc10x) = "seurat_clusters"
  sc10x = AddMetaData( sc10x, metadata = Idents( sc10x), col.name = "original_seurat_clusters")
  new_names = Idents( sc10x)
  for( cluster_id in levels( Idents( sc10x))){
    new_names[ which( Idents( sc10x) == cluster_id)] = RENAME_CLUSTERS[[ SAMPLE_NAME]][[ cluster_id]]
  }
  sc10x = AddMetaData( sc10x, metadata = new_names, col.name = "seurat_clusters")
  
}

Idents( sc10x) = "seurat_clusters"

# Show number of cells in each cluster
clustersCount = as.data.frame( table( Cluster = sc10x[[ "seurat_clusters" ]]), responseName = "CellCount");

# Define a set of colors for clusters (based on ggplot default)
clustersColor = hue_pal()( nlevels( Idents( sc10x)));
names( clustersColor) = levels( Idents( sc10x));


# Save cells cluster identity as determined with 'FindClusters'
write.table( data.frame(sc10x[["numID"]], identity = Idents(sc10x)), 
             file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsClusterIdentity.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# Also save cluster color attribution for reference
# Save cells cluster identity as determined with 'FindClusters'
write.table( clustersColor, 
             file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "clustersColor.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = FALSE,
             sep="\t");


# Create datatable
datatable( clustersCount,
           class = "compact",
           rownames = FALSE,
           colnames = c("Cluster", "Nb. Cells"),
           options = list(dom = "<'row'rt>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          columnDefs = list( # Center all columns
                                            list( targets = 0:(ncol(clustersCount)-1),
                                            className = 'dt-center')),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          paging = FALSE, # Disable pagination (show all)
                          processing = TRUE,
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE)) %>%
  # Add color from cluster
  formatStyle( columns = "Cluster",
               color = styleEqual( names(clustersColor), clustersColor),
               fontWeight = 'bold')


# Compute a matrix of average expression value by cluster (for each gene)
geneExpByCluster = do.call( rbind, 
                            apply( as.matrix( GetAssayData( sc10x)), # expression values
                                   1,                                # by rows
                                   tapply,                           # apply by group
                                   INDEX = Idents( sc10x),           # clusters IDs
                                   mean,                             # summary function
                                   simplify = FALSE));               # do not auto coerce

# Save it as 'tsv' file
write.table( geneExpByCluster, 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "normExpressionByCluster.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");



# ..........................................................................................................
## @knitr heterogeneity_identifyClusters_splitStats
# ..........................................................................................................

# Show UMIs Genes Mitochondrial and Ribosomal content split by cluster

# Gather data to be visualized together (cell name + numID + metrics)
Idents( sc10x) = "seurat_clusters"
cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                   sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = if(length( mito.genes)) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                   "percent.ribo" = if(length( ribo.genes)) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                   "Batch" = sc10x[["orig.ident"]],
                   "Cluster" = Idents( sc10x));

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

# Define size for panels (or assembled figure when using subplot)
panelWidth = 90 * nlevels(cellsData[["Cluster"]]);
panelHeight = 800;

# Generate plotly violin/jitter panels for #umis, #genes, and %mitochondrial stats
lypanel_umis  = plotViolinJitter( cellsData,
                                  xAxisFormula = ~as.numeric( Cluster),
                                  yAxisFormula = ~nCount_RNA,
                                  colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                  yAxisTitle = "# UMIs",
                                  hoverText = hoverText,
                                  traceName = ~paste( "Cluster", Cluster),
                                  xTicklabels = levels(cellsData[["Cluster"]]),
                                  panelWidth = panelWidth,
                                  panelHeight = panelHeight);

lypanel_genes = plotViolinJitter( cellsData,
                                  xAxisFormula = ~as.numeric( Cluster),
                                  yAxisFormula = ~nFeature_RNA,
                                  colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                  yAxisTitle = "# Genes",
                                  hoverText = hoverText,
                                  traceName = ~paste( "Cluster", Cluster),
                                  xTicklabels = levels(cellsData[["Cluster"]]),
                                  panelWidth = panelWidth,
                                  panelHeight = panelHeight);

lypanel_mitos = if(length( mito.genes)) plotViolinJitter( cellsData,
                                                          xAxisFormula = ~as.numeric( Cluster),
                                                          yAxisFormula = ~percent.mito,
                                                          colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                          yAxisTitle = "% Mito",
                                                          hoverText = hoverText,
                                                          traceName = ~paste( "Cluster", Cluster),
                                                          xTicklabels = levels(cellsData[["Cluster"]]),
                                                          panelWidth = panelWidth,
                                                          panelHeight = panelHeight) else NULL;

lypanel_ribos = if(length( ribo.genes)) plotViolinJitter( cellsData,
                                                          xAxisFormula = ~as.numeric( Cluster),
                                                          yAxisFormula = ~percent.ribo,
                                                          colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                          yAxisTitle = "% Ribo",
                                                          hoverText = hoverText,
                                                          traceName = ~paste( "Cluster", Cluster),
                                                          xTicklabels = levels(cellsData[["Cluster"]]),
                                                          panelWidth = panelWidth,
                                                          panelHeight = panelHeight) else NULL;

# Set panels as a list, define plotly config, and remove legend
panelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos, lypanel_ribos);
panelsList = lapply( panelsList, config, displaylogo = FALSE,
                     toImageButtonOptions = list( format='svg'),
                     modeBarButtons = list( list('toImage'),
                                            list( 'zoom2d', 'pan2d', 'resetScale2d')));

# Group plotly violin/jitter panels so we can synchronise axes and use highlight on the full set
plotPanels = layout( subplot( panelsList,
                              nrows = 4,
                              shareX = TRUE,
                              titleY = TRUE),
                     xaxis = list(title = "Seurat Cluster",
                                  showgrid = TRUE,
                                  tickvals = seq(nlevels(cellsData[["Cluster"]])),
                                  ticktext = levels(cellsData[["Cluster"]])),
                     showlegend = FALSE, # Remove eventual legends (does not mix well with subplot and highlight)
                     autosize = TRUE);

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
     div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
          div(plotPanels, style = paste("flex : 0 0 auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
#)







## @knitr heterogeneity_dimReduc_interactivePlot
# Interactive plot using plotly
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# # Using Seurat::Dimplot does not allow to control tooltip properly (neither using HoverLocator on ggplot, or plotly tooltip + added aestetic in ggplot). Mostly because 'label' and 'combine' change the type of returned object)
# # Compute tSNE dimensional reductionsand plot projection
# plotDimReduc = ggplotly( suppressWarnings( DimPlot( sc10x,
#                                                     reduction = ,
#                                                     label=TRUE,
#                                                     label.size = 6)) +
#                                              theme( legend.position = "none",
#                                                     plot.margin = margin(0, 0, 0, 0, "cm")),
#                          height = 600) %>%
#   config( displaylogo = FALSE,
#           toImageButtonOptions = list(format='svg'),
#           modeBarButtons = list(list('toImage'),
#                                 list('zoom2d', 'pan2d', 'resetScale2d')));

## Do it in pure plotly
dimReducData = cbind(as.data.frame(Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"))),
                     Cluster = Idents( sc10x));

# Reusing 'hoverText' summarizing cell info (generated for split violin plots)
plotDimReduc = plot_ly(dimReducData,
                       x = as.formula( paste( "~", colnames( dimReducData)[1])),
                       y = as.formula( paste( "~", colnames( dimReducData)[2])),
                       #                       name = ~paste("Cluster", Cluster),
                       color = ~Cluster,
                       colors = clustersColor,
                       height = 800) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = 5,
                            opacity = 0.6),
             text = hoverText,
             hoverinfo = "text+name") %>%
  # Compute median of coordinates for each group and use it as new data for Cluster text annotations (type='scatter', mode='text')
  add_trace( data = as.data.frame( do.call( rbind, by( dimReducData, dimReducData[["Cluster"]], function(x){ return( data.frame( lapply( x[1:2], median), "Cluster" = x[1, "Cluster"]));}))),
             x = as.formula( paste( "~", colnames( dimReducData)[1])),
             y = as.formula( paste( "~", colnames( dimReducData)[2])),
             type = "scattergl",
             mode = "text",
             text = ~Cluster,
             textfont = list( color = '#000000', size = 20),
             hoverinfo = 'skip',
             showlegend = FALSE) %>%
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





## @knitr rmd_heterogeneity_dimReduc_plot_byHTO

#### Dimensionality reduction 2D (ggplot + plotly conversion)
dimReducWidth  = 800;
dimReducHeight = 800;
initialPointSize = 5;


# Combine Dimreduc data and cells informations (keep coordinates in first cols)
# + add batch and HTO information for facetting
dimReducData = cbind(as.data.frame(Embeddings(sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"))),
                     sc10x[[c("orig.ident", "HTO_classification", "seurat_clusters")]]);


# Plot a thumbnail highlighting cluster cells
ggFigure = ggplot(  dimReducData,
                    aes_string( x = colnames( dimReducData)[1],
                                y = colnames( dimReducData)[2],
                                color = "seurat_clusters")) +
  geom_point( stroke = 0) + # Just ignore point size and alpha here as it will be overriden after plotly conversion
  facet_wrap(~HTO_classification)

ggFigure = ggFigure +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    #axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    #legend.position = "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));

plotDimReducGG =  ggplotly( ggFigure, 
                            width = dimReducWidth, 
                            height = dimReducHeight) %>% 
  style( type = "scattergl",   # trick to use webGL for rendering the scatterplot since toWebGL does not have width/height arguments, and use of width/height in layout is deprecated
         marker.size = initialPointSize, # Set point size and opacity with plotly so we can synchronize them with slider (ggplot and plotly units don't match)
         marker.opacity = 0.3) %>% 
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

cat("\n")
plotDimReducGG
cat("\n")



## @knitr compareClustersHTO_table

kable(table(sc10x[[c("seurat_clusters", "HTO_classification")]]))



## @knitr heterogeneity_dimReduc_thumbnail
# Non-interactive with large labels for generating thumbnails when analyzing monitored and marker genes expression
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)
# If contains several values, loop on them so we don't have to call chunk separately (which creates a new <p>aragraph in html)

# Replot the projection with colored clusters and large labels (add title if several plotted in the loop)
reductionVector = if(exists("useReduction")) useReduction else "umap";
for(currentReduction in reductionVector)
{
  ggFigure = suppressWarnings( DimPlot( sc10x, reduction = currentReduction, label=TRUE, label.size = 10)) +
    scale_color_manual( name = "Cluster", values = clustersColor) +  # Fix legend title + set colors
    theme( axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           legend.position = "none",
           plot.margin = margin( 0, 0, 0, 0, "cm"),
           plot.title = element_text( face = "bold",
                                      size = rel( 16/14),
                                      hjust = 0.5,
                                      vjust = 1,
                                      margin = margin( b = 7)));
  
  if( length( reductionVector) > 1) ggFigure = ggFigure + ggtitle( label = currentReduction);
  print( ggFigure);
}





## @knitr heterogeneity_dimReduc_with_clusters

# Plot the map with clusters with ggplot
print(
DimPlot( sc10x, 
         reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
         group.by = "seurat_clusters",
         label = TRUE) + 
  ggtitle( "Map of cells with clusters")
)

# Identify and plot the cells by their defined cluster gourp (cell type)
if( exists( "CLUSTER_GROUP_LIST")){
  
  cluster_group_annotation = rep( NA, length( Cells( sc10x)))
  names( cluster_group_annotation) = Cells( sc10x)
  
  for( cluster_group_id in names( CLUSTER_GROUP_LIST)){
    
    # Get the cells in the cluster group
    Idents( sc10x) = "seurat_clusters"
    cells_in_group <- WhichCells( sc10x, idents = CLUSTER_GROUP_LIST[[ cluster_group_id]])
    cluster_group_annotation[ cells_in_group] = cluster_group_id
    
    # Highlight the cells in UMAP
    print(
      DimPlot( sc10x, 
               label=TRUE, label.size = 10,
               cells.highlight= list(cells_in_group), 
               cols.highlight = c( SAMPLE_COLOR[ SAMPLE_NAME]), 
               cols= "grey") +
        ggtitle( paste( 'UMAP of Cells in', cluster_group_id)) +
        theme( axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = "none",
               plot.margin = margin( 0, 0, 0, 0, "cm"),
               plot.title = element_text( face = "bold",
                                          size = rel( 16/14),
                                          hjust = 0.5,
                                          vjust = 1,
                                          margin = margin( b = 7)))
    )
  }
  
  sc10x = AddMetaData( sc10x, metadata = cluster_group_annotation, col.name = "cluster_group")
  
  print(
    DimPlot( sc10x, 
             label=TRUE, label.size = 10, group.by = "cluster_group") +
             ggtitle( paste( "UMAP of Cells by cluster group (cell type)")) +
             theme( axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   legend.position = "none",
                   plot.margin = margin( 0, 0, 0, 0, "cm"),
                   plot.title = element_text( face = "bold",
                                              size = rel( 16/14),
                                              hjust = 0.5,
                                              vjust = 1,
                                              margin = margin( b = 7)))
  )
  
}



