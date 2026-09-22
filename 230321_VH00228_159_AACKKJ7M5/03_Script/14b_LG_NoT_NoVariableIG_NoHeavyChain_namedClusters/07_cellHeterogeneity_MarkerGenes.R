# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################

# MARKER GENES
##############

# ..............................................................................
## @knitr heterogeneity_markerGenes
# ..............................................................................

# Identify marker genes
markers = FindAllMarkers( object          = sc10x,
                          test.use        = FINDMARKERS_METHOD,
                          only.pos        = FINDMARKERS_ONLYPOS,
                          min.pct         = FINDMARKERS_MINPCT,
                          logfc.threshold = FINDMARKERS_LOGFC_THR,
                          return.thresh   = FINDMARKERS_PVAL_THR,
                          random.seed     = SEED,
                          verbose         = .VERBOSE);

# Save markers list as 'tsv' table
write.table( markers,
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "MarkerGenes.tsv")),
             quote = FALSE,
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers_table = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP_TABLE)) x else head( x, n = FINDMARKERS_SHOWTOP_TABLE));
});

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers_heatmap = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP_HEATMAP)) x else head( x, n = FINDMARKERS_SHOWTOP_HEATMAP));
});

# Merge marker genes in a single data.frame
topMarkers_tableDF = do.call( rbind, topMarkers_table);
topMarkers_heatmapDF = do.call( rbind, topMarkers_heatmap);

# Select and order columns to be shown in datatable
topMarkers_tableDT = topMarkers_tableDF[c("gene", "cluster", "avg_log2FC", "p_val_adj")]


## Same strategy but removing IG genes from markers (for heatmap)


markersIgVariable = grepl("^ig.*v.*", markers[['gene']], ignore.case = TRUE)
markersIgHeavy = grepl("Igh(d|m|a|g1|g2a|g2b|g3|e)$", markers[['gene']], ignore.case = TRUE)

markers_noIG = markers[ !(markersIgVariable | markersIgHeavy) ,]

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers_table_noIG = by( markers_noIG, markers_noIG[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP_TABLE)) x else head( x, n = FINDMARKERS_SHOWTOP_TABLE));
});

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers_heatmap_noIG = by( markers_noIG, markers_noIG[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP_HEATMAP)) x else head( x, n = FINDMARKERS_SHOWTOP_HEATMAP));
});


# Merge marker genes in a single data.frame
topMarkers_tableDF_noIG = do.call( rbind, topMarkers_table_noIG);
topMarkers_heatmapDF_noIG = do.call( rbind, topMarkers_heatmap_noIG);

# Select and order columns to be shown in datatable
topMarkers_tableDT_noIG = topMarkers_tableDF_noIG[c("gene", "cluster", "avg_log2FC", "p_val_adj")]


# ..............................................................................
## @knitr heterogeneity_markerGenes_table
# ..............................................................................

# Create datatable
datatable( topMarkers_tableDT,
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Cluster", "Avg. LogFC", "Adj. Pvalue"),
           caption = paste(ifelse( is.null( FINDMARKERS_SHOWTOP_TABLE), "All", paste("Top", FINDMARKERS_SHOWTOP_TABLE)), "marker genes for each cluster"),
           extensions = c('Buttons', 'Select'),
           options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          buttons = exportButtonsListDT,
                          columnDefs = list(
                            list( # Center all columns except first one
                              targets = 1:(ncol( topMarkers_tableDT)-1),
                              className = 'dt-center'),
                            list( # Set renderer function for 'float' type columns (LogFC)
                              targets = ncol( topMarkers_tableDT)-2,
                              render = htmlwidgets::JS("function ( data, type, row ) {return type === 'export' ? data : data.toFixed(4);}")),
                            list( # Set renderer function for 'scientific' type columns (PValue)
                              targets = ncol( topMarkers_tableDT)-1,
                              render = htmlwidgets::JS( "function ( data, type, row ) {return type === 'export' ? data : data.toExponential(4);}"))),
                          #fixedColumns = TRUE, # Does not work well with filter on this column
                          #fixedHeader = TRUE, # Does not work well with 'scrollX'
                          lengthMenu = list(c( 10, 50, 100, -1),
                                            c( 10, 50, 100, "All")),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          processing = TRUE,
                          #search.regex= TRUE, # Does not work well with 'search.smart'
                          search.smart = TRUE,
                          select = TRUE, # Enable ability to select rows
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE)) %>%
  # Add bar relative to logFC
  formatStyle( columns = "avg_log2FC",
               background = styleColorBar( data = range( topMarkers_tableDT[["avg_log2FC"]]), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  # Add color from cluster
  formatStyle( columns = "cluster",
               backgroundColor = styleEqual( names(clustersColor),
                                             scales::alpha(clustersColor, 0.3)));




# ..............................................................................
## @knitr heterogeneity_markerGenes_heatmap_prepare_noIG
# ..............................................................................

# Get the matrix of expression and associated clusters from Seurat object
expMat = as.matrix( GetAssayData( sc10x));
clusterID = Idents( sc10x);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkers_heatmapDF_noIG[["gene"]];
clusterOrdering = order( clusterID);

expMat = expMat[topMarkersGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Prefix marker gene name with cluster name to prevent eventual duplicated row names
rownames(expMat) = paste0( topMarkers_heatmapDF_noIG[["cluster"]], ".", rownames(expMat))

# Compute the mean of top markers in clusters and produce matrix with the result
meanExp_ByCluster = by( as.data.frame( t(expMat)), 
                        clusterID, #as.factor(clusterID), 
                        colMeans)
# Flatten to a numeric matrix
meanExpMat = do.call( cbind, meanExp_ByCluster)



# ..............................................................................
## @knitr heterogeneity_markerGenes_heatmap_mean_noIG
# ..............................................................................

cat("\n \n")
pheatmap( meanExpMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          scale = "row",
          annotation_row = data.frame(Markers = topMarkers_heatmapDF_noIG[ rownames( meanExpMat), "cluster"], stringsAsFactors = FALSE, row.names = rownames( meanExpMat) ),
          annotation_col = data.frame(Cluster = factor(colnames(meanExpMat), levels = levels(clusterID)), stringsAsFactors = FALSE, row.names = colnames(meanExpMat)),
          annotation_colors = list( Markers = clustersColor,
                                    Cluster = clustersColor),
          show_colnames = TRUE,
          fontsize_row = 5,
          main = "Scaled mean normalized expression\nby cluster of top markers");



# ..............................................................................
## @knitr heterogeneity_markerGenes_heatmap_prepare
# ..............................................................................

# Get the matrix of expression and associated clusters from Seurat object
expMat = as.matrix( GetAssayData( sc10x));
clusterID = Idents( sc10x);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkers_heatmapDF[["gene"]];
clusterOrdering = order( clusterID);

expMat = expMat[topMarkersGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Prefix marker gene name with cluster name to prevent eventual duplicated row names
rownames(expMat) = paste0( topMarkers_heatmapDF[["cluster"]], ".", rownames(expMat))

# Compute the mean of top markers in clusters and produce matrix with the result
meanExp_ByCluster = by( as.data.frame( t(expMat)), 
                        clusterID, #as.factor(clusterID), 
                        colMeans)
# Flatten to a numeric matrix
meanExpMat = do.call( cbind, meanExp_ByCluster)

# ..............................................................................
## @knitr heterogeneity_markerGenes_heatmap_mean
# ..............................................................................

cat("\n \n")
pheatmap( meanExpMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          scale = "row",
          annotation_row = data.frame(Markers = topMarkers_heatmapDF[ rownames( meanExpMat), "cluster"], stringsAsFactors = FALSE, row.names = rownames( meanExpMat) ),
          annotation_col = data.frame(Cluster = factor(colnames(meanExpMat), levels = levels(clusterID)), stringsAsFactors = FALSE, row.names = colnames(meanExpMat)),
          annotation_colors = list( Markers = clustersColor,
                                    Cluster = clustersColor),
          show_colnames = TRUE,
          fontsize_row = 5,
          main = "Scaled mean normalized expression\nby cluster of top markers");


# ..............................................................................
## @knitr heterogeneity_markerGenes_heatmap_dotplot
# ..............................................................................

cat("\n \n")
DotPlot( sc10x, 
         features = unique( topMarkersGenes),
         cols = "RdBu") +
  theme( axis.text.x = element_text( angle = 45, hjust = 1)) +
  ggtitle( "Scaled normalized mean expression of topmarker genes in clusters")



# ..............................................................................
## @knitr heterogeneity_markerGenes_heatmap_full
# ..............................................................................

# Plot a heatmap of the matrix of mean expression of top marker genes in clusters
cat("\n \n")
pheatmap( expMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(80),
          breaks = seq( -4.2, 4.2, 0.1),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          scale = "row",
          annotation_row = data.frame(Markers = topMarkers_heatmapDF[ rownames( meanExpMat), "cluster"], stringsAsFactors = FALSE, row.names = rownames( meanExpMat) ),
          annotation_col = data.frame(Cluster = clusterID, stringsAsFactors = FALSE, row.names = colnames( expMat)),
          annotation_colors = list( Markers = clustersColor,
                                    Cluster = clustersColor),
          show_colnames = FALSE,
          #fontsize_row = 5,
          main = "Scaled normalized expression by cell\nof top markers");



# ..............................................................................
## @knitr heterogeneity_markerGenes_heatmap_by_organ
# ..............................................................................

# Get the matrix of expression and associated clusters from Seurat object
expMat = as.matrix( GetAssayData( sc10x));
clusterID = Idents( sc10x);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkers_heatmapDF[["gene"]];
clusterOrdering = order( clusterID);

expMat = expMat[topMarkersGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Prefix marker gene name with cluster name to prevent eventual duplicated row names
rownames(expMat) = paste0( topMarkers_heatmapDF[["cluster"]], ".", rownames(expMat))

Tissue = sc10x[["HTOTISSUE_classification", drop = TRUE]][clusterOrdering] #HTOMICE_classification

meanExp_ByTissue_ByCluster = by( as.data.frame( t(expMat)),
                               list(Tissue=Tissue, Cluster=clusterID), 
                               colMeans)

# Convert as a list of matrices by Tissue, and bind as final matrix
meanExp_groupedTissue = apply(meanExp_ByTissue_ByCluster, 1, do.call, what=cbind, simplify = FALSE)
meanExp_ByTissue_ByCluster_matrix = do.call(cbind, meanExp_groupedTissue)

# Get Tissue name and Cluster name
TissueNames = rep(names(meanExp_groupedTissue), sapply(meanExp_groupedTissue, ncol))
ClusterNames = colnames(meanExp_ByTissue_ByCluster_matrix)

# Desambiguate row names in resulting matrix
#rownames(meanExp_ByTissue_ByCluster_matrix) = paste(TissueNames, rownames(meanExp_ByTissue_ByCluster_matrix), sep = " | ")

# Plot a heatmap of the matrix of mean expression of top marker genes in clusters
cat("\n \n")
# pheatmap( meanExp_ByTissue_ByCluster_matrix,
#           #color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(80)),
#           #breaks = seq( -4.2, 4.2, 0.1),
#           cluster_rows = FALSE,
#           cluster_cols = FALSE,
#           scale = "row",
#           #annotation_col = data.frame(Cluster = ClusterNames, 
#           #                            Tissue = TissueNames, 
#           #                            stringsAsFactors = FALSE,
#           #                            row.names = colnames(meanExp_ByTissue_ByCluster_matrix)), 
#           #annotation_colors = list( Markers = clustersColor,
#           #                          Cluster = clustersColor),
#           show_colnames = TRUE,
#           #fontsize_row = 5,
#           main = "Scaled normalized expression of\ntop markers by cluster & tissue");


Heatmap( t( scale( t( meanExp_ByTissue_ByCluster_matrix))), 
         rect_gp = gpar(col = "white", lwd =0.5),
         name = "Expression",
         col = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(3),
         cluster_columns = F, 
         cluster_rows = F,
         column_split = TissueNames)


# ..............................................................................
## @knitr heterogeneity_markerGenes_heatmap_by_mouse
# ..............................................................................

# Get the matrix of expression and associated clusters from Seurat object
expMat = as.matrix( GetAssayData( sc10x));
clusterID = Idents( sc10x);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkers_heatmapDF[["gene"]];
clusterOrdering = order( clusterID);

expMat = expMat[topMarkersGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Prefix marker gene name with cluster name to prevent eventual duplicated row names
rownames(expMat) = paste0( topMarkers_heatmapDF[["cluster"]], ".", rownames(expMat))

Mouse = sc10x[["HTOMICE_classification", drop = TRUE]][clusterOrdering] #HTOMICE_classification

meanExp_ByCluster_ByMouse = by( as.data.frame( t(expMat)),
                                list(Cluster=clusterID, Mouse=Mouse), 
                                colMeans)

# Convert as a list of matrices by Tissue, and bind alltogether as final matrix
meanExp_groupedMouse = apply(meanExp_ByCluster_ByMouse, 1, do.call, what=cbind, simplify = FALSE)
meanExp_ByCluster_ByMouse_matrix = do.call(cbind, meanExp_groupedMouse)

# Get Tissue name and Cluster name
ClusterNames = factor(rep(names(meanExp_groupedMouse), sapply(meanExp_groupedMouse, ncol)), levels = levels(Idents(sc10x)))
MouseNames = colnames(meanExp_ByCluster_ByMouse_matrix)

# Desambiguate row names in resulting matrix
#rownames(meanExp_ByTissue_ByCluster_matrix) = paste(TissueNames, rownames(meanExp_ByTissue_ByCluster_matrix), sep = " | ")

# Plot a heatmap of the matrix of mean expression of top marker genes in clusters
cat("\n \n")
# pheatmap( meanExp_ByTissue_ByCluster_matrix,
#           #color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(80)),
#           #breaks = seq( -4.2, 4.2, 0.1),
#           cluster_rows = FALSE,
#           cluster_cols = FALSE,
#           scale = "row",
#           #annotation_col = data.frame(Cluster = ClusterNames, 
#           #                            Tissue = TissueNames, 
#           #                            stringsAsFactors = FALSE,
#           #                            row.names = colnames(meanExp_ByTissue_ByCluster_matrix)), 
#           #annotation_colors = list( Markers = clustersColor,
#           #                          Cluster = clustersColor),
#           show_colnames = TRUE,
#           #fontsize_row = 5,
#           main = "Scaled normalized expression of\ntop markers by cluster & tissue");


Heatmap( t( scale( t( meanExp_ByCluster_ByMouse_matrix))), 
         rect_gp = gpar(col = "white", lwd =0.5),
         name = "Expression",
         col = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(3),
         cluster_columns = FALSE, 
         cluster_rows = FALSE,
         column_split = ClusterNames,
         row_split = topMarkers_heatmapDF[["cluster"]])



# ..............................................................................
## @knitr heterogeneity_markerGenes_heatmap_by_mouse_123
# ..............................................................................

sc10x_M123 = subset(sc10x, cells = which(sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3")))

# Get the matrix of expression and associated clusters from Seurat object
expMat = as.matrix( GetAssayData( sc10x_M123));
clusterID = Idents( sc10x_M123);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkers_heatmapDF[["gene"]];
clusterOrdering = order( clusterID);

expMat = expMat[topMarkersGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Prefix marker gene name with cluster name to prevent eventual duplicated row names
rownames(expMat) = paste0( topMarkers_heatmapDF[["cluster"]], ".", rownames(expMat))

Mouse = sc10x_M123[["HTOMICE_classification", drop = TRUE]][clusterOrdering] #HTOMICE_classification

meanExp_ByCluster_ByMouse = by( as.data.frame( t(expMat)),
                                list(Cluster=clusterID, Mouse=Mouse), 
                                colMeans)

# Convert as a list of matrices by Tissue, and bind alltogether as final matrix
meanExp_groupedMouse = apply(meanExp_ByCluster_ByMouse, 1, do.call, what=cbind, simplify = FALSE)
meanExp_ByCluster_ByMouse_matrix = do.call(cbind, meanExp_groupedMouse)

# Get Tissue name and Cluster name
ClusterNames = factor(rep(names(meanExp_groupedMouse), sapply(meanExp_groupedMouse, ncol)), levels = levels(Idents(sc10x)))
MouseNames = colnames(meanExp_ByCluster_ByMouse_matrix)

# Desambiguate row names in resulting matrix
#rownames(meanExp_ByTissue_ByCluster_matrix) = paste(TissueNames, rownames(meanExp_ByTissue_ByCluster_matrix), sep = " | ")

# Plot a heatmap of the matrix of mean expression of top marker genes in clusters
cat("\n \n")


Heatmap( t( scale( t( meanExp_ByCluster_ByMouse_matrix))), 
         rect_gp = gpar(col = "white", lwd =0.5),
         name = "Expression",
         col = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(3),
         cluster_columns = F, 
         cluster_rows = F,
         column_split = ClusterNames,
         row_split = topMarkers_heatmapDF[["cluster"]])




# ..............................................................................
## @knitr heterogeneity_markerGenes_heatmap_by_mouse_456
# ..............................................................................

sc10x_M456 = subset(sc10x, cells = which(sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M4", "M5", "M6")))

# Get the matrix of expression and associated clusters from Seurat object
expMat = as.matrix( GetAssayData( sc10x_M456));
clusterID = Idents( sc10x_M456);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkers_heatmapDF[["gene"]];
clusterOrdering = order( clusterID);

expMat = expMat[topMarkersGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Prefix marker gene name with cluster name to prevent eventual duplicated row names
rownames(expMat) = paste0( topMarkers_heatmapDF[["cluster"]], ".", rownames(expMat))

Mouse = sc10x_M456[["HTOMICE_classification", drop = TRUE]][clusterOrdering] #HTOMICE_classification

meanExp_ByCluster_ByMouse = by( as.data.frame( t(expMat)),
                                list(Cluster=clusterID, Mouse=Mouse), 
                                colMeans)

# Convert as a list of matrices by Tissue, and bind alltogether as final matrix
meanExp_groupedMouse = apply(meanExp_ByCluster_ByMouse, 1, do.call, what=cbind, simplify = FALSE)
meanExp_ByCluster_ByMouse_matrix = do.call(cbind, meanExp_groupedMouse)

# Get Tissue name and Cluster name
ClusterNames = factor(rep(names(meanExp_groupedMouse), sapply(meanExp_groupedMouse, ncol)), levels = levels(Idents(sc10x)))
MouseNames = colnames(meanExp_ByCluster_ByMouse_matrix)

# Desambiguate row names in resulting matrix
#rownames(meanExp_ByTissue_ByCluster_matrix) = paste(TissueNames, rownames(meanExp_ByTissue_ByCluster_matrix), sep = " | ")

# Plot a heatmap of the matrix of mean expression of top marker genes in clusters
cat("\n \n")


Heatmap( t( scale( t( meanExp_ByCluster_ByMouse_matrix))), 
         rect_gp = gpar(col = "white", lwd =0.5),
         name = "Expression",
         col = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(3),
         cluster_columns = F, 
         cluster_rows = F,
         column_split = ClusterNames,
         row_split = topMarkers_heatmapDF[["cluster"]])




# ..............................................................................
## @knitr heterogeneity_markerGenes_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)
# ..............................................................................

# Plot expression values of marker genes on dimreduc figures for each cluster (TODO: message if list empty)
invisible( lapply( names( topMarkers_heatmap), function(clusterName)
{
  cat("#### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");

  # Highlight cells of current cluster on a dimreduc plot
  highlightClusterPlot(clusterName, seuratObject = sc10x, reduction = ifelse( exists("useReduction"), useReduction, "umap"));

  # Plots expression on projected cells
  invisible( lapply( topMarkers_heatmap[[clusterName]][["gene"]], function(featureName)
    {
      print( FeaturePlot( sc10x, features = featureName, reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE) +
               theme( axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      legend.position = "none"));
    }));

  cat(" \n \n"); # Required for '.tabset'
}));




# ..............................................................................
## @knitr heterogeneity_markerGenes_expression_violin
# ..............................................................................

# Plot expression values of marker genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( topMarkers_heatmap), function(clusterName)
{
  cat("#### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");

  # Remind cluster name in an empty figure to keep consistent alignment of panels between tabs
  plot( c( 0, 1), c( 0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
  text( x = 0.5, y = 0.5, paste( "Cluster", clusterName), cex = 2, col = clustersColor[clusterName]);

  # Violinplot for expression value of marker genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( topMarkers_heatmap[[clusterName]][["gene"]], violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor));

  cat(" \n \n"); # Required for '.tabset'
}))

