# #########################################################
# This script creates figures in high resolution with final
# names and colors for publication
# #########################################################

# UMAP
######


# cat(" \n \n"); # Required for '.tabset'






# MARKER GENES
##############

# Mostly imported from 14c 07_cellHeterogeneity_MarkerGenes.R
# variable 'markers' replaced by reading from 14c result (csv) 

# ..............................................................................
## @knitr heterogeneity_markerGenes
# ..............................................................................

# Identify marker genes
markers = read.csv(PATH_MARKERS_LIST, sep = '\t', row.names = 1)

# Save markers list as 'tsv' table
# write.table( markers,
#              file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "MarkerGenes.tsv")),
#              quote = FALSE,
#              row.names = TRUE, 
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");


# Update cluster names to last nomenclature 
markers[["cluster"]] = factor(markers[["cluster"]])

# Update cluster names to last nomenclature 
# NOTE : one should take this in account in analysisParams value 'ORDER_CLUSTER'
if(!is.null(CLUSTER_FACTOR_RECODE)) 
{
  markers[["cluster"]] = do.call(fct_recode, c(list('.f'=markers[["cluster"]]), CLUSTER_FACTOR_RECODE))
}

# Eventual reordering of cluster factor for figures
if(!is.null( ORDER_CLUSTER)) 
{
  markers[["cluster"]] = factor(markers[["cluster"]], levels = ORDER_CLUSTER)
}




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

# Same filters used to filter variable genes in previous steps 
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





# MODULES / SIGNATURES ANALYSIS (STRATIFIED)
############################################

## @knitr heterogeneity_modules

# Just remind the warning for genes names not in object, or modules that were transfered to individual monitoring of genes
for(currentModule in names(MODULES_GENES))
{
  if(any( is.na( matchModulesGenes[[currentModule]])))
  {
    warning( paste0( "Following gene(s) from modules list '", currentModule, "' could not be found in experimental data and will be ignored: ",
                     paste( modulesGenesNotFound[[currentModule]], collapse=" - ")));
  }
}



## @knitr heterogeneity_modules_scoring

# Compute the score of the cells according to group of monitored genes
for(currentModule in names(MODULES_GENES))
{
  for( listName in names( MODULES_GENES[[currentModule]]))
  {
    message(paste(currentModule, listName, sep="_"));
    if( length( MODULES_GENES[[currentModule]][[listName]]) == 0)
    {
      warning( paste0( "List of genes in module '", listName, "' is empty, ignoring..."));
    } else
    {
      sc10x <- AddModuleScore( object = sc10x,
                               features = MODULES_GENES[[currentModule]][ listName], # Must be a list
                               ctrl = MODULES_CONTROL_SIZE,           #length(MODULES_GENES[[ listName]]),
                               name = paste(currentModule, listName, sep="_"),
                               seed = SEED);
    }
  }
}



## @knitr heterogeneity_modules_dotplot

# Define an eventual subselection on object
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6") &
#                                                sc10x[["HTOTISSUE_classification", drop = TRUE]] %in% c("LN", "LG")])

for(currentModule in names(MODULES_GENES))
{
  cat("\n#### ", currentModule, "\n", sep="")
  p = DotPlot(scObject, 
              features = paste0(paste(currentModule, names(MODULES_GENES[[currentModule]]), sep="_"), 1), 
              cols = "RdBu", dot.scale = 15 ) + 
    coord_flip() +
    scale_y_discrete(guide = guide_axis(n.dodge = 2)) # Dodge x axis labels
  
  print(p)
  
  cat(" \n \n"); # Required for '.tabset'
  
}
cat(" \n \n"); # Required for '.tabset'



## @knitr heterogeneity_modules_dotplot_splitTissue

# Define an eventual subselection on object
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6") &
#                                                sc10x[["HTOTISSUE_classification", drop = TRUE]] %in% c("LN", "LG")])

# Override temporary "identity", splitting clusters by Tissue 
Idents(scObject) = factor(paste(scObject[["HTOTISSUE_classification", drop = TRUE]], 
                                scObject[["Identity", drop = TRUE]], 
                                sep = "_"),
                          levels = paste(c("LG", "LN"), rep(ORDER_CLUSTER, each=2), sep = "_"))

for(currentModule in names(MODULES_GENES))
{
  cat("\n#### ", currentModule, "\n", sep="")
  
  p = DotPlot(scObject, 
              features = paste0(paste(currentModule, names(MODULES_GENES[[currentModule]]), sep="_"), 1), 
              cols = "RdBu", dot.scale = 15 ) + 
    coord_flip() +
    scale_y_discrete(guide = guide_axis(n.dodge = 2)) # Dodge x axis labels
  
  print(p)
  
  cat(" \n \n"); # Required for '.tabset'
  
}
cat(" \n \n"); # Required for '.tabset'



## @knitr heterogeneity_modules_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Define an eventual subselection on object
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6")])

# Plot scores of modules on dimensionality reduction figures
for(currentModule in names(MODULES_GENES))
{
  cat("\n#### ", currentModule, "\n", sep="")
  
  invisible( lapply( names(MODULES_GENES[[currentModule]]), function(listName)
  {
    print( FeaturePlot(scObject, features = paste0(paste(currentModule, listName, sep="_"), "1"), reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE) +
             #ggtitle( label = listName) +
             theme( axis.title.x = element_blank(),
                    axis.title.y = element_blank() ));
                    #legend.position = "none"));
  }));
  cat(" \n \n"); # Required for '.tabset'
  
}
cat(" \n \n"); # Required for '.tabset'



## @knitr heterogeneity_modules_expression_projection_splitTissue
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Define an eventual subselection on object (only Mice, Use 'split.by' for Tissue to keep a consistent color scale)
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6")])

# Plot scores of modules on dimensionality reduction figures
for(currentModule in names(MODULES_GENES))
{
  cat("\n#### ", currentModule, "\n", sep="")
  
  invisible( lapply( names(MODULES_GENES[[currentModule]]), function(listName)
  {
    print( FeaturePlot(scObject, features = paste0(paste(currentModule, listName, sep="_"), "1"), reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE, split.by = "HTOTISSUE_classification") +
             #ggtitle( label = listName) +
             theme( axis.title.x = element_blank(),
                    axis.title.y = element_blank(),# ));
                    legend.position = "right"));
  }));
  cat(" \n \n"); # Required for '.tabset'
  
}
cat(" \n \n"); # Required for '.tabset'




## @knitr heterogeneity_modules_expression_violin
# Violinplot for module values by cluster

# Define an eventual subselection on object
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6") &
#                                              sc10x[["HTOTISSUE_classification", drop = TRUE]] %in% c("LN", "LG")])

# Violinplot of modules
for(currentModule in names(MODULES_GENES))
{
  cat("\n#### ", currentModule, "\n", sep="")
  
  invisible( lapply( paste0(paste(currentModule, names(MODULES_GENES[[currentModule]]), sep="_"), 1), 
                     violinFeatureByCluster, 
                     seuratObject = scObject, 
                     yJitter = 0,
                     clustersColor = clustersColor, 
                     yLabel = "Score", 
                     addStats = FALSE, 
                     trimTitle = 1, 
                     symmetricY=TRUE));
  cat(" \n \n"); # Required for '.tabset'
  
}
cat(" \n \n"); # Required for '.tabset'







## @knitr heterogeneity_modules_expression_violin_splitTissue
# Violinplot for module values by cluster

# Define an eventual subselection on object
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6") &
#                                                sc10x[["HTOTISSUE_classification", drop = TRUE]] %in% c("LN", "LG")])

# Override temporary "identity", splitting clusters by Tissue 
Idents(scObject) = factor(paste(scObject[["HTOTISSUE_classification", drop = TRUE]], 
                                scObject[["Identity", drop = TRUE]], 
                                sep = "_"),
                          levels = paste(c("LG", "LN"), rep(ORDER_CLUSTER, each=2), sep = "_"))
# And define colors based on cluster and Tissue
clustersColor_byTissue = rep(clustersColor, each=2)
names(clustersColor_byTissue) = paste(c("LG", "LN"), rep(names(clustersColor), each=2), sep = "_")


# Violinplot of modules
for(currentModule in names(MODULES_GENES))
{
  cat("\n#### ", currentModule, "\n", sep="")
  
  invisible( lapply( paste0(paste(currentModule, names(MODULES_GENES[[currentModule]]), sep="_"), 1), 
                     violinFeatureByCluster, 
                     seuratObject = scObject, 
                     yJitter = 0,
                     clustersColor = clustersColor_byTissue, 
                     yLabel = "Score", 
                     addStats = FALSE, 
                     trimTitle = 1, 
                     symmetricY=TRUE));
  cat(" \n \n"); # Required for '.tabset'
  
}
cat(" \n \n"); # Required for '.tabset'
















# MONITORED GENES
#################

## @knitr heterogeneity_monitoredGenes

# Just remind the warning for genes names not in object
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));



# 
# ## @knitr heterogeneity_monitoredGenes_heatmap
# 
# # DoHeatmap replaced by use of iHeatmapr after testing several options
# #htmltools::tagList( ggplotly( DoHeatmap(object = sc10x, features = unlist(MONITORED_GENES))));
# 
# # Get the matrix of expression and associated clusters from Seurat object
# expMat = as.matrix( GetAssayData( sc10x));
# clusterID = Idents( sc10x);
# 
# # Select monitored genes and reorder cells to group clusters together
# monitoredGenes = unlist( MONITORED_GENES);
# clusterOrdering = order( clusterID);
# 
# expMat = expMat[monitoredGenes, clusterOrdering];
# clusterID = clusterID[ clusterOrdering];
# 
# # Prepare rows and columns annotation bars (monitored group and cluster respectively)
# rowsAnnot = data.frame( Monitored = fct_rev(fct_inorder(rep( names( MONITORED_GENES), sapply( MONITORED_GENES, length)))));
# colsAnnot = data.frame( Cluster = clusterID);
# 
# # Prepare unique rows and cols names for pheatmap (annotation rows) and match with rowAnnots and colAnnots row names
# originalRowNames = rownames( expMat);
# originalColNames = colnames( expMat);
# rownames( expMat) = make.unique( originalRowNames);
# colnames( expMat) = make.unique( originalColNames);
# rownames( rowsAnnot) = rownames( expMat);
# rownames( colsAnnot) = colnames( expMat);
# 
# # Prepare colors of monitored genes groups for pheatmap (requires named vector matching data factor levels)
# monitoredColors = rainbow( nlevels( rowsAnnot[["Monitored"]]), s = 0.8);
# names( monitoredColors) = levels( rowsAnnot[["Monitored"]]);
# 
# pheatmap( expMat,
#           color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
#           cluster_rows = FALSE,
#           cluster_cols = FALSE,
#           annotation_row = rowsAnnot,
#           annotation_col = colsAnnot,
#           labels_row = originalRowNames,
#           annotation_colors = list( Monitored = monitoredColors,
#                                     Cluster = clustersColor),
#           show_colnames = FALSE);


## @knitr heterogeneity_monitoredGenes_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Define an eventual subselection on object
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6") &
#                                              sc10x[["HTOTISSUE_classification", drop = TRUE]] %in% c("LN", "LG")])

# Plot expression values of monitored genes (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("####", monitoredGroup, "\n");

  # Plots expression on projected cells (or error message if feature not found)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], function(featureName)
  {
    print(
      tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
                FeaturePlot( scObject, features = featureName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE) +
                  theme( axis.title.x = element_blank(),
                         axis.title.y = element_blank() )),
                         #legend.position = "none")),
                error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
  }));

  cat(" \n \n"); # Required for '.tabset'
}));



## @knitr heterogeneity_monitoredGenes_expression_projection_splitTissue
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Define an eventual subselection on object (only Mice, Use 'split.by' for Tissue to keep a consistent color scale)
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6")])


# Plot expression values of monitored genes (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("####", monitoredGroup, "\n");
  
  # Plots expression on projected cells (or error message if feature not found)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], function(featureName)
  {
    print(
      tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
        FeaturePlot(scObject, features = featureName, reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE, split.by = "HTOTISSUE_classification") +
          theme( axis.title.x = element_blank(),
                 axis.title.y = element_blank(),# ));
                 legend.position = "right")),
        error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
  }));
  
  cat(" \n \n"); # Required for '.tabset'
}));



## @knitr heterogeneity_monitoredGenes_expression_violin

# Define an eventual subselection on object
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6") &
#                                              sc10x[["HTOTISSUE_classification", drop = TRUE]] %in% c("LN", "LG")])

# Plot expression values of monitored genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("####", monitoredGroup, "\n");

  # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], violinFeatureByCluster, seuratObject = scObject, slot = "data", clustersColor = clustersColor));

  cat(" \n \n"); # Required for '.tabset'
}));


## @knitr heterogeneity_monitoredGenes_expression_violin_splitTissue
# Violinplot for module values by cluster

# Define an eventual subselection on object
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6") &
#                                                sc10x[["HTOTISSUE_classification", drop = TRUE]] %in% c("LN", "LG")])

# Override temporary "identity", splitting clusters by Tissue 
Idents(scObject) = factor(paste(scObject[["HTOTISSUE_classification", drop = TRUE]], 
                                scObject[["Identity", drop = TRUE]], 
                                sep = "_"),
                          levels = paste(c("LG", "LN"), rep(ORDER_CLUSTER, each=2), sep = "_"))
# And define colors based on cluster and Tissue
clustersColor_byTissue = rep(clustersColor, each=2)
names(clustersColor_byTissue) = paste(c("LG", "LN"), rep(names(clustersColor), each=2), sep = "_")

# Plot expression values of monitored genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("####", monitoredGroup, "\n");
  
  # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], violinFeatureByCluster, seuratObject = scObject, slot = "data", clustersColor = clustersColor_byTissue));
  
  cat(" \n \n"); # Required for '.tabset'
}));







# SCATTER_GENES 
##################

## @knitr scatter_genes_prepare
# Create a scatterplot of cells for two selected features
# Can define parentLevel as the number of '#' in parent rmarkdown title so the 
# loop adds paragraph title in the correct (parentLevel+1) level.

# Function creating a scatterplot of selected features scores colored by identity.
# Returns a list with "scatter" ggplot object, and a "tablePos" dataframe giving count and pct of positive cells for each identity
# DEBUG: seuratObject = sc10x;featureNames = SCATTER_GENES[[1]]; assay = "RNA"; slot = "data"; main = ""
genesTableScatter = function(seuratObject, featureNames, assay = "RNA", slot = "data", main = "")
{
  expScores = as.data.frame( t( GetAssayData( object = seuratObject,
                                              slot = slot,
                                              assay = assay)[featureNames[1:2],]))
  if(!(dim(expScores)[2]==2)) stop(paste("Could not retrieve expression values to create scatterplot for selected features:", paste(featureNames, collapse = " - ") ))
  
  #expPos = expScores>0
  #colnames(expPos) = paste0(colnames(expScores), "_pos")
  
  dataPlot = cbind(expScores, ident = Idents(seuratObject))
  
  scatter = ggplot( data = dataPlot,
                    aes( x = !!sym(colnames(expScores)[1]),
                         y = !!sym(colnames(expScores)[2]),
                         color = ident,
                         fill = ident)) +
    geom_point(alpha=0.3) +
    scale_color_manual( values = clustersColor) +
    scale_fill_manual( values = clustersColor) +
    ggtitle(main)
  
  # Compute counts of positive cells by identity
  identSorted = factor(Idents(seuratObject), levels = sort(levels(Idents(sc10x))))  
  posCount = do.call(rbind, by(expScores, identSorted, function(x){colSums(x>0)}))
  
  # Count nb of cells by identity
  ident_nbcells = as.numeric(table(identSorted)[rownames(posCount)])
  
  # Use previous results to compute percentage of positive cell of each identity
  posPct = round((posCount/ident_nbcells)*100, digits = 2)
  colnames(posPct) = paste0(colnames(posCount), "_pct_ident")
  
  # Combine stats in a single table
  tablePos = cbind(posCount, posPct, ident_nbcells)
  
  #ggplot( data = reshape2::melt(posCount),
  #        aes( x = "", y = value, fill = as.character(Var1))) +
  #  geom_bar(stat = "identity") +
  #  coord_polar("y", start = 0) +
  #  facet_grid(rows = Var2, 
  
  #ggMarginal( p,
  #            color = NA,
  #            groupColour = FALSE,
  #            groupFill = TRUE,
  #            type = "histogram",
  #            alpha = 1)
  
  return(list("scatter" = scatter,
              "tablePos" = tablePos))
}


## @knitr scatter_genes_plot_and_table

# Define an eventual subselection on object
#scObject = subset(sc10x, cells = Cells(sc10x)[sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M1", "M2", "M3", "M4", "M5", "M6") &
#                                              sc10x[["HTOTISSUE_classification", drop = TRUE]] %in% c("LN", "LG")])


# For each list of 2 genes, create a scatterplot and a table in tabs
for(scatterName in names(SCATTER_GENES))
{
  cat("\n\n", rep('#', if(exists("parentLevel")) parentLevel+1 else 2), " ", scatterName, " {.tabset .tabset-pills .tabset-fade}\n", sep = "")
  message( paste("Scatter comparison:", scatterName))
  
  ts = genesTableScatter( seuratObject = scObject,
                          featureNames = SCATTER_GENES[[scatterName]],
                          main = scatterName)
  
  cat("\n\n", rep('#', if(exists("parentLevel")) parentLevel+2 else 3), " Scatterplot\n", sep = "")
  print(ts[["scatter"]] + theme_classic())
  
  cat("\n\n", rep('#', if(exists("parentLevel")) parentLevel+2 else 3), " By cluster\n", sep = "")
  print(ts[["scatter"]] + 
          #facet_wrap(vars(factor(ident, levels = sort(levels(ident))))) +
          facet_wrap(vars(factor(ident, levels = levels(Idents(scObject))))) +
          theme_pubr() +
          theme(legend.position = "none"))
  
  cat("\n\n", rep('#', if(exists("parentLevel")) parentLevel+2 else 3), " Table\n", sep = "")
  print( kbl(ts[["tablePos"]], caption = paste0(scatterName, " - positive cells")) %>%
           kable_paper("hover", full_width = F))
  cat("\n")
}


