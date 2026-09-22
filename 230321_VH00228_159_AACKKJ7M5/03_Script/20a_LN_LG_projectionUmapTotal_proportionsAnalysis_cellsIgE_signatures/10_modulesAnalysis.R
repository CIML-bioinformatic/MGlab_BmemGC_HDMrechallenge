

# MODULES ANALYSIS
##################

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
                    axis.title.y = element_blank(),
                    legend.position = "none"));
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
                    axis.title.y = element_blank(),
                    legend.position = "none"));
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




## @knitr heterogeneity_modules_expression_projection_pngFile
# Plot expression values of individual module genes as png files (not in report)
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Compute 'pixel' dimensions to match figures in report (based on chunk params)
defaultDim = c( 800, 800);  # Fallback values if chunk params not available
defaultDpi = 72;            # (e.g. executing directly from R)

chunkDim    = knitr::opts_current$get( "fig.dim");
chunkDpi    = knitr::opts_current$get( "dpi");

figDim = if( is.null( chunkDim) || is.null( chunkDpi) ) defaultDim else chunkDim * chunkDpi;
figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;


invisible( lapply( names( MODULES_GENES), function(moduleName)
{
  message(moduleName); # Just for tracking progress in console

  # Create subfolder for current module to store all png files
  pathCurrentModule = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "modulesExpressionIndividual_", ifelse(exists("useReduction"), useReduction, "umap")), moduleName);
  dir.create(pathCurrentModule, showWarnings = FALSE, recursive = TRUE);

  # Plots expression on projected cells (or error message if feature not found)
  invisible( lapply( MODULES_GENES[[moduleName]], function(featureName)
  {
    # Create png file
    png( file.path( pathCurrentModule, paste0( featureName, ".png")), 
         width = figDim[1],
         height = figDim[2],
         res = figDpi);

    print(FeaturePlot( sc10x, features = featureName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE)  +
      theme( axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none"))

    dev.off(); # Close file descriptor
  }));

}));



## @knitr heterogeneity_modules_expression_violin_pngFile
# Plot expression values of monitored genes as violinplot in a png file for each cluster (TODO: message if list empty)

# Compute 'pixel' dimensions to match figures in report (based on chunk params)
defaultDim = c( 800, 800);  # Fallback values if chunk params not available
defaultDpi = 72;            # (e.g. executing directly from R)

chunkDim    = knitr::opts_current$get( "fig.dim");
chunkDpi    = knitr::opts_current$get( "dpi");

figDim = if( is.null( chunkDim) || is.null( chunkDpi) ) defaultDim else chunkDim * chunkDpi;
figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;


invisible( lapply( names( MODULES_GENES), function(moduleName)
{
  message(moduleName); # Just for tracking progress in console

  # Create subfolder for current module to store all png files
  pathCurrentModule = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "modulesExpressionIndividual_violin"), moduleName);
  dir.create(pathCurrentModule, showWarnings = FALSE, recursive = TRUE);

  # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( MODULES_GENES[[moduleName]], function(featureName)
  {
    # Create png file
    png( file.path( pathCurrentModule, paste0( featureName, ".png")), 
         width = figDim[1],
         height = figDim[2],
         res = figDpi);

    violinFeatureByCluster(featureName, seuratObject = sc10x, clustersColor = clustersColor);

    dev.off(); # Close file descriptor
  }))

}))


