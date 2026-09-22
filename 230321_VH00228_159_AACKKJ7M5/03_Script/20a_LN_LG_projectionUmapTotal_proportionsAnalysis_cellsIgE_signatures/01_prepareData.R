# #########################################
# This script reads and filters sc10x  data
# #########################################


###############
# CREATE OUTPUT
###############

# ..............................................................................
## @knitr createOutputFolders
# ..............................................................................

### Create output directory
message(paste("Creating 'output' folder:", PATH_ANALYSIS_OUTPUT))
dir.create( PATH_ANALYSIS_OUTPUT, recursive = TRUE, showWarnings = FALSE);

PATH_ANALYSIS_EXTRA_OUTPUT = file.path( PATH_ANALYSIS_OUTPUT, "Extra")
message(paste("Creating 'extra output' folder:", PATH_ANALYSIS_EXTRA_OUTPUT))
dir.create( PATH_ANALYSIS_EXTRA_OUTPUT, recursive = TRUE, showWarnings = FALSE);




###########
# READ DATA
###########

# ..............................................................................
## @knitr loadData
# ..............................................................................

#### Seurat object

# Load Seurat from previously saved binary RDS file (must contain numID, )
sc10x = readRDS(PATH_RDS_SEURAT_OBJECT);

cat( paste0( "\n<br>Successfuly loaded Seurat object: ", ncol( sc10x), " Cells x ", nrow( sc10x)," Genes."));


# Save current Ident as a proper metadata column
sc10x[["Identity"]] = Idents(sc10x)

# Filtering cells with Negative HTO (first or second)
negative_hto = grepl("Negative", sc10x[["factorHTO", drop = TRUE]])
if(any(negative_hto))
{
  cat( paste0( "\n<br>Removing ", sum( negative_hto), " cells with negative HTO..."));
  sc10x = sc10x[,!negative_hto]
  cat( paste0( "\n<br>Using Seurat object: ", ncol( sc10x), " Cells x ", nrow( sc10x)," Genes."));
  sc10x[["factorHTO"]] = factor(sc10x[["factorHTO", drop = TRUE]]) 
}

# Extract coordinates of interest (umap), as opposed to other scripts allowing
# external coordinates or other dimreduc
cellsCoordinates = Embeddings( sc10x, reduction = "umap")


# Simplify cluster names for figures layout (remove parenthesis information)
# NOTE : one should take this in account in analysisParams value 'ORDER_CLUSTER'
if(SIMPLIFY_CLUSTER_NAMES) 
{
  levels(Idents(sc10x)) = gsub( " \\(.*\\)", "", levels(Idents(sc10x)))
}


# Eventual reordering of cluster factor for figures
if(!is.null( ORDER_CLUSTER)) 
{
  Idents(sc10x) = factor(Idents(sc10x), levels = ORDER_CLUSTER)
}



# Read eventual TSV file defining a color for each cluster (attribute ggplot defaults otherwise)
clustersColor = NULL;
if(file.exists( EXTERNAL_CLUSTERSCOLOR_PATH))
{
  clustersColor =  as.matrix( read.table( EXTERNAL_CLUSTERSCOLOR_PATH, 
                                          sep = "\t", 
                                          header = FALSE, 
                                          comment.char= "", 
                                          row.names = 1, 
                                          stringsAsFactors=FALSE))[,1]; # Convert to matrix to get a named vector when extracting column
  
  
  # Check consistency between colors file and clusters
  clusterInColorFile = levels( Idents( sc10x)) %in% names( clustersColor);
  if(!all( clusterInColorFile)) stop( paste0( "Following cluster name(s) could not be found in file defining clusters color: ", paste( levels( Idents( sc10x))[!clusterInColorFile], collapse = " - "), "."));
  
  clusterNameExists = names( clustersColor) %in% levels( Idents( sc10x));
  if(!all( clusterNameExists)) warning( paste0( "Following cluster name from file defining colors could not be found in loaded data: ", paste( names( clustersColor)[!clusterNameExists], collapse = " - "), "."));
  
} else
{
  # Define a set of colors for clusters (based on ggplot default)
  clustersColor = hue_pal()( nlevels( Idents( sc10x)));
  names( clustersColor) = levels( Idents( sc10x));
}




# Prepare the metadata column (as factor) that will be used for HTO reference
if(!is.null(HTO_FACTOR_LEVELS)) # Eventually reorder levels for plots (needs a separate call as setting 'levels' to null is not what we want)
{
  sc10x[["factorHTO"]] = factor( sc10x[[HTO_METADATA_COLUMN, drop = TRUE]], levels = HTO_FACTOR_LEVELS)
} else
{
  sc10x[["factorHTO"]] = factor( sc10x[[HTO_METADATA_COLUMN, drop = TRUE]])
}



# Read eventual TSV file defining a color for each HTO (attribute ggplot defaults otherwise)
HTOsColor = NULL;
if(file.exists( EXTERNAL_HTOSCOLOR_PATH))
{
  HTOsColor =  as.matrix( read.table( EXTERNAL_HTOSCOLOR_PATH, 
                                      sep = "\t", 
                                      header = FALSE, 
                                      comment.char= "", 
                                      row.names = 1, 
                                      stringsAsFactors=FALSE))[,1]; # Convert to matrix to get a named vector when extracting column
  
  
  # Check consistency between colors file and HTOs
  HTOsInColorFile = levels( sc10x[["factorHTO"]]) %in% names( HTOsColor);
  if(!all( HTOsInColorFile)) stop( paste0( "Following HTO name(s) could not be found in file defining HTOs color: ", paste( levels( sc10x[["factorHTO", drop = TRUE]])[!HTOsInColorFile], collapse = " - "), "."));
  
  HTOsNameExists = names(  HTOsColor) %in% levels( sc10x[["factorHTO", drop = TRUE]]);
  if(!all( HTOsNameExists)) warning( paste0( "Following HTO name from file defining colors could not be found in loaded data: ", paste( names( HTOsColor)[!HTOsNameExists], collapse = " - "), "."));
  
} else
{
  # Define a set of colors for HTOs (based on ggplot default)
  HTOsColor = hue_pal()( nlevels( sc10x[["factorHTO", drop = TRUE]]));
  names( HTOsColor) = levels( sc10x[["factorHTO", drop = TRUE]]);
}



# !!! MANUAL OVERRIDE to define colors for HTOs !!!
#HTOsColor = tail(RColorBrewer::brewer.pal(10, "Paired"), 4)
#names(HTOsColor) = HTO_FACTOR_LEVELS



#### Make some expression and counts computations used in further analyses 

# Compute the number of cells by cluster
clustersCount = as.data.frame( table( Cluster = Idents(sc10x)), responseName = "Total");

# Compute the cluster counts sliced by HTO condition
clustersCountByHTO = table(cbind("Identity" = Idents( sc10x), "HTO" = sc10x[["factorHTO"]]))

# Compute a matrix of average expression value by cluster (for each gene)
geneExpByCluster = do.call( rbind, 
                            apply( as.matrix( GetAssayData( sc10x)), # expression values
                                   1,                                # by rows
                                   tapply,                           # apply by group
                                   INDEX = Idents( sc10x),           # clusters IDs
                                   mean,                             # summary function
                                   simplify = FALSE));               # do not auto coerce

# # Save it as 'tsv' file
# write.table( geneExpByCluster, 
#              file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "normExpressionByCluster.tsv")), 
#              quote = FALSE, 
#              row.names = TRUE, 
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");

# Compute a matrix of average expression value by cluster (for each gene)
geneExpByHTO = do.call( rbind, 
                        apply( as.matrix( GetAssayData( sc10x, 
                                                        slot = "data", 
                                                        assay="RNA")),  # normalized expression values
                               1,                                       # by rows
                               tapply,                                  # apply by group
                               INDEX = sc10x[["factorHTO", drop = TRUE]],                  # clusters IDs
                               mean,                                    # summary function
                               simplify = FALSE));                      # do not auto coerce

# # Save it as 'tsv' file
# write.table( geneExpByHTO, 
#              file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "normExpressionByHTO.tsv")), 
#              quote = FALSE, 
#              row.names = TRUE, 
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");

# Compute a matrix of average expression value by cluster and HTOs (for each gene)
geneExpByClusterAndHTO = do.call( rbind, 
                                  apply( as.matrix( GetAssayData( sc10x, 
                                                                  slot = "data", 
                                                                  assay="RNA")),  # normalized expression values
                                         1,                                       # by rows
                                         tapply,                                  # apply by group
                                         INDEX = paste( Idents( sc10x),
                                                        sc10x[["factorHTO", drop = TRUE]], 
                                                        sep = "_"),               # clusters + HTOs IDs
                                         mean,                                    # summary function
                                         simplify = FALSE));                      # do not auto coerce
# 
# # Save it as 'tsv' file
# write.table( geneExpByClusterAndHTO, 
#              file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "normExpressionByClusterAndHTO.tsv")), 
#              quote = FALSE, 
#              row.names = TRUE, 
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");
# 
# 
# 

########################################
# CHECK MODULES GENES LISTS 
########################################

modulesGenesNotFound = list()
matchModulesGenes = list()
for(currentModule in names(MODULES_GENES))
{
  ### Remove eventual NULL (empty) list elements from list of genes in modules
  modulesGroupEmpty = sapply( MODULES_GENES[[currentModule]], is.null);
  if(any( modulesGroupEmpty)) warning( paste("Following module(s) of genes will be ignored because empty:", paste( names(modulesGroupEmpty)[modulesGroupEmpty], collapse=" - ")));
  MODULES_GENES[[currentModule]] = MODULES_GENES[[currentModule]][! modulesGroupEmpty];
  
  # Check whether genes in MODULES_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
  #matchModulesGenes[[currentModule]] = match( toupper( unlist( MODULES_GENES)), toupper( rownames( GetAssayData( sc10x))));
  matchModulesGenes[[currentModule]] = match( ( unlist( MODULES_GENES[[currentModule]])), ( rownames( GetAssayData( sc10x))));
  modulesGenesNotFound[[currentModule]] = unique( unlist( MODULES_GENES[[currentModule]])[is.na( matchModulesGenes[[currentModule]])]);
  if(any( is.na( matchModulesGenes[[currentModule]]))) warning( paste0( "Following gene(s) from modules list '", currentModule, "' could not be found in experimental data and will be ignored: ", paste( paste0("'", modulesGenesNotFound[[currentModule]], "'"), collapse=" - ")));
  # Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
  MODULES_GENES[[currentModule]] = relist( rownames( GetAssayData( sc10x))[ matchModulesGenes[[currentModule]] ], skeleton = MODULES_GENES[[currentModule]]); # Does not work with NULL list elements (removed earlier)
  # Finally remove names that did not match
  MODULES_GENES[[currentModule]] = lapply( MODULES_GENES[[currentModule]], na.omit);
  
  ### Remove eventual NULL (empty) list elements from list of genes in modules
  modulesGroupEmpty = sapply( MODULES_GENES[[currentModule]], is.null);
  if(any( modulesGroupEmpty)) warning( paste("Following module(s) of genes will be ignored because empty:", paste( names(modulesGroupEmpty)[modulesGroupEmpty], collapse=" - ")));
  MODULES_GENES[[currentModule]] = MODULES_GENES[[currentModule]][! modulesGroupEmpty];
}

