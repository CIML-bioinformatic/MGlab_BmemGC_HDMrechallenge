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




#### External data files

# Read eventual TSV file to define groups of interest (and override eventual 
# clustering contained in seurat object).
if(file.exists( EXTERNAL_CLUSTERING_PATH))
{
  externalClusters = read.csv( EXTERNAL_CLUSTERING_PATH, sep = "\t", row.names = 1, stringsAsFactors = TRUE);
  
  if(!any( rownames( externalClusters) %in% colnames( sc10x)))
  {
    stop("Cells barcode from initial clustering TSV file do not match any barcode in Seurat object...");
  }
  
  if(!all( rownames( externalClusters) %in% colnames( sc10x))) 
  {
    warning( "Some cells barcode from initial clustering TSV file are not found in Seurat object...");
  }
  
  if(!all( colnames( sc10x) %in% rownames( externalClusters) ))
  {
    warning( "Some barcodes from Seurat object are not found in TSV file for initial clustering, they will be grouped as 'unknown' class...");
    # Restitute all barcodes of seurat object in 'externalClusters' and atribute them level 'unknown'
    externalClusters = externalClusters[colnames( sc10x),]; # Not found are NA
    rownames(externalClusters) = colnames( sc10x);
    externalClusters[["identity"]] = fct_explicit_na( externalClusters[["identity"]], 
                                                      na_level = "unknown"); # Replace NA factor values by an actual level
  }
  
  # Sort and select cells identity according to seurat object content
  externalClusters = externalClusters[colnames( sc10x),];
  
  # Set cells identity directly into Seurat object 
  Idents( sc10x) = factor(externalClusters[["identity"]]);
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




# Read eventual TSV file to define cells coordinates
cellsCoordinates = NULL; # This will be a matrix/df of cell coordinates used for custom/manual dimreduc plots (plotly/ggplot)
if(file.exists( CELLS_COORDINATES))
{
  cellsCoordinates = read.csv( CELLS_COORDINATES, sep = "\t", row.names = 1, stringsAsFactors = TRUE);
  
  if(!all( colnames( sc10x) %in% rownames( cellsCoordinates) ))
  {
    stop( "Some barcodes from Seurat object are not found in TSV file for reference coordinates...");
  }
  
  if(!all( rownames( cellsCoordinates) %in% colnames( sc10x))) 
  {
    warning( "Some barcodes from cells coordinates TSV file are not found in Seurat object...");
  }
  
  # Sort and select cells identity according to seurat object content
  cellsCoordinates = cellsCoordinates[colnames( sc10x),];
  
  # Update the Seurat object with a new Reduction slot for plotting with Seurat methods
  sc10x@reductions$external = Seurat::CreateDimReducObject( key = "EXT_", 
                                                            embeddings = cellsCoordinates)
  # Instruct plotting chunks to use this slot (they check for existence of this variable)
  useReduction = "external"
  
} else if((!RECOMPUTE_PREPROCESSING) && CELLS_COORDINATES %in% Reductions( sc10x) ) # We want to use the precomputed dimreduc
{
  cellsCoordinates = Embeddings(Reductions(sc10x, slot = CELLS_COORDINATES))
  
  # PCA can have arbitrary number of dimensions, use requested ones (2d plot)
  if(CELLS_COORDINATES == "pca")
  {
    cellsCoordinates = cellsCoordinates[, PCA_DIMS[1:2]];
  }
  
  # Instruct plotting chunks to use this slot (they check for existence of this variable)
  useReduction = CELLS_COORDINATES
  
} else if(!(RECOMPUTE_PREPROCESSING && (CELLS_COORDINATES %in% c("pca", "tsne", "umap")))) 
{
  stop("Parameter 'CELLS_COORDINATES' must refer to a DimReduc stored in loaded Seurat object, or contain a path to a valid file containing 2D cells coordinates, or be either 'pca', 'tsne', 'umap' when RECOMPUTE_PREPROCESSING is TRUE...");
}
# At this stage 'cellsCoordinates' is either a set of coordinates (from TSV or previously computed object), or NULL (to be replaced by newly computed values)



#### Eventual filtering of cells to exclude from further analyses

cellsToRemove = logical( ncol( sc10x));

# By cluter name
if(! is.null( FILTER_CLUSTERS))
{
  cat( "\n<br><br>Selecting cluster(s) to filter before further analyses...")
  # Process each cluster name separately (not using %in%) to report specific warning messages
  for(currentFilter in FILTER_CLUSTERS)
  {
    cat( paste0( "\n<br>Cluster '", currentFilter, "'... "));
    cellsSelection = as.character( Idents( sc10x)) == as.character( currentFilter);
    cat( paste0( sum( cellsSelection) , " cells."));
    
    if(! any( cellsSelection)) warning( paste0( "Argument FILTER_CLUSTERS is specified but value '", currentFilter, "' did not identify any cell from Seurat object..."));
    cellsToRemove = cellsToRemove | cellsSelection;
  }
}


# Using a 'csv' file containing cell(s) barcode in a column named 'Cell'
if(! is.null( FILE_FILTER_CELLS) && all( file.exists( FILE_FILTER_CELLS)))
{
  cat( "\n<br><br>Selecting individual cell(s) to filter (from file) before further analyses...")
  # Process each cluster name separately (not using %in%) to report specific warning messages
  for(currentFile in FILE_FILTER_CELLS)
  {
    cat( paste0( "\n<br>File '", basename( currentFile), "'... "));
    
    fileContent = read.table( currentFile, 
                              sep = ",", 
                              header = TRUE, 
                              stringsAsFactors = FALSE)[["Cell"]]
    
    cellsSelection = colnames( sc10x) %in% fileContent;
    cat( paste0( sum( cellsSelection) , " cells."));
    
    if(! any( cellsSelection)) warning( paste0( "Argument FILE_FILTER_CELLS is specified but file '", basename( currentFile), "' did not identify any cell from Seurat object..."));
    cellsToRemove = cellsToRemove | cellsSelection;
  }
}


# Do the actual filtering
if(all( cellsToRemove)) 
{
  stop("Not enough cells remaining after filtering step...")
}

if(any( cellsToRemove))
{
  cat( paste0( "\n<br><br>Removing a total of ", sum( cellsToRemove), " cells (", sum( !cellsToRemove), " remaining)... "));
  sc10x = sc10x[, !cellsToRemove];
  
  # Eventually update identity factor levels
  Idents( sc10x) = factor( Idents( sc10x)) # Calling factor to recompute levels
  
  # And update loaded cells coordinates
  cellsCoordinates = cellsCoordinates[!cellsToRemove, ]
}




########################################
# CHECK MONITORED AND MODULE GENES LISTS 
########################################

### Remove eventual NULL (empty) list elements from list of genes to be monitored
monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];

# Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchMonitoredGenes = match( toupper( unlist( MONITORED_GENES)), toupper( rownames( GetAssayData( sc10x))));
matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( sc10x))));
monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( paste0("'", monitoredGenesNotFound, "'"), collapse=" - ")));
# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MONITORED_GENES = relist( rownames( GetAssayData( sc10x))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)
# Finally remove names that did not match
MONITORED_GENES = lapply( MONITORED_GENES, na.omit);



### Remove eventual NULL (empty) list elements from list of genes in modules
modulesGroupEmpty = sapply( MODULES_GENES, is.null);
if(any( modulesGroupEmpty)) warning( paste("Following module(s) of genes will be ignored because empty:", paste( names(modulesGroupEmpty)[modulesGroupEmpty], collapse=" - ")));
MODULES_GENES = MODULES_GENES[! modulesGroupEmpty];

# Check whether genes in MODULES_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchModulesGenes = match( toupper( unlist( MODULES_GENES)), toupper( rownames( GetAssayData( sc10x))));
matchModulesGenes = match( ( unlist( MODULES_GENES)), ( rownames( GetAssayData( sc10x))));
modulesGenesNotFound = unique( unlist( MODULES_GENES)[is.na( matchModulesGenes)]);
if(any( is.na( matchModulesGenes))) warning( paste( "Following gene(s) from modules list could not be found in experimental data and will be ignored:", paste( paste0("'", modulesGenesNotFound, "'"), collapse=" - ")));
# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MODULES_GENES = relist( rownames( GetAssayData( sc10x))[ matchModulesGenes ], skeleton = MODULES_GENES); # Does not work with NULL list elements (removed earlier)
# Finally remove names that did not match
MODULES_GENES = lapply( MODULES_GENES, na.omit);



### Transfer genes in very small modules (MODULES_GENES) to be analyzed individually (MONITORED_GENES)
modulesToTransfer = sapply(MODULES_GENES, length) < MONITORED_GENES_SMALL_MODULE
if(any(modulesToTransfer))
{
    warning( paste0( "Following Module(s) contained very few genes (<", MONITORED_GENES_SMALL_MODULE, "). These genes were transfered to 'Monitored genes' to be analyzed individually: ", paste( names(modulesToTransfer)[modulesToTransfer], collapse = " - ")));
    # Alter name so they can be recognized in list of Monitored genes
    names( MODULES_GENES)[modulesToTransfer] = paste0( "MOD_", names( MODULES_GENES)[modulesToTransfer]);
    MONITORED_GENES = c( MONITORED_GENES, MODULES_GENES[modulesToTransfer]);
    MODULES_GENES = MODULES_GENES[!modulesToTransfer];
}




#################
# Heat Shock data
#################

# Load the list of 512 genes associated with heat shock and stress in article from Campbell et al.
hs_stress_genes_df = read.table( file = PATH_HS_STRESS_MARKER_GENES_TABLE_FILE, header = TRUE, sep=",", quote = NULL)

# Convert list of genes from human symbol to mouse symbol
converted_genes_df = convertHumanGeneList( hs_stress_genes_df$gene_symbol)
hs_stress_genes_df = merge( hs_stress_genes_df, converted_genes_df, by.x = "gene_symbol", by.y="human_symbol")
hs_stress_genes_df = hs_stress_genes_df[ order( hs_stress_genes_df$PValue, hs_stress_genes_df$logFC, decreasing = c(FALSE,TRUE)), ]

DT::datatable( hs_stress_genes_df, caption = "List of HS and Stress associated genes in human")

MODULES_GENES[[ "HeatShock"]] = unique( hs_stress_genes_df[ , "symbol"])
MODULES_GENES[[ "HeatShock_Top40"]] = unique( hs_stress_genes_df[ 1:40, "symbol"])




###############
# SCATTER GENES
###############

# Check that all genes names for scatter plots comparisons are in Seurat object
validScatter = sapply( names(SCATTER_GENES), 
                  function(scatterName)
                  {
                    validGenes = SCATTER_GENES[[scatterName]] %in% rownames(sc10x)
                    if(all( validGenes)) return( TRUE)
                    warning( paste0( "Some genes selected for scatterplot comparisons '", scatterName, "' are not available in Seurat object: ", paste( SCATTER_GENES[[scatterName]][!validGenes], collapse = " - ")))
                    return( FALSE)
                  })
SCATTER_GENES = SCATTER_GENES[validScatter]




#################
# NORMALIZE DATA
#################

# ..............................................................................
## @knitr normalizeData
# ..............................................................................

if(RECOMPUTE_PREPROCESSING)
{
  cat("\n\n(Re-)Performing Normalization & Scaling...\n")
  sc10x = NormalizeData( object = sc10x,
                         normalization.method = DATA_NORM_METHOD,
                         scale.factor = DATA_NORM_SCALEFACTOR,
                         verbose = .VERBOSE)
  
  sc10x = ScaleData( object    = sc10x,
                     do.center = DATA_CENTER,
                     do.scale  = DATA_SCALE,
                     vars.to.regress = DATA_VARS_REGRESS,
                     verbose = .VERBOSE)
} else
{
  cat("\n\nUsing precomputed Normalization & Scaling...\n")
}
