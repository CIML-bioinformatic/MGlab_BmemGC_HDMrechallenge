# ##########################################
# This script reads a precomputed Seurat object, and eventual external data from
# TSV files for clustering results and 2D representation of cells such as result
# of dimensionality reduction algorithms. 
#
# ##########################################



## @knitr loadDataTest
# 
# binaryResultFilename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "binarySeuratWithLoom.h5Seurat"));
# 
# # Increase max mem size for future (parallel) to 2048MB
# options(future.globals.maxSize= 2048*1024^2)
# 
# 
# ldat <- ReadVelocity(file = PATH_LOOM)
# bm <- as.Seurat(x = ldat)
# bm[["RNA"]] <- bm[["spliced"]]
# bm <- SCTransform(bm)
# bm <- RunPCA(bm)
# bm <- RunUMAP(bm, dims = 1:20)
# bm <- FindNeighbors(bm, dims = 1:20)
# bm <- FindClusters(bm)
# DefaultAssay(bm) <- "RNA"
# SaveH5Seurat(bm, filename = binaryResultFilename)
# Convert(binaryResultFilename, dest = "h5ad")



## @knitr loadData

#### Seurat object

# Resolve path to appropriate RDS file for Seurat object (with ClusteringRes)
PATH_RDS_SEURAT_OBJECT = dir( file.path( BASE_PATH_RDS_SEURAT_OBJECT, 
                                         paste0("ClusteringRes_", FINDCLUSTERS_RESOLUTION)), # Use variable defined in environment by launchReportCompilation
                              pattern = ".*_seuratObject_final.RDS",
                              full.names = TRUE);

# Load Seurat from previously saved binary RDS file
sc10x = readRDS(PATH_RDS_SEURAT_OBJECT);

cat( paste0( "\n<br>Successfuly loaded Seurat object: ", ncol( sc10x), " Cells x ", nrow( sc10x)," Genes.\n"));



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

# Save current identity (as metadata = stash idents) so it is included during the conversion for scVelo
sc10x[["identity"]] = Idents(object = sc10x) # Metadata can be loaded back into current identity with Idents(sc10X) = 'identity'




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

# Convert named vector to a named list for conserving names during onversion to python by reticulate (named list -> dict)
clustersColorList = as.list(clustersColor)




# Read eventual TSV file to define cells coordinates
cellsCoordinates = NULL;
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
  
} else if( CELLS_COORDINATES %in% Reductions( sc10x) ) 
{
  cellsCoordinates = Embeddings(Reductions(sc10x, slot = CELLS_COORDINATES))
  
  # PCA can have arbitrary number of dimensions, use requested ones (2d plot)
  if(CELLS_COORDINATES == "pca")
  {
    cellsCoordinates = cellsCoordinates[, PCA_DIMS[1:2]];
  }
  
} else stop("Parameter 'CELLS_COORDINATES' must refer to a DimReduc stored in loaded Seurat object, or contain a path to a valid file containing 2D cells coordinates...");

# Add these coordinates as additional DimReduc into Seurat object so it is included during the conversion for scVelo
if(!is.null(cellsCoordinates))
{
  sc10x[['dimreduc']] = CreateDimReducObject(embeddings = as.matrix(cellsCoordinates), assay = Seurat::DefaultAssay(sc10x), key = "DIMRED_", global = TRUE)
}




#### Eventual filtering of cells to exclude from further analyses (the actual point of this report)

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
}




#### Loom data

# Prepare a base filename for conversion of Seurat object to binary file
binaryBaseFilename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "binarySeuratWithLoom"));


# Load loom files (with spliced/unspliced counts) computed by 'velocyto' tool
loomDataList = lapply(PATH_LOOM, ReadVelocity);
#seuratLoom = as.Seurat( loomData)

# Reformat velocyto barcodes to match Seurat ones in each matrix
loomDataList = lapply( loomDataList, lapply, function(x)
            {
              colnames( x) = gsub( "x$", "-1", gsub( ".*:", "", colnames( x))) # Replace ending 'x' by '-1', and remove ':' and everything before 
              return( x);
            })

# Merge each matrix from all objects together (MUST HAVE SAME OOWS/ANNOTATIONS)
loomData = Reduce(function(x, y)
                  {
                    Map(cbind, x, y)
                  },
                  loomDataList)
# Generate a separate binary file for each condition/HTO to be used in scvelo
#for(currentCondition in levels(sc10x[["factorHTO", drop = TRUE]]))
{
  # Subset Seurat object with cells selected for current condition/HTO
  #currentSC = sc10x[,sc10x[["factorHTO", drop = TRUE]] %in% currentCondition]
  #binaryResultFilename = paste0(binaryBaseFilename, "_HTO_", currentCondition, ".h5Seurat")
  currentSC = sc10x
  binaryResultFilename = paste0(binaryBaseFilename, ".h5Seurat")
  
  # Get index of barcodes from Seurat object (colnames) in Loom matrix (colnames)
  barcodesMatchInd = match( colnames( currentSC), colnames( loomData[[1]]));
  barcodesNotFound = is.na ( barcodesMatchInd);
  if( any( barcodesNotFound)) stop( paste0( sum( barcodesNotFound), " barcode(s) from Seurat object could not be found in Loom barcodes, please check consistency: ", paste( colnames( sc10x)[barcodesNotFound], collapse = " - "), "."));

  # Get index of features from Seurat object (rownames) in Loom matrix (rownames)
  featuresMatchInd = match( rownames( currentSC), rownames( loomData[[1]]));
  featuresNotFound = is.na ( featuresMatchInd);
  if( any( featuresNotFound)) 
  {
    # prepare a warning/error message
    errorMessage = paste0( sum( featuresNotFound), " feature(s) from Seurat object could not be found in Loom features, please check consistency: ", paste( rownames( sc10x)[featuresNotFound], collapse = " - "), ".");
    
    if(FORCE_MERGE_LOOM_SEURAT) # Force and remove missing genes from Seurat object ?
    {
      # Filter features from Seurat object (not recommended, see analysisParams.R)
      currentSC = currentSC[!featuresNotFound, ];
      warning( paste( "Some genes were removed from analysis because", errorMessage));
    }
    else stop( errorMessage);
  }
  
  #all(rownames(loomData[["spliced"]])==rownames(currentSC))
  
  # Add unspliced/spliced/ambiguous matrices into existing Seurat object
  currentSC[["spliced"]]   = CreateAssayObject( loomData[["spliced"]]  [rownames( currentSC), colnames( currentSC)]);
  currentSC[["unspliced"]] = CreateAssayObject( loomData[["unspliced"]][rownames( currentSC), colnames( currentSC)]);
  currentSC[["ambiguous"]] = CreateAssayObject( loomData[["ambiguous"]][rownames( currentSC), colnames( currentSC)]);

  # Set spliced counts as RNA abundance measure
  currentSC[["RNA"]] = currentSC[["spliced"]] # Replace existing RNA assay beacause 'scale.data' shape is altered (variable genes), and it will prevent inclusion of other layers (spliced/unspliced...) into scVelo file (see ?Convert).
  DefaultAssay(currentSC) = "RNA" # Make sure RNA is the selected default assay

  # Convert factor metadata to actual string values
  factorMetadata = colnames(currentSC[[]])[sapply(currentSC[[]], is.factor)] # Get the columns name of metadata stored as factors
  if(length(factorMetadata))
  {
    currentSC[[factorMetadata]] = lapply(currentSC[[factorMetadata]], as.character)
  }

  # Create a temporary file in local filesystem as Writing sometimes fails on remote filesystem mounted
  tempFilename = tempfile(paste0( "binarySeuratWithLoom_"), fileext = ".h5Seurat")
  
  

  # Save Seurat  object as 'H5' binary file and convert it for reading with Python
  SaveH5Seurat( currentSC, filename = tempFilename, overwrite = TRUE);
  Convert( tempFilename, dest = "h5ad", overwrite = TRUE);

  # Then copy both files from temp to result folder
  stopifnot( file.copy( tempFilename, binaryResultFilename))
  stopifnot( file.copy(gsub("h5Seurat", "h5ad", tempFilename), gsub("h5Seurat", "h5ad", binaryResultFilename)))
  
  
}
