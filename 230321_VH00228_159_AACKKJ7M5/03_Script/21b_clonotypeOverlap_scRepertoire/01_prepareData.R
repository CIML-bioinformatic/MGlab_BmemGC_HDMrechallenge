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

# Load Seurat from previously saved binary RDS file
sc10x = readRDS(PATH_RDS_SEURAT_OBJECT);

# Report Idents to a dedicated metadata slot 
sc10x[["Identity"]] = Idents(sc10x)

# Quick fix for the Seurat object created in previous version of Seurat without project name
sc10x@project.name = "Project"


cat( paste0( "\nSuccessfuly loaded Seurat object: ", ncol( sc10x), " Cells x ", nrow( sc10x)," Genes.\n"));

# Get gene names of IGs based on base (to be compared with Ighe)
#IgGenesSelection = unlist(lapply(IgGenesSelectionBase, grep, rownames(sc10x), value = TRUE))
#names(IgGenesSelection)=IgGenesSelection

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


# SUBSET current object to selected mouse
cat("Subsetting Seurat object for selected sample :", CURRENT_MOUSE)
sc10x = subset(sc10x, subset = HTOMICE_maxID == CURRENT_MOUSE)

print(sc10x)
