###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#


ANALYSIS_STEP_NAME = "25_IgE_HighRes_updatedClusterNamesColors"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)


# Path to Seurat object (from previous analysis steps)
PATH_RDS_SEURAT_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, 
                                    "14c_ALL_NoT_NoVariableIG_NoHeavyChain_namedClusters",
                                    "AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_seuratObject_final.RDS")



# TSV file giving an hexadecimal-coded color for each cluster.
# Format: no header, cluster names as first column, colors in second column. All
# existing clusters must be defined in file (match by name, error otherwise). A
# warning is raised if more colors are declared than than existing clusters. 
# Empty string (or non-existing file path) for automatic coloring. 
# Ignored if RECOMPUTE_CLUSTERING is TRUE.
EXTERNAL_CLUSTERSCOLOR_PATH = file.path( PATH_EXPERIMENT_REFERENCE,
                                         "02_Analyses",
                                         "clustersColor.tsv")


# TSV file giving an hexadecimal-coded color for each HTO.
# Format: no header, HTOs names as first column, colors in second column. All
# existing HTOs must be defined in file (match by name, error otherwise). A
# warning is raised if more colors are declared than than existing HTOs. 
# Empty string (or non-existing file path) for automatic coloring. 
EXTERNAL_HTOSCOLOR_PATH = ""


# Dimreduc coordinates from Seurat object to be used for plot: "umap" "tsne" or 
# "pca". 
# Alternatively, a valid path to a TSV file containing cells coordinates from a
# previously computed 2-dimensions representation (dimensionality reduction). An
# error occurs for cells of Seurat object not found in file. A warning is raised
# for cells of file not found in Seurat object (cells ignored). As a consequence
# eventual subsetting MUST be made on Seurat object.
# Format: barcodes as rownames and two (named) columns for x and y coordinates.
CELLS_COORDINATES = "umap" # From loaded seurat object


# If 'CELLS_COORDINATES' is "pca", define which dimensions are used for 2D plots
PCA_DIMS = c( 1, 2);



#### Markers selection (for plotting only, results csv imported from 14c)
FINDMARKERS_PVAL_THR        = 0.001 # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP_TABLE   = 100   # Number of marker genes to show in report tables (NULL for all)
FINDMARKERS_SHOWTOP_HEATMAP = 5     # Number of marker genes to show in repot heatmaps (NULL for all)

#### HTO parameters

HTO_METADATA_COLUMN = "HTO_classification" # Name of the metadata column to use for HTO
HTO_FACTOR_LEVELS = NULL

#### Eventual clusters renaming and reordering for figures (by name, NULL to Ignore)

# Argument to function 'fct_recode'
CLUSTER_FACTOR_RECODE = list( "i-MBC" = "B mem 2",
                              "MBC" = "B mem 1",
                              "Act-MBC" = "Activated/GC",
                              "GC" = "GC",
                              "Pre-PC" = "Proliferating",
                              "PC1" = "Plasma PC3",
                              "PC2" = "Plasma PC1",
                              "PC3" = "Plasma PC2")


# Specify an ordering for cluster names (NULL to ignore)

ORDER_CLUSTER = c( "i-MBC", "MBC", "Act-MBC", "GC", "Pre-PC", "PC1", "PC2", "PC3")

#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (for Seurat using 'future')
NBCORES = 4;



