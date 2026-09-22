###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "09d_RNAvelocity_from06b_LG_NoT_NoVariableIG_NoHeavyChain"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "RNA velocity conversion"

# Path to the loom files (computed in earlier step by 'velocyto' command line tool)
PATH_LOOM = dir( file.path( PATH_EXPERIMENT_OUTPUT, "08_velocyto"),
                 pattern = "*.loom",
                 full.names = TRUE);

# Base path to Seurat object (from previous analysis steps)
# Modified from other scripts we resolve the full path in prepareData script
# so we can use FINDCLUSTERS_RESOLUTION to get the proper one.
BASE_PATH_RDS_SEURAT_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, 
                                         "06b_LG_GlobalHeterogeneity_NoT_NoVariableIG_NoHeavyChain")


# WARNING: This is a bypass in case loom and cellranger annotations don't match.
# Annotations not matching are REMOVED, NOT RECOMMENDED for normal use.
FORCE_MERGE_LOOM_SEURAT = TRUE; # FALSE recommended (error if features don't match)


# TSV file listing cell identities to be used for clustering representation.
# Format: barcodes as rownames and one column named "identity" (other columns
# allowed but ignored).
# Barcodes not found in Seurat object are ignored (with a warning), missing ones
# are grouped as class 'unknown' (with warning).
# Invalid path (or empty string) to use eventual 'Idents' from Seurat object.
EXTERNAL_CLUSTERING_PATH = "" 


# TSV file giving an hexadecimal-coded color for each cluster.
# Format: no header, cluster names as first column, colors in second column. All
# existing clusters must be defined in file (match by name, error otherwise). A
# warning is raised if more colors are declared than than existing clusters. 
# Empty string (or non-existing file path) for automatic coloring. 
EXTERNAL_CLUSTERSCOLOR_PATH = ""


# Dimreduc coordinates from Seurat object to be used for plot: "umap" "tsne" or 
# "pca". Must already be computed and stored in loaded Seurat object, an error 
# is raised otherwise.
# Alternatively, a valid path to a TSV file containing cells coordinates from a
# previously computed 2-dimensions representation (dimensionality reduction). An
# error occurs for cells of Seurat object not found in file. A warning is raised
# for cells of file not found in Seurat object (cells ignored). As a consequence
# eventual subsetting MUST be made on Seurat object.
# Format: barcodes as rownames and two (named) columns for x and y coordinates.
CELLS_COORDINATES = "umap"

# If 'CELLS_COORDINATES' is "pca", define which dimensions are used for 2D plots
PCA_DIMS = c( 1, 2);




#### Eventual root cells and end points for LATENT TIME

# A numeric ID defining eventual root cell and end points to be used fot latent-time estimation 
# ('integer(0)' to ignore, makes an empty list in python, more convenient to handle than R 'NULL').
ROOT_CELLS_ID = integer(0)
END_POINTS_ID = integer(0)

# For this sample, we load the list of end points from a mnually identified set of cells (from velocities)
#END_POINTS_ID = as.numeric( readLines( file.path( PATH_PROJECT_EXTERNALDATA, 
#                                      "05_rootCells_endPoints_latentTime", 
#                                      "2022_01_24_endpoints_j_MO_onlyLysoMAC.txt")));




#### Filtering cells (to be excluded before current and further analyses)

# Cluster(s) name to remove, NULL to ignore
FILTER_CLUSTERS = NULL;

# File(s) (csv) containing cells barcode to remove (in a column named 'Cell', as 
# in 'cellsExplorer' exports), NULL to ignore
FILE_FILTER_CELLS = NULL;




#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (using 'future' for Seurat3, mclapply for
# other loops)
NBCORES = 4;

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 10000;




#### Loading, normalization and analysis parameters 

# filter_and_normalize (Filtering, normalization and log transform)
FILTER_NORMALIZE_MIN_SHARED_COUNTS = 20;   # Minimum number of counts (both unspliced and spliced) required for a gene.
FILTER_NORMALIZE_MIN_SHARED_CELLS  = 20;   # Minimum number of cells required to be expressed (both unspliced and spliced).
FILTER_NORMALIZE_N_TOP_GENES       = 5000; # Number of genes to keep.



# moments (Computes moments for velocity estimation)
MOMENTS_N_PCS = 30;       # Number of principal components to use. If not specified, the full space is used of a pre-computed PCA, or 30 components are used when PCA is computed internally.
MOMENTS_N_NEIGHBORS = 10; # Number of neighbors to use (default: 30)

# recover_dynamics (Recovers the full splicing kinetics of specified genes, i.e. fit to model)
RECOVER_DYNAMICS_MAX_ITER = 25; # Maximal iterations in the EM-Algorithm for fitting (default: 10)



# Number of genes to show in reports
SCVELO_LIKELIHOOD_TOP  = 150; # Also define number shown on heatmap and used for differential kinetics
SCVELO_DRIVERGENES_TOP = 50;  # (Re-)Identification of cluster driver genes by scVelo added information (differential velocity and partial likelihood)
SCVELO_DRIVERGENES_FUNCTIONAL_TOP = SCVELO_DRIVERGENES_TOP; # Number of top driver genes to use for functional analyses (KEGG/GO)




#### Plot options

## Select base represntation for most plots based on dimensional reduction
PLOT_DIMREDUC_BASIS = "dimreduc"; # Dimreduc object from Seurat object (umap/tsne). Object in slot "dimreduc" is created with selected coordinates (see CELLS_COORDINATES) before conversion of Seurat object for python/scVelo.
PLOT_WIDTH  = 8;
PLOT_HEIGHT = 8;
PLOT_DPI = 120;

## Plots showing velocities on selected embedding
PLOT_VELOCITY_EMBEDDING_BASIS       = PLOT_DIMREDUC_BASIS;
# Individual cells
PLOT_VELOCITY_EMBEDDING_ARROWSIZE   = 2.5;
PLOT_VELOCITY_EMBEDDING_ARROWLENGTH = 5.0;
# Streams
PLOT_VELOCITY_EMBEDDING_STREAM_DENSITY     = 2;      # Amount of velocities to show (float 0=none to default 1=all as per documentation, can be higher), linewidth adjusted accordingly (1/value)
PLOT_VELOCITY_EMBEDDING_STREAM_SMOOTH      = 0.2;    # Scale factor for gaussian kernel around grid point (default 0.5)
PLOT_VELOCITY_EMBEDDING_STREAM_MINMASS     = 0.2;    # Minimum threshold for mass to be shown (0=all velocities to 5=large velocities only)
PLOT_VELOCITY_EMBEDDING_STREAM_LINEWIDTH   = 1.5/PLOT_VELOCITY_EMBEDDING_STREAM_DENSITY;
PLOT_VELOCITY_EMBEDDING_STREAM_POINT_SIZE  = 20;
PLOT_VELOCITY_EMBEDDING_STREAM_POINT_ALPHA = 0.3;
# Grid
PLOT_VELOCITY_EMBEDDING_GRID_DENSITY     = 2;      # Density of velocities stream lines (float 0=none to default 1=all as per documentation, can be higher)
PLOT_VELOCITY_EMBEDDING_GRID_SMOOTH      = 0.2;    # Scale factor for gaussian kernel around grid point (default 0.5)
PLOT_VELOCITY_EMBEDDING_GRID_MINMASS     = 0.2;    # Minimum threshold for mass to be shown (0=all velocities to 100=large velocities)
PLOT_VELOCITY_EMBEDDING_GRID_SCALE       = 1.0;    # Length of velocities in the embedding
PLOT_VELOCITY_EMBEDDING_GRID_ARROWSIZE   = 2.0;
PLOT_VELOCITY_EMBEDDING_GRID_ARROWLENGTH = 2.0;
PLOT_VELOCITY_EMBEDDING_GRID_POINT_SIZE  = 20;
PLOT_VELOCITY_EMBEDDING_GRID_POINT_ALPHA = 0.3;


## Velocities graph and connectivity
PLOT_VELOCITY_GRAPH_NN        = 10;   # Number of neighbors to be included
PLOT_VELOCITY_GRAPH_ARROWSIZE = 1.0;  # Size of arrow heads (0 to disable arrows)

## PAGA (based on previously computed seurat clusters)
PLOT_PAGA_BASIS             = PLOT_DIMREDUC_BASIS;
PLOT_PAGA_MIN_SPANNING_TREE = TRUE;   # Path from A-to-B is removed if a more confident one exists
PLOT_PAGA_USE_TIME_PRIOR    = TRUE;   # Obs key for pseudo-time values. If TRUE, ‘velocity_pseudotime’ is used if available.
PLOT_PAGA_THRESHOLD         = 0.01;   # Do not draw edges for weights below this threshold. 0 for all edges.
PLOT_PAGA_LAYOUT            = "fr";   # Plotting layout that computes positions (igraph: fr/fa/rt/eq_tree/...).
PLOT_PAGA_NODE_SCALE        = 1;
PLOT_PAGA_NODE_ARROWSIZE    = 30;

## Heatmap of top likelihood genes
PLOT_HEATMAP_NCONVOLVE = 30;    # Size of averaging window along 'latent-time' (default 30)
PLOT_HEATMAP_COLORBAR  = TRUE;  # Whether to show colorbar.
PLOT_HEATMAP_WIDTH     = 10;
PLOT_HEATMAP_HEIGHT    = 20;



#### Genes of interest (for individual velocity/phase/time plots)

MONITORED_GENES = list();




#### Functional enrichment analyses

FUNCTIONAL_ORGANISM_LIBRARY = "org.Mm.eg.db" # "org.Hs.eg.db" for Homo Sapiens. Bioconductor annotation for current species (convert gene symbol to entrez id)
FUNCTIONAL_ORGANISM_CODE = "mmu" # "hsa" for Homo Sapiens (see http://www.genome.jp/kegg/catalog/org_list.html) 

## KEGG functionnal analysis (clusterprofiler)
KEGG_PADJUST_METHOD   = "BH"
KEGG_PVALUE_CUTOFF    = 0.05
KEGG_TOPTERMS_SUMMARY = 5    # Number of top terms (based on 'p.adjust') to show for each 'sample' (summary figures with genes names)
KEGG_TOPTERMS_FIGURE  = 20   # Number of top terms (based on 'p.adjust') to show for each 'sample' (individual sample/cluster figures)

## GO functionnal analysis (clusterprofiler)
GO_CATEGORIES         = list( "Biological Process"   = "BP", # GO categories to check enrichments for (remove unused elements).
                              "Cellular Compartment" = "CC", # Elements names shown on the report (can be changed).
                              "Molecular Function"   = "MF") # Elements values used for enrichment function (do not change).
GO_PADJUST_METHOD     = "BH"
GO_PVALUE_CUTOFF      = 0.05
GO_SIMPLIFY_CUTOFF    = 0.7 # Cutoff for eventual simplification of redundant GO terms in results (NULL to ignore simplification)
GO_TOPTERMS_SUMMARY   = 5
GO_TOPTERMS_FIGURE    = 20

#KEGG_UNIVERSE_IS_UNIONGENES = TRUE # Controlled by Rmd file (both values). Define if the background for enrichment analysis must be defined as the union of genes from all clusters (all genes otherwise)
#GO_UNIVERSE_IS_UNIONGENES = TRUE # Controlled by Rmd file (both values). Define if the background for enrichment analysis must be defined as the union of genes from all clusters (all genes otherwise)


