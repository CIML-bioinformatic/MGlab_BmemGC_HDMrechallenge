###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "14c_ALL_NoT_NoVariableIG_NoHeavyChain_namedClusters"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)




#Path to Seurat object (from previous analysis steps)
PATH_RDS_SEURAT_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, 
                                    "06c_ALL_GlobalHeterogeneity_NoT_NoVariableIG_NoHeavyChain",
                                    "ClusteringRes_0.4",
                                    "AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_seuratObject_final.RDS")



# TSV file listing cell identities to be used for clustering representation.
# Format: barcodes as rownames and one column named "identity" (other columns
# allowed but ignored).
# Barcodes not found in Seurat object are ignored (with a warning), missing ones
# are grouped as class 'unknown' (with warning).
# Invalid path (or empty string) to use eventual 'Idents' from Seurat object.
# Ignored if RECOMPUTE_CLUSTERING is TRUE.
EXTERNAL_CLUSTERING_PATH = file.path( PATH_EXPERIMENT_OUTPUT, 
                                    "13c_ALL_GroupClusters",
                                    "AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_cellsClusterIdentity.tsv")



# TSV file giving an hexadecimal-coded color for each cluster.
# Format: no header, cluster names as first column, colors in second column. All
# existing clusters must be defined in file (match by name, error otherwise). A
# warning is raised if more colors are declared than than existing clusters. 
# Empty string (or non-existing file path) for automatic coloring. 
# Ignored if RECOMPUTE_CLUSTERING is TRUE.
EXTERNAL_CLUSTERSCOLOR_PATH = ""



# TSV file giving an hexadecimal-coded color for each HTO.
# Format: no header, HTOs names as first column, colors in second column. All
# existing HTOs must be defined in file (match by name, error otherwise). A
# warning is raised if more colors are declared than than existing HTOs. 
# Empty string (or non-existing file path) for automatic coloring. 
EXTERNAL_HTOSCOLOR_PATH = ""



# Dimreduc coordinates from Seurat object to be used for plot: "umap" "tsne" or 
# "pca". If RECOMPUTE_PREPROCESSING is FALSE, the corresponding reduction must
# already be computed and stored in loaded Seurat object, an error is raised
# otherwise.
# Alternatively, a valid path to a TSV file containing cells coordinates from a
# previously computed 2-dimensions representation (dimensionality reduction). An
# error occurs for cells of Seurat object not found in file. A warning is raised
# for cells of file not found in Seurat object (cells ignored). As a consequence
# eventual subsetting MUST be made on Seurat object.
# Format: barcodes as rownames and two (named) columns for x and y coordinates.
# TSV overrides newly computed dim reduc when RECOMPUTE_PREPROCESSING is TRUE.
CELLS_COORDINATES = "umap" # Newly computed umap (see RECOMPUTE_PREPROCESSING)

# If 'CELLS_COORDINATES' is "pca", define which dimensions are used for 2D plots
PCA_DIMS = c( 1, 2);




#### Recomputing options

# Whether to recompute preprocessing steps (in case cells were filtered/Selected
# and one wants to use the new population only).
# Performs data Scaling/Centering, Normalization, detect variable genes, and
# compute PCAs and dimreduc with parameters defined here.
# Set to FALSE to use Seurat objects slots computed in previous steps. 
RECOMPUTE_PREPROCESSING = FALSE

# Whether to recompute a clustering. If TRUE, predefined identify and associated
# colors (EXTERNAL_CLUSTERING_PATH/EXTERNAL_CLUSTERSCOLOR_PATH) will be ignored.
RECOMPUTE_CLUSTERING = FALSE




#### HTO parameters

HTO_METADATA_COLUMN = "HTO_classification" # Name of the metadata column to use for HTO
HTO_FACTOR_LEVELS = NULL
#HTO_FACTOR_LEVELS = c("D21-PBS", "D21-Panth", "D29-PBS", "D29-Panth"); # Order of levels to use for HTOs (for plots, must match with values in selected metadata column). NULL to ignore.

#### Eventual clusters renaming and reordering for figures (by name, NULL to Ignore)

# Simplify clusters names (remove parenthesis information) #
SIMPLIFY_CLUSTER_NAMES = TRUE # One should take this in account in  'ORDER_CLUSTER'

# Specify an ordering for cluster names (NULL to ignore)
ORDER_CLUSTER = c( "B mem 1", "Activated/GC", "Plasma PC2", "GC", "Proliferating", "Plasma PC3", "Plasma PC1", "B mem 2")


#### Filtering variable genes

# Removing Immunoglobulin variable genes from set of HVF to prevent them from driving all analyses 
HVF_REMOVE_IG_VARIABLE = TRUE

# Eventually also remove IG heavy chain genes
HVF_REMOVE_IG_HEAVY = TRUE



#### Filtering cells (to be excluded before current and further analyses)

# Remove cells matching pattern in "factorHTO" (computed from HTO_METADATA_COLUMN)
FILTER_CELLS_HTO_REGEX = ".*Negative.*"

# Cluster(s) name to remove, NULL to ignore
FILTER_CLUSTERS = NULL;#"Dead";

# File(s) (csv) containing cells barcode to remove (in a column named 'Cell', as 
# in 'cellsExplorer' exports), NULL to ignore
FILE_FILTER_CELLS = NULL



# Path to the Cell cycle gene lists
CELL_CYCLE_SPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "02_Analyses", "01_CellCycle", "S_phase_genes.csv")))
CELL_CYCLE_G2MPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "02_Analyses","01_CellCycle", "G2M_phase_genes.csv")))

# Path to Heat Shock stress genes list
PATH_HS_STRESS_MARKER_GENES_TABLE_FILE = file.path( PATH_EXPERIMENT_REFERENCE, "02_Analyses", "02_HeatShock", "coregene_df-FALSE-v3.csv")

#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (for Seurat using 'future')
NBCORES = 4;


#### Plot options


# Scatterplot for group of cells (clusters/highlights) on dimreduc (ggplot)
PLOT_DIMREDUC_GROUPS_ALPHA     = 0.4;  # Alpha value (0-1) for points, 1 to mimic FeaturePlot and get maximum contrasts, less to prevent masking overlaid cells
PLOT_DIMREDUC_GROUPS_POINTSIZE = 1.5     # Size of points (0 for using Seurat:::AutoPointSize internal function)

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 20000;




#### Filtering / Normalization (in case reprocessing, see RECOMPUTE_PREPROCESSING)

# Filtering has eventually been done in previous steps

# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize";
DATA_NORM_SCALEFACTOR = 10000;

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE;
DATA_SCALE        = FALSE;
DATA_VARS_REGRESS = NULL;  # c("nCount_RNA") for UMIs (NULL to ignore)




#### Analysis parameters (in case reprocessing, see RECOMPUTE_PREPROCESSING)

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report

# Nearest-neighbor graph construction
FINDNEIGHBORS_K = 30

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION     = 0.8;
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
FINDCLUSTERS_ALGORITHM      = 1;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden


# Rename clustering results (only applies if RECOMPUTE_CLUSTERING is TRUE)
# RENAME_CLUSTERS = list()
# RENAME_CLUSTERS[[ SAMPLE_ALL]] = list()
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "0"]] = 0
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "1"]] = 5
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "2"]] = 1
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "3"]] = 3
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "4"]] = 6
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "5"]] = 2
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "6"]] = 4
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "7"]] = 7


# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results




#### Markers / Differential expression between (groups of) cells

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP_TABLE     = 100;       # Number of marker genes to show in report tables (NULL for all)
FINDMARKERS_SHOWTOP_HEATMAP   = 5;       # Number of marker genes to show in repot heatmaps (NULL for all)

# Parameter for enrichment analysis in GO Terms
ENRICHMENT_GO_PVALUECUTOFF = 0.05
ENRICHMENT_GO_QVALUECUTOFF = 0.05



#### List of group of clusters of interest

#CLUSTER_GROUP_LIST = list( T_NK_cell = c( 0, 3, 5, 7, 11, 14, 17, 19),
#                			     Neutrophil_cell = c( 12),
#                			     B_cell = c( 4,6,13,18),
#                			     Myeloid_cell = c( 1,2, 8,9, 10, 16, 20, 15))


#### Lists of genes of interest

# Number of genes under which genes of a module (MODULES_GENES) are transferred to be analyzed individually (MONITORED_GENES)
MONITORED_GENES_SMALL_MODULE = 0; # Disable
MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES = NULL

## Genes monitored individually (tsv file, one column for each group of genes)
#MONITORED_GENES = list()
GENES_SERGIO = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                                          "02_Analyses",
                                                          "03_MonitoredGenes",
                                                          "230925_Sergio.tsv"),
                                               sep = "\t",
                                               header = TRUE,
                                               stringsAsFactors = FALSE,
                                               row.names = NULL, 
                                               fill = TRUE));

# Read genes list from CUI sup table 4
GENES_CUI = read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                   "02_Analyses",
                                   "03_MonitoredGenes",
                                   "2024_11_12_41586_2023_6816_MOESM5_ESM_BCELLS.csv"),
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        row.names = NULL, 
                        fill = TRUE);

GENES_CUI[["Cytokine"]] = GENES_CUI[["Cytokine_Str"]]

CUI_pval_colname = "FDR"



# # Read genes list from CUI sup table 4
# GENES_CUI = read.table( file.path( PATH_EXPERIMENT_REFERENCE,
#                                    "02_Analyses",
#                                    "03_MonitoredGenes",
#                                    "2024_11_07_Cui_41586_2023_6816_MOESM6_ESM.tsv"),
#                         sep = "\t",
#                         header = TRUE,
#                         stringsAsFactors = FALSE,
#                         row.names = NULL, 
#                         fill = TRUE);
# 
# # Select categories of interest only
# GENES_CUI = GENES_CUI[GENES_CUI[["Celltype"]] == "B cell", ]
# 
#CUI_pval_colname = "P_val"

# Create R compatible names
GENES_CUI[["Cytokine"]] = make.names(GENES_CUI[["Cytokine"]])

# Split by cytokine
GENES_CUI_BY_CYTOKINE = lapply(split.data.frame(GENES_CUI, GENES_CUI[["Cytokine"]]), '[[', 'Gene')
names(GENES_CUI_BY_CYTOKINE) = paste0("CUI_All_", names(GENES_CUI_BY_CYTOKINE), "_", sapply(GENES_CUI_BY_CYTOKINE, length), "_")

# Get a top list of genes based on PValue
GENES_CUI_BY_CYTOKINE_TOP10_PVAL = by(GENES_CUI, GENES_CUI[["Cytokine"]], function(x){return(head(x[order(x[[CUI_pval_colname]]),][['Gene']], 10))})
names(GENES_CUI_BY_CYTOKINE_TOP10_PVAL) = paste0("CUI_Pval_Top10_", names(GENES_CUI_BY_CYTOKINE_TOP10_PVAL), "_", sapply(GENES_CUI_BY_CYTOKINE_TOP10_PVAL, length), "_")
GENES_CUI_BY_CYTOKINE_TOP50_PVAL = by(GENES_CUI, GENES_CUI[["Cytokine"]], function(x){return(head(x[order(x[[CUI_pval_colname]]),][['Gene']], 50))})
names(GENES_CUI_BY_CYTOKINE_TOP50_PVAL) = paste0("CUI_Pval_Top50_", names(GENES_CUI_BY_CYTOKINE_TOP50_PVAL), "_", sapply(GENES_CUI_BY_CYTOKINE_TOP50_PVAL, length), "_")

# Get a top list of genes based on Fold Change
GENES_CUI_BY_CYTOKINE_TOP10_FC = by(GENES_CUI, GENES_CUI[["Cytokine"]], function(x){return(head(x[order(x[["Avg_log2FC"]], decreasing = TRUE),][['Gene']], 10))})
names(GENES_CUI_BY_CYTOKINE_TOP10_FC) = paste0("CUI_FC_Top10_", names(GENES_CUI_BY_CYTOKINE_TOP10_FC), "_", sapply(GENES_CUI_BY_CYTOKINE_TOP10_FC, length), "_")
GENES_CUI_BY_CYTOKINE_TOP50_FC = by(GENES_CUI, GENES_CUI[["Cytokine"]], function(x){return(head(x[order(x[["Avg_log2FC"]], decreasing = TRUE),][['Gene']], 50))})
names(GENES_CUI_BY_CYTOKINE_TOP50_FC) = paste0("CUI_FC_Top50_", names(GENES_CUI_BY_CYTOKINE_TOP50_FC), "_", sapply(GENES_CUI_BY_CYTOKINE_TOP50_FC, length), "_")


# Read genes list from Ota et al. sup table 1
GENES_OTA = read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                   "02_Analyses",
                                   "03_MonitoredGenes",
                                   "2025_08_20_Ota_NIHMS1975647-supplement-Dataset_1.tsv"),
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        row.names = NULL, 
                        fill = TRUE);

# Get genes names from cluster5 only and convert to mouse IDs
GENES_OTA_CL5_MOUSE = babelgene::orthologs(GENES_OTA[GENES_OTA[["cluster"]]==5, "gene"], human = TRUE, species = "mouse")

GENES_OTA = list(OTA_CL5_TOMOUSE = GENES_OTA_CL5_MOUSE[["symbol"]])


# Additional genes for reviewing process
GENES_REVIEW = list("REVIEW_01" = c("Mki67", "Tfap4", "Foxo1", "Nr4a1"),
                    "REVIEW_02" = c("Fcgr1", "Mertk", "Adgre1", "Itgam"))


#names( GENES_HALLMARK_IFNA) = paste0( "HALLMARK_IFNA_", names(GENES_HALLMARK_IFNA)) # Single column with proper name
#names( GENES_LUETAL) = paste0( "LuEtAl_", names(GENES_LUETAL)) 

MONITORED_GENES = c( GENES_SERGIO, GENES_OTA, GENES_REVIEW)
# Remove empty strings
MONITORED_GENES = Map( '[', MONITORED_GENES, lapply(MONITORED_GENES, function(x){ which( nchar( x)>0)}));



## Genes monitored as modules (tsv file, one column for each group of genes)
#MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
#                                                "05_Modules",
#                                                "Modules.csv"),
#                                        sep = ",",
#                                        header = TRUE,
#                                        stringsAsFactors = FALSE,
#                                        row.names = NULL, fill = TRUE));
#MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings
MODULES_GENES = list()
MODULES_GENES[[ "S_PHASE"]] = CELL_CYCLE_SPHASE_GENELIST
MODULES_GENES[[ "G2M_PHASE"]] = CELL_CYCLE_G2MPHASE_GENELIST
# Also add following lists already added to monitored genes
MODULES_GENES = c(MODULES_GENES, 
                  GENES_SERGIO, 
                  GENES_OTA,
                  GENES_CUI_BY_CYTOKINE, 
                  GENES_CUI_BY_CYTOKINE_TOP10_PVAL, 
                  GENES_CUI_BY_CYTOKINE_TOP50_PVAL,
                  GENES_CUI_BY_CYTOKINE_TOP10_FC,
                  GENES_CUI_BY_CYTOKINE_TOP50_FC)

# Remove empty strings
MODULES_GENES = Map( '[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)}));


## Couple of genes to be compared on a 2D scatterplot

SCATTER_GENES = NULL; #list("CBC_Stem" = c("Lgr5", "Atoh1"))



