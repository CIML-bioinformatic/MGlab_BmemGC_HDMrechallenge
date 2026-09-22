###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#

# CURRENT_MOUSE = "M4" Defined dynamically in launching script "launch_reports_compilation.R" loop to generate one report for each mouse

ANALYSIS_STEP_NAME = "24_clonotypeOverlap_scRepertoire_updatedClusterNamesColors"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)





#Path to Seurat object (from previous analysis steps)
PATH_RDS_SEURAT_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, 
                                    "14c_ALL_NoT_NoVariableIG_NoHeavyChain_namedClusters",
                                    "AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_seuratObject_final.RDS")

# For scRepertoire
PATH_VDJ_FILE_FILTERED  = file.path( PATH_EXPERIMENT_OUTPUT,
                                     "01_CellRanger_FeatureBarcoding",
                                     "mm10",
                                     "outs",
                                     "per_sample_outs",
                                     "mm10",
                                     "vdj_b",
                                     "filtered_contig_annotations.csv")




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


