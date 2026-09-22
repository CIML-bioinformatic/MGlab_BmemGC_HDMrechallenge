###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#


ANALYSIS_STEP_NAME = "23_exportFiguresHighRes"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)


# Path to Seurat object (from previous analysis steps)
PATH_RDS_SEURAT_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, 
                                    "14c_ALL_NoT_NoVariableIG_NoHeavyChain_namedClusters",
                                    "AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_seuratObject_final.RDS")



# Markers list from previous analysis step (to be replotted as heatmap without IG genes, and updated clusters names)
PATH_MARKERS_LIST = file.path( PATH_EXPERIMENT_OUTPUT, 
                               "14c_ALL_NoT_NoVariableIG_NoHeavyChain_namedClusters",
                               "Extra",
                               "AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_MarkerGenes.tsv")

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


#### Plotting options

# Scatterplot for group of cells (clusters/highlights) on dimreduc (ggplot)
PLOT_DIMREDUC_GROUPS_ALPHA     = 0.4;  # Alpha value (0-1) for points, 1 to mimic FeaturePlot and get maximum contrasts, less to prevent masking overlaid cells
PLOT_DIMREDUC_GROUPS_POINTSIZE = 1.5     # Size of points (0 for using Seurat:::AutoPointSize internal function)

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 20000;



#### List of interesting signatures to subselect from CUI et al (from Seurat object) to add to newly loaded signatures (MODULE_GENES)

#SIGNATURES_NAMES = c("CUI_All_IL21_191_1", "CUI_All_CD40L_205_1", "CUI_All_IL4_722_1", "CUI_Pval_Top50_IL4_50_1", "CUI_FC_Top50_IL4_50_1")
# Replaced by direct subselection in CUI genes list loaded below, when creating MODULE_GENES list

#### Lists of genes of interest (for MODULES_GENES)

MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES = NULL


# Path to the Cell cycle gene lists
CELL_CYCLE_SPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "02_Analyses", "01_CellCycle", "S_phase_genes.csv")))
CELL_CYCLE_G2MPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "02_Analyses","01_CellCycle", "G2M_phase_genes.csv")))


# Heat Shock stress genes list
HEATSHOCK_GENELIST = list('heatshock_full' = readLines(file.path(PATH_EXPERIMENT_REFERENCE, "02_Analyses", "02_HeatShock", "MOUSE_heatshock_full.txt")),
                          'heatshock_top40' = readLines(file.path(PATH_EXPERIMENT_REFERENCE, "02_Analyses", "02_HeatShock", "MOUSE_heatshock_top40.txt")));

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

# Read list of genes from Laurine Binet
GENES_LAURINE_DF = list()

# Get as a list of genes with associated Pvalue and fold change
GENES_LAURINE_DF[["PrePC_Act"]] = read.csv( file.path( PATH_EXPERIMENT_REFERENCE,
                                                       "02_Analyses",
                                                       "03_MonitoredGenes",
                                                       "Laurine",
                                                       "PrePC.Act.csv"),
                                            sep = ";",
                                            header = TRUE,
                                            stringsAsFactors = FALSE,
                                            row.names = NULL, 
                                            fill = TRUE);

# Get genes as DF from DEG result
GENES_LAURINE_DF[["PrePC_ERStressLow"]] = read.csv( file.path( PATH_EXPERIMENT_REFERENCE,
                                                               "02_Analyses",
                                                               "03_MonitoredGenes",
                                                               "Laurine",
                                                               "PrePC.ERStresslo.csv"),
                                                    sep = ";",
                                                    header = TRUE,
                                                    stringsAsFactors = FALSE,
                                                    row.names = NULL, 
                                                    fill = TRUE);

# Add an eventual label to top genes by PValue in each category (positive and negative enrichment)
GENES_LAURINE_DF[["PrePC_ERStressLow"]] = do.call(rbind, by(GENES_LAURINE_DF[["PrePC_ERStressLow"]], GENES_LAURINE_DF[["PrePC_ERStressLow"]][["avg_log2FC"]]>0, function(x){x[["label"]] = ifelse(rank(x[["p_val"]])<15, x[["Genes"]], ""); return(x)}))
GENES_LAURINE_DF[["PrePC_Act"]] = do.call(rbind, by(GENES_LAURINE_DF[["PrePC_Act"]], GENES_LAURINE_DF[["PrePC_Act"]][["avg_log2FC"]]>0, function(x){x[["label"]] = ifelse(rank(x[["p_val"]])<15, x[["Genes"]], ""); return(x)}))


# Show a volcano plot of DEG result
#ggplot(GENES_LAURINE_DF[["PrePC_ERStressLow"]], aes(x=avg_log2FC, y=-log10(p_val))) + geom_point() + ggrepel::geom_label_repel(aes(label = label))
#ggplot(GENES_LAURINE_DF[["PrePC_Act"]], aes(x=avg_log2FC, y=-log10(p_val))) + geom_point() + ggrepel::geom_label_repel(aes(label = label))

# Select positive only from list



# Filter based on adjusted PValue threshold and positive only "PrePC_ERStressLow"
GENES_LAURINE_DF = lapply(GENES_LAURINE_DF, function(x)
{
  return(x[x[["p_val_adj"]]<0.05 & x[["avg_log2FC"]]>0, ] )
})

# Get genes names only (based on various classification)
GENES_LAURINE_TOP10_PVAL = lapply(GENES_LAURINE_DF, function(x)
{
  head(x[ order( x[["p_val_adj"]]), "Genes"], 10)
})

GENES_LAURINE_TOP50_PVAL = lapply(GENES_LAURINE_DF, function(x)
{
  head(x[ order( x[["p_val_adj"]]), "Genes"], 50)
})

GENES_LAURINE_TOP10_FC = lapply(GENES_LAURINE_DF, function(x)
{
  head(x[ order( x[["avg_log2FC"]], decreasing = TRUE), "Genes"], 10)
})

GENES_LAURINE_TOP50_FC = lapply(GENES_LAURINE_DF, function(x)
{
  head(x[ order( x[["avg_log2FC"]], decreasing = TRUE), "Genes"], 50)
})




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

# 
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
# CUI_pval_colname = "P_val"

# Create R compatible names
GENES_CUI[["Cytokine"]] = make.names(GENES_CUI[["Cytokine"]])

# Split by cytokine
GENES_CUI_BY_CYTOKINE = lapply(split.data.frame(GENES_CUI, GENES_CUI[["Cytokine"]]), '[[', 'Gene')
#names(GENES_CUI_BY_CYTOKINE) = paste0("CUI_All_", names(GENES_CUI_BY_CYTOKINE), "_", sapply(GENES_CUI_BY_CYTOKINE, length), "_")

# Get a top list of genes based on PValue
GENES_CUI_BY_CYTOKINE_TOP10_PVAL = by(GENES_CUI, GENES_CUI[["Cytokine"]], function(x){return(head(x[order(x[[CUI_pval_colname]]),][['Gene']], 10))})
#names(GENES_CUI_BY_CYTOKINE_TOP10_PVAL) = paste0("CUI_Pval_Top10_", names(GENES_CUI_BY_CYTOKINE_TOP10_PVAL), "_", sapply(GENES_CUI_BY_CYTOKINE_TOP10_PVAL, length), "_")
GENES_CUI_BY_CYTOKINE_TOP50_PVAL = by(GENES_CUI, GENES_CUI[["Cytokine"]], function(x){return(head(x[order(x[[CUI_pval_colname]]),][['Gene']], 50))})
#names(GENES_CUI_BY_CYTOKINE_TOP50_PVAL) = paste0("CUI_Pval_Top50_", names(GENES_CUI_BY_CYTOKINE_TOP50_PVAL), "_", sapply(GENES_CUI_BY_CYTOKINE_TOP50_PVAL, length), "_")

# Get a top list of genes based on Fold Change
GENES_CUI_BY_CYTOKINE_TOP10_FC = by(GENES_CUI, GENES_CUI[["Cytokine"]], function(x){return(head(x[order(x[["Avg_log2FC"]], decreasing = TRUE),][['Gene']], 10))})
#names(GENES_CUI_BY_CYTOKINE_TOP10_FC) = paste0("CUI_FC_Top10_", names(GENES_CUI_BY_CYTOKINE_TOP10_FC), "_", sapply(GENES_CUI_BY_CYTOKINE_TOP10_FC, length), "_")
GENES_CUI_BY_CYTOKINE_TOP50_FC = by(GENES_CUI, GENES_CUI[["Cytokine"]], function(x){return(head(x[order(x[["Avg_log2FC"]], decreasing = TRUE),][['Gene']], 50))})
#names(GENES_CUI_BY_CYTOKINE_TOP50_FC) = paste0("CUI_FC_Top50_", names(GENES_CUI_BY_CYTOKINE_TOP50_FC), "_", sapply(GENES_CUI_BY_CYTOKINE_TOP50_FC, length), "_")


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
GENES_REVIEW = list( REVIEW_01 = c( "Mki67", "Tfap4", "Foxo1", "Nr4a1"),
                     REVIEW_02 = c( "Fcgr1", "Mertk", "Adgre1", "Itgam"),
                     REVIEW_03 = c( "Nt5e"))

# Additional genes for scatterplot analysis
SCATTER_GENES = list( CD80_PDL2 = c( "Cd80", "Pdcd1lg2"))


#### SIGNATURES analysis

# Same approach as in other scripts (modulescore) but also compares conditions
# Warning : as opposed to other scripts, here MODULES_GENES is a 'list of list'
# to allow grouping of signatures

MODULES_GENES = list()

MODULES_GENES = list( "CYCLE" = list("S_PHASE" = CELL_CYCLE_SPHASE_GENELIST, 
                                     "G2M_PHASE" = CELL_CYCLE_G2MPHASE_GENELIST),
                      "SERGIO" = GENES_SERGIO,
                      # "OTA_CL5" = GENES_OTA,
                      # "CUI_SELECTION" = c( GENES_CUI_BY_CYTOKINE[c("IL21", "CD40L", "IL4")],
                      #                      list("IL4_top50_FC" = GENES_CUI_BY_CYTOKINE_TOP50_FC[["IL4"]]),
                      #                      list("IL4_top50_PV" = GENES_CUI_BY_CYTOKINE_TOP50_PVAL[["IL4"]])),
                      "GENES_REVIEW" = GENES_REVIEW,
                      "HEATSHOCK" = HEATSHOCK_GENELIST,
                      "SCATTER_GENES" = SCATTER_GENES)
                      # "CUI_ALL" = GENES_CUI_BY_CYTOKINE,
                      # "CUI_T10_PV" = GENES_CUI_BY_CYTOKINE_TOP10_PVAL,
                      # "CUI_T50_PV" = GENES_CUI_BY_CYTOKINE_TOP50_PVAL,
                      # "CUI_T10_FC" = GENES_CUI_BY_CYTOKINE_TOP10_FC,
                      # "CUI_T50_FC"= GENES_CUI_BY_CYTOKINE_TOP50_FC,
                      # "LAU_T10_PV" = GENES_LAURINE_TOP10_PVAL,
                      # "LAU_T50_PV" = GENES_LAURINE_TOP50_PVAL,
                      # "LAU_T10_FC" = GENES_LAURINE_TOP10_FC,
                      # "LAU_T50_FC" = GENES_LAURINE_TOP50_FC)
                      

# Remove empty strings
for(currentModule in names(MODULES_GENES)) MODULES_GENES[[currentModule]] = Map( '[', MODULES_GENES[[currentModule]], lapply(MODULES_GENES[[currentModule]], function(x){ which( nchar( x)>0)}));


#### List of features of interest (featureplot)

#MONITORED_GENES = list("SERGIO" = c("Cd19", "Cd38", "Cd83", "Ccr6", "Ccr7", "Aicda", "Bcl6", "S1pr2", "Cd86", "Fas", "Mki67", "Jchain", "Prdm1", "Sdc1", "Xbp1", "Il13ra1", "Fcer2a", "Il4ra"))

# Create a list of monitored genes from lists of MODULES_GENES (concatenate family and group name)
MONITORED_GENES_FROM_MODULES = list();
# Compute the score of the cells according to group of monitored genes
for(currentModule in names(MODULES_GENES))
{
  for( listName in names( MODULES_GENES[[currentModule]]))
  {
    #message(paste(currentModule, listName, sep="_"));
    if( length( MODULES_GENES[[currentModule]][[listName]]) == 0)
    {
      warning( paste0( "List of genes in module '", listName, "' is empty, ignoring..."));
    } else
    {
      MONITORED_GENES_FROM_MODULES[[paste(currentModule, listName, sep="_")]] = MODULES_GENES[[currentModule]][[listName]]
    }
  }
}

MONITORED_GENES = c(MONITORED_GENES_FROM_MODULES)

# Here we arbitrarily cap the max number of MONITOREDGENES in each list to prevent signatures overflowing the report with individual featureplots
MONITORED_MAX = 50
maxedCategories = sapply(MONITORED_GENES, length) > MONITORED_MAX
if(any(maxedCategories)) 
  {
    warning(paste("At least one category in MONITORED_GENES was arbitrarily capped to a maximum of", MONITORED_MAX, "elements to prevent report bloating. Consider checking your list of genes."))
    MONITORED_GENES = lapply(MONITORED_GENES, head, MONITORED_MAX)
  }

