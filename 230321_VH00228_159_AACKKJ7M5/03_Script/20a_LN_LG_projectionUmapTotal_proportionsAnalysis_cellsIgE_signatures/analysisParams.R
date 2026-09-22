###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "20a_LN_LG_projectionUmapTotal_proportionsAnalysis_cellsIgE_signatures"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)




#Path to Seurat object (from previous analysis steps)
PATH_RDS_SEURAT_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, 
                                    "14c_ALL_NoT_NoVariableIG_NoHeavyChain_namedClusters",
                                    "AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_seuratObject_final.RDS")


PATH_VDJ_FILE = file.path( PATH_EXPERIMENT_OUTPUT,
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
EXTERNAL_CLUSTERSCOLOR_PATH = ""



# TSV file giving an hexadecimal-coded color for each HTO.
# Format: no header, HTOs names as first column, colors in second column. All
# existing HTOs must be defined in file (match by name, error otherwise). A
# warning is raised if more colors are declared than than existing HTOs. 
# Empty string (or non-existing file path) for automatic coloring. 
EXTERNAL_HTOSCOLOR_PATH = ""


#### HTO parameters

HTO_METADATA_COLUMN = "HTO_classification" # Name of the metadata column to use for HTO
HTO_FACTOR_LEVELS = NULL
#HTO_FACTOR_LEVELS = c("D21-PBS", "D21-Panth", "D29-PBS", "D29-Panth"); # Order of levels to use for HTOs (for plots, must match with values in selected metadata column). NULL to ignore.


#### Eventual clusters renaming and reordering for figures (by name, NULL to Ignore)

# Simplify clusters names (remove parenthesis information) #
SIMPLIFY_CLUSTER_NAMES = TRUE # One should take this in account in  'ORDER_CLUSTER'

# Specify an ordering for cluster names (NULL to ignore)
ORDER_CLUSTER = c( "B mem 1", "Activated/GC", "Plasma PC2", "GC", "Proliferating", "Plasma PC3", "Plasma PC1", "B mem 2")



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



#### Lists of genes of interest

MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES = NULL


# Path to the Cell cycle gene lists
CELL_CYCLE_SPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "02_Analyses", "01_CellCycle", "S_phase_genes.csv")))
CELL_CYCLE_G2MPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "02_Analyses","01_CellCycle", "G2M_phase_genes.csv")))


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




#### SIGNATURES analysis

# Same approach as in other scripts (modulescore) but also compares conditions
# Warning : as opposed to other scripts, here MODULES_GENES is a 'list of list'
# to allow grouping of signatures

MODULES_GENES = list()

MODULES_GENES = list( "CYCLE" = list("S_PHASE" = CELL_CYCLE_SPHASE_GENELIST, 
                                     "G2M_PHASE" = CELL_CYCLE_G2MPHASE_GENELIST),
                      "SERGIO" = GENES_SERGIO,
                      "OTA_CL5" = GENES_OTA,
                      "CUI_ALL" = GENES_CUI_BY_CYTOKINE,
                      "CUI_T10_PV" = GENES_CUI_BY_CYTOKINE_TOP10_PVAL,
                      "CUI_T50_PV" = GENES_CUI_BY_CYTOKINE_TOP50_PVAL,
                      "CUI_T10_FC" = GENES_CUI_BY_CYTOKINE_TOP10_FC,
                      "CUI_T50_FC"= GENES_CUI_BY_CYTOKINE_TOP50_FC,
                      "LAU_T10_PV" = GENES_LAURINE_TOP10_PVAL,
                      "LAU_T50_PV" = GENES_LAURINE_TOP50_PVAL,
                      "LAU_T10_FC" = GENES_LAURINE_TOP10_FC,
                      "LAU_T50_FC" = GENES_LAURINE_TOP50_FC)


# Remove empty strings
for(currentModule in names(MODULES_GENES)) MODULES_GENES[[currentModule]] = Map( '[', MODULES_GENES[[currentModule]], lapply(MODULES_GENES[[currentModule]], function(x){ which( nchar( x)>0)}));





