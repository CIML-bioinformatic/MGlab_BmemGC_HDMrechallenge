# Merge clusters from TSV file (taken from 01_QC) to be used in downstream analyses

library( funr)
library(forcats);


# Get path of current script
WORKING_DIR = dirname( sys.script())

# Get globalParams variables (project paths)
source( file.path( WORKING_DIR, "..", "globalParams.R"));

# Get path to original clustering output and read as data.frame
inputTSV = file.path( PATH_EXPERIMENT_OUTPUT,
                      "06b_LG_GlobalHeterogeneity_NoT_NoVariableIG_NoHeavyChain",
                      "ClusteringRes_0.4",
                      "Extra",
                      "AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_cellsClusterIdentity.tsv")

idents = read.csv( inputTSV, sep = "\t", row.names = 1);

# Recode factor levels
idents[["identity"]] = fct_recode( factor( idents[["identity"]]), 
                                   "Plasma PC1" = "0", 
                                   "Plasma PC2" = "3", 
                                   "Activated/GC" = "2", 
                                   "B mem 1 (Activated)" = "1",
                                   "B mem 2" = "5",
                                   "Proliferating" = "4");

# Write results
outputTSV = file.path( PATH_EXPERIMENT_OUTPUT,
                       "13b_LG_GroupClusters",
                       "AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_cellsClusterIdentity.tsv")

# Make sure output directory exists
dir.create( dirname( outputTSV))

# Save cells cluster identity as determined with 'FindClusters'
write.table( idents, 
             file = outputTSV, 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");




