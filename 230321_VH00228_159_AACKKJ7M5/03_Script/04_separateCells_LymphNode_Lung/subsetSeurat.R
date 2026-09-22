# Load seurat object and subset it based on cell identity as specified in tsv file
# Used for downstream analysis when clustering is not taken as resulting from 
# subclustering script (which gives a subset Seurat object too).


library( funr)
library( forcats);
library( Seurat);


# Get path of current script
WORKING_DIR = dirname( sys.script())
# debug: WORKING_DIR = getwd()


# Get globalParams variables (project paths)
source( file.path( WORKING_DIR, "..", "globalParams.R"));

# Get path to original seurat object (take it from subclustering results)
inputSeurat = file.path( PATH_EXPERIMENT_OUTPUT,
                         "02_Demux",
                         "ClusteringRes_0.6",
                         "AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_seuratObject_final.RDS") # Could have taken the one with all cells from 01_QC but this one has new umap coordinates inside (no need to use external csv file for coordinates)

# Load Seurat object
sc10x = readRDS(inputSeurat)



# In this version we only use HTO metadata to filter cells, we don' load file of cell identity as usually dobe for subsetting 


### Define subgroups of interest
identsList = split(colnames(sc10x), sc10x[["HTOTISSUE_classification", drop = TRUE]]);





###
# Subset seurat object and save result

# Make sure output directory exists
outputFolder = file.path( PATH_EXPERIMENT_OUTPUT,
                          "04_separateCells_LymphNode_Lung")
dir.create( outputFolder)

for(currentSetName in names(identsList))
{
  message(currentSetName)
  # Subset the original object using cell names
  resultSubset = sc10x[, identsList[[currentSetName]]]
  
  # Save as binary file for downstream analyses
  saveRDS( object = resultSubset,
           file = file.path(outputFolder, paste0("AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge_230321_VH00228_159_AACKKJ7M5_seuratObject_subset_", currentSetName, ".RDS")))
  
}



