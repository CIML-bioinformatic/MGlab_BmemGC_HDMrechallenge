# ##############################################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ##############################################################################


##################################################
# COMPUTATION OF DIMENSIONAL REDUCTION (TSNE/UMAP)
##################################################

# ..............................................................................
## @knitr heterogeneity_dimReduc
# ..............................................................................

#### Compute dimensionality reduction (UMAP and TSNE)

nbPC_dimreduc=DIMREDUC_USE_PCA_NBDIMS
if(RECOMPUTE_PREPROCESSING)
{
  cat("\n\n(Re-)computing tSNE and UMAP...\n")
  
  if(DIMREDUC_USE_PCA_NBDIMS>nbPC)
  {
    warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'dimreduc' (", DIMREDUC_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
    nbPC_dimreduc = nbPC
  }
  
  sc10x = RunUMAP( sc10x, dims = 1:nbPC_dimreduc);
  sc10x = RunTSNE( sc10x, dims = 1:nbPC_dimreduc);
  
  # Update cellsCoordinates for plots that don't rely on Seurat methods (custom plotly & ggplot) 
  if(is.null(cellsCoordinates)) # It has not been overridden by TSV file (and not by preexisting dimreduc in previous object because RECOMPUTE_PREPROCESSING is TRUE)
  {
    cellsCoordinates = Embeddings(Reductions(sc10x, slot = CELLS_COORDINATES))
    
    # PCA can have arbitrary number of dimensions, use requested ones (2d plot)
    if(CELLS_COORDINATES == "pca")
    {
      cellsCoordinates = cellsCoordinates[, PCA_DIMS[1:2]];
    }
    
    useReduction = CELLS_COORDINATES
  }
  
} else
{
  cat("\n\nNot recomputing tSNE and UMAP...\n")
}




#### Save results as tsv files

write.table( Embeddings(sc10x, reduction = "umap"), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_umap.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

write.table( Embeddings(sc10x, reduction = "tsne"), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_tsne.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");
