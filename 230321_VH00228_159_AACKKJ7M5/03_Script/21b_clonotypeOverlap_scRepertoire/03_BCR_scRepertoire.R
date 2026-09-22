# ##################################################
# BCR analysis (scRepertoire)
# ##################################################



# ..............................................................................
## @knitr screpertoire_load
# ..............................................................................

# For scRepertoire
contigs_csv = loadContigs( dirname(PATH_VDJ_FILE_FILTERED), format = "10X")


# During 'combineBCR()' group.by = "sample" on the full file  was not successful 
# Seurat object has been subset to contain only single sample (mouse)
# Then createHTOContigList() selects from 'contigs_csv' the barcodes also showing in filtered Seurat object, so unrelated clonotypes cannot get grouped, and we can use "group.by = NULL" 

# Create a list expected for downstream analyses, split by Sample (Mouse HTO)
# => Treats mice clonotypes separately
contig.list <- createHTOContigList( contig = contigs_csv[[1]], 
                                    sc10x, 
                                    group.by = "HTOMICE_maxID")

#names(contig.list) # Should reflect an eventual subselection of samples/mice after loading Seurat object

# Manually add sample (Mouse) variable so we we can use group.by argument to compute clusters by sample
#contig.list = Map(function(x, y) {cbind(x, sample = y)}, contig.list, names(contig.list))
#contig.list = addVariable(contig.list, variable.name = "sample", names(contig.list))


# Combine information in 10X csv by clonotype 
# !!! TODO : try to use group.by to process samples separately
combinedBCR = combineBCR(contig.list, 
                         samples = names(contig.list), # As split by mouse above
                         ID=NULL, #names(contig.list), 
                         chain = "both", # The chain to use for clustering when call.related.clones = TRUE. Passed to clonalCluster(). Default is "both".
                         sequence = "nt", # clonotypes based on nucleotidic sequence (instead of amino acids with "aa")	
                         call.related.clones = TRUE, # If TRUE, uses clonalCluster() to identify related clones based on sequence similarity. If FALSE, defines clones by the exact V-gene and CDR3 amino acid sequence.
                         group.by = NULL, #"sample", # The column header used for to group clones. If ('NULL“), clusters will be calculated across samples. Value "sample" ensures that sequences from different samples are never compared or clustered together, even if they are identical.
                         threshold = 0.85, # The similarity threshold passed to clonalCluster() if call.related.clones = TRUE. See ?clonalCluster for details.
                         cluster.method = "components", # The clustering algorithm to use. Defaults to "components", which finds connected subgraphs.
                         use.V = TRUE, # Logical. If TRUE, sequences must share the same V gene to be clustered together.
                         use.J = TRUE, # Logical. If TRUE, sequences must share the same J gene to be clustered together.
                         removeNA = FALSE, # This will remove any chain without values.
                         removeMulti = FALSE, # Logical. If TRUE, removes cells that have more than one distinct heavy or light chain after processing.
                         filterMulti = TRUE, # Logical. If TRUE, filters multi-chain cells to retain only the most abundant IGH and IGL/IGK chains.
                         filterNonproductive = TRUE) # Logical. If TRUE, removes non-productive contigs from the analysis.
                         
# Correct barcode after combination to keep matching with seurat object...
combined = lapply(combinedBCR, function(x){x[["barcode"]] = gsub(".*_", "", x[["barcode"]]); return(x)})
#head(combined[[1]][["barcode"]])


# clonalQuant(combined, 
#             cloneCall="strict", 
#             chain = "both", 
#             scale = TRUE)



# ..............................................................................
## @knitr screpertoire_clonalAbundance
# ..............................................................................

clonalAbundance(combined, 
                cloneCall = "gene", 
                scale = FALSE)

# clonalAbundance(combined, 
#                 cloneCall = "gene", 
#                 scale = TRUE)



# ..............................................................................
## @knitr screpertoire_clonalLength
# ..............................................................................

clonalLength(combined, 
             cloneCall="nt", 
             chain = "both") 



# ..............................................................................
## @knitr screpertoire_clonalCompare
# ..............................................................................

# clonalCompare(combined, 
#               top.clones = 3, 
#               #samples = c("P17B", "P17L"), 
#               cloneCall="nt", 
#               graph = "alluvial")

# clonalScatter(combined, 
#               cloneCall ="gene", 
#               x.axis = "M4", 
#               y.axis = "M6",
#               dot.size = "total",
#               graph = "proportion")



# ..............................................................................
## @knitr screpertoire_clonalHomeostasis
# ..............................................................................

clonalHomeostasis(combined, 
                  cloneCall = "gene")

# clonalHomeostasis(combined, 
#                  cloneCall = "gene",
#                  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, 
#                                Large = 0.3, Hyperexpanded = 1))

# clonalProportion(combined, 
#                  cloneCall = "gene+nt") 

# clonalProportion(combined,
#                  cloneCall = "nt") 

# vizGenes(combined,
#         x.axis = "TRBV",
#         y.axis = "TRBJ",
#         plot = "heatmap",
#         summary.fun = "percent")
  
# percentGenes(combined,
#              chain = "IGH",
#              gene = "Vgene",
#              summary.fun = "percent")

# percentVJ(combined,
#           chain = "TRB",
#           summary.fun = "percent")


# percentAA(combined, 
#           chain = "IGH", 
#           aa.length = 20)



# ..............................................................................
## @knitr screpertoire_positionalEntropy
# ..............................................................................

positionalEntropy(combined, 
                  chain = "IGH", 
                  aa.length = 20)

# clonalDiversity(combined, 
#                 cloneCall = "gene")

# clonalRarefaction(combined,
#                   plot.type = 1,
#                   hill.numbers = 0,
#                   n.boots = 2)

# clonalRarefaction(combined,
#                   plot.type = 2,
#                   hill.numbers = 0,
#                   n.boots = 2)

# clonalRarefaction(combined,
#                   plot.type = 3,
#                   hill.numbers = 0,
#                   n.boots = 2)


# clonalOverlap(combined, 
#               cloneCall = "strict", 
#               method = "morisita")


# clonalOverlap(combined, 
#               cloneCall = "strict", 
#               method = "raw")



# ..............................................................................
## @knitr screpertoire_combineGexBcr
# ..............................................................................

# Create a grouping variable from metadata of interest
sc10x[["clonotypeOverlap_combinedCategories"]] = Reduce( function(x, y){paste(x, y, sep = "_")}, 
                                                         cbind( sc10x[[c("HTOTISSUE_classification")]], 
                                                                cluster = Idents(sc10x)))

# Here "group.by = NULL" as the Seurat abject has been subset to contain only single sample, so clonotypes cannot get mixed
# group.by = "sample" on the full object was not successful 
sc10x = combineExpression( combined, # The product of combineTCR(), combineBCR() or a list of both c(combineTCR(), combineBCR()).
                           sc10x,
                           cloneCall = "strict", # Defines the clonal sequence grouping. Accepted values are: gene (VDJC genes), nt (CDR3 nucleotide sequence), aa (CDR3 amino acid sequence), or strict (VDJC + nt). A custom column header can also be used.
                           chain = "both", # The TCR/BCR chain to use. Use both to include both chains (e.g., TRA/TRB). Accepted values: TRA, TRB, TRG, TRD, IGH, IGL (for both light chains), both.
                           group.by = NULL, # A column header in lists to group the analysis by (e.g., "sample", "treatment"). If NULL, will be based on the list element.
                           proportion = TRUE, # Whether to proportion (TRUE) or total frequency (FALSE) of the clone based on the group.by variable.
                           filterNA = FALSE, # Method to subset Seurat/SCE object of barcodes without clone information
                           cloneSize = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1),
                           addLabel = FALSE) # This will add a label to the frequency header, allowing the user to try multiple group.by variables or recalculate frequencies after subsetting the data.

#FeaturePlot(sc10x, features = c("cloneSize"))


sc10x[["hasVdjInfo"]] = !is.na(sc10x[["CTstrict"]]) # For each cell, use one of the inserted columns to flag whether we have VDJ/colonotype info
pander(table(sc10x[["hasVdjInfo"]]), caption = "Number of cells with clonotype info")


# From IG genes, get information about heavy and light chains main 'isotypes'
# Add it to metadata as in previous manual analysis for comparison
sc10x[["heavy_c_gene"]] = gsub(".*\\.", "", gsub( "_.*", "", sc10x[["CTgene", drop = TRUE]], perl = TRUE), perl = TRUE)
sc10x[["light_c_gene"]] = gsub( ".*\\.", "", gsub( ".*_", "", sc10x[["CTgene", drop = TRUE]], perl = TRUE), perl = TRUE)




# ..............................................................................
## @knitr bcr_clone_size
# ..............................................................................

Seurat::DimPlot(sc10x, group.by = "cloneSize")


# ..............................................................................
## @knitr bcr_chains_stats
# ..............................................................................

pander(addmargins(as.matrix(table(list(factor(sc10x[["heavy_c_gene", drop = TRUE]]), factor(sc10x[["light_c_gene", drop = TRUE]]))))), caption = "Nb cells for each combination of heavy & light chain constant gene")




# ..............................................................................
## @knitr bcr_constant_chains_combination_projection
# ..............................................................................


# Get the values taken by heavy chain in data and filter empty ones  
heavy_chain_values = levels(factor(sc10x[["heavy_c_gene", drop = TRUE]]))
heavy_chain_values = heavy_chain_values[nchar(heavy_chain_values)>0]


for(currentHeavy in c("any_heavy_chain", heavy_chain_values))
{
  cat(paste0("\n### ", currentHeavy, " {.tabset .tabset-fade .tabset-pills}\n\n"))
  
  # Get the values taken by heavy chain in data and filter empty ones  
  light_chain_values = levels(factor(sc10x[["light_c_gene", drop = TRUE]]))
  light_chain_values = light_chain_values[nchar(light_chain_values)>0]
  
  for(currentLight in c("any_light_chain", light_chain_values))
  {
    
    selectHeavy = if(currentHeavy=="any_heavy_chain") rep(TRUE, dim(sc10x)[2]) else grepl(currentHeavy, sc10x[["heavy_c_gene", drop = TRUE]])
    selectLight = if(currentLight=="any_light_chain") rep(TRUE, dim(sc10x)[2]) else grepl(currentLight, sc10x[["light_c_gene", drop = TRUE]])
    
    currentCells = selectHeavy & selectLight
    
    if(sum(currentCells)>0 && (!all(currentCells)))
    {
      cat(paste0("\n#### ", currentLight, "\n"))
      
      # Highlight cells of interest on umap (with count in title)
      p = DimPlot(sc10x, cells.highlight = which(currentCells)) + 
        ggtitle(paste0(currentHeavy, " - ", currentLight," : ", sum(currentCells), " cell(s)")) +
        NoLegend()
      
      cat("\n \n ")
      print(p)
      
      # Same plot but split by HTO
      p = DimPlot(sc10x, cells.highlight = which(currentCells), split.by = "factorHTO", ncol=3) + 
        ggtitle(paste0(currentHeavy, " - ", currentLight," : ", sum(currentCells), " cell(s)")) +
        NoLegend()
      
      cat("\n \n ")
      print(p)
      
      
      # Show a table of cells count by cluster for this chains combination 
      cat("\n \n ")
      cat( pander(addmargins(t(as.matrix(table(list(Cluster = Idents(sc10x)[currentCells], HTO=sc10x[["factorHTO", drop = TRUE]][currentCells]))))), caption = "Cluster distribution for current chain combination" ) )
      cat("\n \n ")
      
    }
    
  }
}



# ..............................................................................
## @knitr bcr_clonotypes_stats
# ..............................................................................
# Show number of cells by clonotypes and the ID of most abundant ones

# Count number of cells for each colonotype (distribution of clonotypes represented in given number of cells)
countCellsByClonotype = as.data.frame.table( table( clonotype = sc10x[["CTstrict", drop = TRUE]]), 
                                             responseName = "nbCells")


# Plot as histogram
histoCellsClonotype = ggplot( countCellsByClonotype) +
  geom_histogram(aes(x=nbCells)) +
  ylab("Clonotypes count")


print(ggplotly(histoCellsClonotype))

cat("\n \n ")
# Create a table
pander(( as.data.frame.table( table( nbCells = table( colonotype = sc10x[["CTstrict"]])),
                              responseName = "nbClonotypes")))

cat("\n \n ")


# Show a table of the clonotypes in more than 'n' cells (to focus on the right of the distribution)
nbCellsTopClonotypes = 5

# Isolate the clonotype counts represented in multiple cells and sort by number of cells
topClonotypes = countCellsByClonotype[countCellsByClonotype[["nbCells"]]>=nbCellsTopClonotypes,]
topClonotypes = topClonotypes[order(topClonotypes[["nbCells"]], decreasing = TRUE),]

# Isolate each of these 'top clonotypes' in object sc10x and compute distribution based on clusters
topClonotypesDistributionByCluster = lapply(topClonotypes[["clonotype"]], function(currentClonotype)
{
  currentCells = which(sc10x[["CTstrict", drop = TRUE]]==currentClonotype)
  table(Cluster = Idents(sc10x)[currentCells])
})

# Combine all results with previously existing data.frame
topClonotypes = cbind(topClonotypes,
                      do.call(rbind, topClonotypesDistributionByCluster))

# Show as dataframe on reports
topClonotypesDT = datatable( topClonotypes, 
                             caption = paste0( "Top clonotypes (down to ", nbCellsTopClonotypes, " cells)"),
                             rownames = FALSE, 
                             height = 800, 
                             width = 800)
htmltools::tagList(print(topClonotypesDT))


# ..............................................................................
## @knitr bcr_top_clonotypes_projection
# ..............................................................................

cat("\n \n ")

for(currentClonotype in as.character(topClonotypes[["clonotype"]]))
{
  currentCells = which(sc10x[["CTstrict", drop = TRUE]]==currentClonotype)
  
  p = DimPlot(sc10x, cells.highlight = currentCells) + 
    ggtitle(currentClonotype) +
    NoLegend()
  
  #table(Idents(sc10x)[currentCells])
  
  print(p)
  
}

cat("\n \n ")



################################################################################
################################################################################
################################################################################
################################################################################



# ..............................................................................
## @knitr clonotype_counts_heatmap
# ..............................................................................

# Eventual subsetting
clonotype_colname = "CTstrict"
scObject = sc10x

# Detect clonotypes showing up in multiple cells
clonoTable = sort( table( scObject[[clonotype_colname]]), decreasing = TRUE)
clonoMultiple = names( clonoTable)[clonoTable>1]

# Inject categories in metadata
scObject[["clonotypeOverlap_combinedCategories"]] = Reduce(function(x, y){paste(x, y, sep = "_")}, cbind(scObject[[c("HTOMICE_classification", "HTOTISSUE_classification")]], cluster = Idents(scObject)))

# Select cells with a defined clonotype to track for overlap in cells groups
cellsWithClonotype = !is.na(scObject[[clonotype_colname, drop = TRUE]])
cellsWithHTO = !grepl("Negative", scObject[["HTO_classification", drop = TRUE]])
cellsWithClonotypeAndHTO_barcode = names(cellsWithClonotype)[cellsWithClonotype & cellsWithHTO]




# Compute the count of each clonotype in each category (as defined earlier)
clonotype_categories = as.data.frame( table( list(category = scObject[["clonotypeOverlap_combinedCategories", drop = TRUE]][cellsWithClonotypeAndHTO_barcode], clonotype = scObject[[clonotype_colname, drop = TRUE]][cellsWithClonotypeAndHTO_barcode])))

# Select colnotypes found at least twice
clonotype_categories_multipleOnly = clonotype_categories[clonotype_categories[["clonotype"]] %in% clonoMultiple, ]

# Convert to wide and move clonotype names to rownames
clonotype_categories_matrix = spread(clonotype_categories_multipleOnly, "category", "Freq")
rownames(clonotype_categories_matrix) = clonotype_categories_matrix[["clonotype"]]
clonotype_categories_matrix[["clonotype"]] = NULL

# Remove categories where no selected clonotypes found
clonotype_categories_matrix = clonotype_categories_matrix[, colSums(clonotype_categories_matrix)>0]


# Reorder to put larger clonotypes on top
#clonotype_categories_matrix = clonotype_categories_matrix[order(rowSums(clonotype_categories_matrix), decreasing = TRUE),]
# Reorder to put clonotype showing up in most categories first
clonotype_categories_matrix = clonotype_categories_matrix[order(rowSums(clonotype_categories_matrix>0), decreasing = TRUE),]

# Replace zeros by NA for graphic representation
clonotype_categories_matrix_NAs = clonotype_categories_matrix
clonotype_categories_matrix_NAs[clonotype_categories_matrix==0] = NA

# Select top clonotypes and remove empty columns
clonotype_categories_matrix_NAs = clonotype_categories_matrix_NAs[1:min(nrow(clonotype_categories_matrix_NAs), 50),]
clonotype_categories_matrix_NAs = clonotype_categories_matrix_NAs[, colSums(clonotype_categories_matrix_NAs, na.rm = TRUE)>0]






anno_bar_left = function(x) {
  max_x = max(x)
  cell_fun_pct = function(i) {
    pushViewport(viewport(xscale = c(0, max_x)))
    grid.roundrect(x = unit(1, "npc"), width = unit(x[i], "native"), 
                   height = unit(1, "npc") - unit(4, "pt"), 
                   just = "right", gp = gpar(fill = "#FF000080", col = NA))
    if(as.numeric(convertWidth(grobWidth(textGrob(x[i], gp = gpar(), rot = 0)), "native"))  + as.numeric(convertWidth(unit(4, "pt"), "native")) > as.numeric(unit(x[i], "native")))
      grid.text(x[i], x = unit(1, "npc") - unit(x[i], "native") - unit(2, "pt"), just = "right")
    else
      grid.text(x[i], x = unit(1, "npc") - unit(x[i], "native") + unit(2, "pt"), just = "left")
    convertWidth(grobWidth(textGrob(x[i], gp = gpar(), rot = 0)), "native")
    
    popViewport()
  }
  
  AnnotationFunction(
    cell_fun = cell_fun_pct,
    var_import = list(max_x, x), 
    which = "row",
    width = max_text_width(x)*5
  )
}

anno_bar_right = function(x) {
  max_x = max(x)
  cell_fun_pct = function(i) {
    pushViewport(viewport(xscale = c(0, max_x)))
    grid.roundrect(x = unit(0, "npc"), width = unit(x[i], "native"), 
                   height = unit(1, "npc") - unit(4, "pt"), 
                   just = "left", gp = gpar(fill = "#0000FF80", col = NA))
    if(as.numeric(convertWidth(grobWidth(textGrob(x[i], gp = gpar(), rot = 0)), "native")) + as.numeric(convertWidth(unit(4, "pt"), "native")) > as.numeric(unit(x[i], "native")))
      grid.text(x[i], x = unit(x[i], "native") + unit(2, "pt"), just = "left")
    else
      grid.text(x[i], x = unit(x[i], "native") - unit(2, "pt"), just = "right")
    popViewport()
  }
  
  AnnotationFunction(
    cell_fun = cell_fun_pct,
    var_import = list(max_x, x), 
    which = "row",
    width = max_text_width(x)*5, 
  )
}




# Prepare rows annotations
rowsAnnotation = HeatmapAnnotation(#clonotype_sumNbCells = rowSums(clonotype_categories_matrix_NAs, na.rm = TRUE),
  #clonotype_sumNbcategories = rowSums(as.matrix(clonotype_categories_matrix_NAs)>0, na.rm = TRUE),
  nbCells_LN = anno_bar_left(rowSums(as.matrix(clonotype_categories_matrix_NAs[,grepl("_LN", colnames(clonotype_categories_matrix_NAs))]), na.rm = TRUE)),
  nbCells = anno_text(rowSums(as.matrix(clonotype_categories_matrix_NAs), na.rm = TRUE), show_name = TRUE, just = "center", location = unit(0.5, 'npc')),
  nbCells_LG = anno_bar_right(rowSums(as.matrix(clonotype_categories_matrix_NAs[,grepl("_LG", colnames(clonotype_categories_matrix_NAs))]), na.rm = TRUE)),
  nbCategories_LN = anno_bar_left(rowSums(as.matrix(clonotype_categories_matrix_NAs[,grepl("_LN", colnames(clonotype_categories_matrix_NAs))])>0, na.rm = TRUE)),
  nbCategories = anno_text(rowSums(as.matrix(clonotype_categories_matrix_NAs)>0, na.rm = TRUE), show_name = TRUE, just = "center", location = unit(0.5, 'npc')),
  nbCategories_LG = anno_bar_right(rowSums(as.matrix(clonotype_categories_matrix_NAs[,grepl("_LG", colnames(clonotype_categories_matrix_NAs))])>0, na.rm = TRUE)),
  which = "row")


colsAnnotation = HeatmapAnnotation(Mouse = gsub("_.*", "", colnames(clonotype_categories_matrix_NAs)),
                                   Tissue = gsub("(.*)_(.*)_.*", "\\2", colnames(clonotype_categories_matrix_NAs)),
                                   col = list(Mouse = MOUSE_COLORS,
                                              Tissue = TISSUE_COLORS))



Heatmap(clonotype_categories_matrix_NAs,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = gsub("(.*)_(.*)_.*", "\\1_\\2", colnames(clonotype_categories_matrix_NAs)),
        cell_fun = function(j, i, x, y, w, h, col){if(!is.na(clonotype_categories_matrix_NAs[i, j])) grid.text(clonotype_categories_matrix_NAs[i, j], x, y)},
        na_col = "grey",
        col = colorRamp2(range(clonotype_categories_matrix_NAs, na.rm = TRUE), hcl_palette = "Mint", reverse = TRUE), 
        column_title_rot = 90,
        left_annotation = rowsAnnotation,
        top_annotation = colsAnnotation)




# ..............................................................................
## @knitr clonotype_overlap_currentmouse
# ..............................................................................

clonotype_colname = "CTstrict"

# Overlap coefficient (Szymkiewicz–Simpson coefficient)
overlapIndex <- function(a, b) {
  a = unique(a)
  b = unique(b)
  intersection = length(intersect(a, b))
  #union = length(a) + length(b) - intersection
  return (intersection/min(c(length(a), length(b))))
}

#for( in rev(sort(unique( sc10x[["HTOMICE_classification", drop = TRUE]]))))
#{

  currentMouse = CURRENT_MOUSE

  #cat("\n####", currentMouse, "\n", sep = " ")
  
  fromCurrentMouse = (sc10x[["HTOMICE_classification", drop = TRUE]] == currentMouse)
  # Actually a simple combined factor is enough as we do all combinatory comparisons
  clonotypeNames_byClusterAndTissue = tapply(sc10x[[clonotype_colname, drop = TRUE]][fromCurrentMouse], 
                                             INDEX = paste(sc10x[["HTOTISSUE_classification", drop = TRUE]][fromCurrentMouse], Idents(sc10x)[fromCurrentMouse], sep = "_"), function(x){ unique(x[!is.na(x)])})
  
  # Remove eventual groups without clonotype defined
  clonotypeNames_byClusterAndTissue = clonotypeNames_byClusterAndTissue[sapply(clonotypeNames_byClusterAndTissue, length)>0]
  
  res = matrix(NA, nrow = length(clonotypeNames_byClusterAndTissue), ncol = length(clonotypeNames_byClusterAndTissue))
  rownames(res ) = names(clonotypeNames_byClusterAndTissue)
  colnames(res ) = names(clonotypeNames_byClusterAndTissue)
  
  for(i in rownames(res ))
    for(j in colnames(res)) res[i,j] = overlapIndex(clonotypeNames_byClusterAndTissue[[i]], clonotypeNames_byClusterAndTissue[[j]] )
  
  print(Heatmap( res,
                 #rect_gp = gpar(col = "white", lwd =0.5),
                 name = "Clonal ovelap index",
                 col = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(3),
                 cluster_columns = T,
                 cluster_rows = T, 
                 column_title = currentMouse))
  
  cat("\n \n")
#}



# ..............................................................................
## @knitr clonotype_overlap_customSergio
# ..............................................................................
# Customized figure with subselection of clusters and mouse 456 only, with 
# clusters renaming on the fly. (individual mice and grouped in same chunk)

clonotype_colname = "CTstrict"

# Save identities as we need to override names temporarily
backupIdents = Idents(sc10x)

# Rename concerned clusters
Idents(sc10x) = fct_recode( Idents(sc10x), 
                            "MBCs" = "B mem 1", 
                            "Mature PCs" = "Plasma PC1",
                            "Immature PC1" = "Plasma PC2",
                            "Immature PC2" = "Plasma PC3");

# Individual mice

#for(currentMouse in c("M4", "M5", "M6"))
#{

  currentMouse = CURRENT_MOUSE

  #cat("\n#### Custom:", currentMouse, "\n", sep = " ")
  
  fromCurrentMouse = (sc10x[["HTOMICE_classification", drop = TRUE]] == currentMouse &
                        Idents(sc10x) %in% c("MBCs", "Mature PCs", "Immature PC1", "Immature PC2", "Proliferating"))
  
  # Actually a simple combined factor is enough as we do all combinatory comparisons
  clonotypeNames_byClusterAndTissue = tapply(sc10x[[clonotype_colname, drop = TRUE]][fromCurrentMouse], 
                                             INDEX = paste(sc10x[["HTOTISSUE_classification", drop = TRUE]][fromCurrentMouse], Idents(sc10x)[fromCurrentMouse], sep = "_"), function(x){ unique(x[!is.na(x)])})
  
  # Remove eventual groups without clonotype defined
  clonotypeNames_byClusterAndTissue = clonotypeNames_byClusterAndTissue[sapply(clonotypeNames_byClusterAndTissue, length)>0]
  
  res = matrix(NA, nrow = length(clonotypeNames_byClusterAndTissue), ncol = length(clonotypeNames_byClusterAndTissue))
  rownames(res ) = names(clonotypeNames_byClusterAndTissue)
  colnames(res ) = names(clonotypeNames_byClusterAndTissue)
  
  for(i in rownames(res ))
    for(j in colnames(res)) res[i,j] = overlapIndex(clonotypeNames_byClusterAndTissue[[i]], clonotypeNames_byClusterAndTissue[[j]] )
  
  print(Heatmap( res,
                 #rect_gp = gpar(col = "white", lwd =0.5),
                 name = "Clonal ovelap index",
                 col = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(3),
                 cluster_columns = T,
                 cluster_rows = T, 
                 column_title = currentMouse))
  
  cat("\n \n")
#}

Idents(sc10x) = backupIdents
  





################################################################################
################################################################################
################################################################################
################################################################################


# ..............................................................................
## @knitr bcr_clonotypes_circlize_ALL
# ..............................................................................

# clonalOverlay(sc10x, 
#               reduction = "umap", 
#               cutpoint = 1, 
#               bins = 10, 
#               facet.by = "HTOTISSUE_classification") + 
#   guides(color = "none")

circles = getCirclize(sc10x, 
                      #group.by = "clonotypeOverlap_combinedCategories",
                      group.by = "Identity", 
                      proportion = FALSE,  
                      cloneCall = "strict", 
                      include.self = TRUE)


chordDiagram(circles, 
             self.link = 2, 
             grid.col = clustersColor)
             #directional = 1)#, 
             #direction.type =  "arrows",
             #link.arr.type = "big.arrow")


# ..............................................................................
## @knitr bcr_clonotypes_circlize_LN
# ..............................................................................

currentTissue = "LN"
scObject = subset(sc10x, subset = HTOTISSUE_maxID == currentTissue)

circles = getCirclize(scObject, 
                      # group.by = "clonotypeOverlap_combinedCategories",
                      group.by = "Identity", 
                      proportion = FALSE,  
                      cloneCall = "strict", 
                      include.self = TRUE)

chordDiagram(circles, 
             self.link = 2, 
             grid.col = clustersColor)#, 
             #directional = 1, 
             #direction.type =  "arrows",
             #link.arr.type = "big.arrow")



# ..............................................................................
## @knitr bcr_clonotypes_circlize_LG
# ..............................................................................

currentTissue = "LG"
scObject = subset(sc10x, subset = HTOTISSUE_maxID == currentTissue)

circles = getCirclize(scObject, 
                      group.by = "Identity",
                      proportion = FALSE,  
                      cloneCall = "strict", 
                      include.self = TRUE)

chordDiagram(circles, 
             self.link = 2, 
             grid.col = clustersColor)#, 
             #directional = 1, 
             #direction.type =  "arrows",
             #link.arr.type = "big.arrow")


# ..............................................................................
## @knitr bcr_clonotypes_clonalNetwork_LN_LG
# ..............................................................................

currentTissue = "LN"
scObject = subset(sc10x, subset = HTOTISSUE_maxID == currentTissue)

clonalNetwork(scObject,
              cloneCall = "strict",
              chain = "both",
              reduction = "umap",
              group.by = "clonotypeOverlap_combinedCategories",
              filter.clones = NULL,
              filter.identity = NULL,
              filter.proportion = NULL,
              filter.graph = TRUE,
              exportClones = FALSE,
              exportTable = FALSE,
              palette = "inferno") + ggtitle(currentTissue)


currentTissue = "LG"
scObject = subset(sc10x, subset = HTOTISSUE_maxID == currentTissue)


clonalNetwork(scObject,
              cloneCall = "strict",
              chain = "both",
              reduction = "umap",
              group.by = "clonotypeOverlap_combinedCategories",
              filter.clones = NULL,
              filter.identity = NULL,
              filter.proportion = NULL,
              filter.graph = TRUE,
              exportClones = FALSE,
              exportTable = FALSE,
              palette = "inferno") + ggtitle(currentTissue)



# ..............................................................................
## @knitr bcr_clonotypes_clonalBias
# ..............................................................................

clonalBias(scObject, 
           cloneCall = "strict", 
           split.by = "HTOTISSUE_maxID", 
           group.by = "Identity",
           n.boots = 10, 
           min.expand =0)



# ..............................................................................
## @knitr bcr_clonotypes_venn
# ..............................................................................

# Split clonotype clustering by tissue
listClonotypeByTissue = split(sc10x[['CTstrict', drop = TRUE]], sc10x[['HTOTISSUE_maxID', drop = TRUE]])

vennTable = venn(lapply(lapply(listClonotypeByTissue, unique), na.exclude))
#attr(vennTable, "intersections")

# Save markers list as 'tsv' table
write.table( vennTable,
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "vennTable.tsv")),
             quote = FALSE,
             row.names = TRUE, 
             #col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");


# ..............................................................................
## @knitr bcr_clonotypes_venn_merge
# ..............................................................................

# Search csv files in other samples results folders (generated by loop in launch_report_compilation.R) 
precomputedVenn_otherSamples_path = dir(path = dirname(PATH_ANALYSIS_OUTPUT), pattern = ".*vennTable.tsv", recursive = TRUE, full.names = TRUE)
# debug: precomputedVenn_otherSamples_path = dir(path = (PATH_ANALYSIS_OUTPUT), pattern = ".*vennTable.tsv", recursive = TRUE, full.names = TRUE)
names(precomputedVenn_otherSamples_path) = gsub('.*_', '', dirname(precomputedVenn_otherSamples_path)) # basename(dirname(precomputedVenn_otherSamples_path))

# Read previously computed results
csvDataList = lapply(precomputedVenn_otherSamples_path, read.table, sep = "\t")
# Rename count column with sample name
csvDataList = Map(function(data, sampleName){ dimnames(data)[[2]][1] = sampleName; return(data)}, 
                  csvDataList, 
                  names(csvDataList))

# Merge data from all list elements
csvDataMerge = Reduce(function(x,y){merge(x,y,by = c("LN", "LG"), )}, csvDataList)

# Select combination matching at least one tissue (remove 0 0)
csvDataMerge = csvDataMerge[csvDataMerge[["LN"]] | csvDataMerge[["LG"]],]

# Give rowname based on Tissue, format metadata, removed useless columns
rownames(csvDataMerge) = gsub("^\\+|\\+$","", paste(ifelse(csvDataMerge[["LN"]], "LN", ""), 
                                                ifelse(csvDataMerge[["LG"]], "LG", ""), 
                                                sep = "+"))
csvDataMerge = csvDataMerge[3:ncol(csvDataMerge)]

# Order tissue levels for plotting
csvDataMerge[["Tissue"]] = factor(rownames(csvDataMerge), levels = c("LN", "LN+LG", "LG"))

# Convert to long plot
csvDataLong = reshape2::melt(csvDataMerge)
ggplot(csvDataLong, aes(x = variable, y=value, fill = Tissue, label=value)) +
  geom_bar(stat = 'identity', alpha = 0.9) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  ylab("Clonal Identity") +
  xlab(NULL) +
  scale_fill_manual(values = TISSUE_COLORS)  

