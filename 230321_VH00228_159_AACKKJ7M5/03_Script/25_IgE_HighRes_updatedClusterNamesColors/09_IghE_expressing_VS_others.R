# ##############################################################################
# Analyze IgE expression as compared to others
# ##############################################################################



# ..............................................................................
## @knitr ige_vs_others
# ..............................................................................

IGsCompared = c("Ighe", "Igha", "Ighg2c", "Ighg2b", "Ighg1", "Ighg3", "Ighd", "Ighm")

IGsColor = hue_pal()(length(IGsCompared))
names(IGsColor) = IGsCompared

IGs_Expression = GetAssayData( scObject, assay = "RNA", layer = "data")[IGsCompared, ]

# Determine a cell classification by looking at the most expressed heavy chain 
cellIdentity_maxIG = apply(IGs_Expression, 2, function(x)
  {
    # Only assign if there is at least a valid value for one of selected transcripts (NA otherwise)
    return(ifelse(max(x), names(x)[which.max(x)], NA))
  })


kable(addmargins(table(list(cellIdentity_maxIG, scObject[["Identity", drop = TRUE]]), useNA = "ifany")), format = "pipe")

# ..............................................................................
## @knitr ige_expression_distribution
# ..............................................................................

IgheData = data.frame(Expression = IGs_Expression['Ighe',], 
                      Identity = scObject[["Identity", drop = TRUE]], 
                      IgIdentity = cellIdentity_maxIG)

# Select only cells ith a positive value
IgheDataPos = IgheData[IgheData[["Expression"]]>0,]

ggplot(IgheDataPos) +
  geom_histogram(aes(x = Expression)) + 
  ggtitle("Ighe normalized expression (Ighe-positive cells)")

# ..............................................................................
## @knitr ige_expression_distribution_colorCluster
# ..............................................................................

ggplot(IgheDataPos) +
  geom_histogram(aes(x = Expression,
                     fill = Identity)) +
  scale_fill_manual(values = clustersColor) +
  ggtitle("Ighe normalized expression (Ighe-positive cells)")



# ..............................................................................
## @knitr ige_expression_distribution_colorIgIdent
# ..............................................................................

ggplot(IgheDataPos) +
  geom_histogram(aes(x = Expression,
                     fill = IgIdentity)) +
  scale_fill_manual(values = IGsColor) +
  ggtitle("Ighe normalized expression (Ighe-positive cells)")



# ..............................................................................
## @knitr ige_expression_distribution_colorCluster_facetIgIdent
# ..............................................................................

ggplot(IgheDataPos) +
  geom_histogram(aes(x = Expression,
                     fill = Identity)) +
  facet_wrap(vars(IgIdentity)) +
  scale_fill_manual(values = clustersColor) +
  ggtitle("Ighe normalized expression (Ighe-positive cells)")



# ..............................................................................
## @knitr ige_vs_others_table
# ..............................................................................

# Count number of ighe positive expression cells in each previously computed cell identity + clusters
IghePos_nbCells = t( tapply( IGs_Expression["Ighe",], 
                             list( cellIdentity_maxIG, 
                                   scObject[["Identity", drop = TRUE]][names(cellIdentity_maxIG)]), 
                             function(x){sum(x>0)}))

# Replace not-happening combination of IG + Cluster with 0 (no cells => no cells expressing ighe)
IghePos_nbCells[is.na(IghePos_nbCells)] = 0

# Same but divide by number of cells in cluster
nbCellsByCluster = table(scObject[["Identity"]])[rownames(IghePos_nbCells)]
IghePos_pctCluster = apply(IghePos_nbCells, 2, "/", nbCellsByCluster)*100

Heatmap( IghePos_pctCluster,
         name = "%",
         border = 1,
         col = c("white", "darkorange"),
         cluster_rows = FALSE,
         cluster_columns = FALSE,
         row_names_side = "left",
         cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
           grid.text(IghePos_nbCells[i, j], x, y)
         }, 
         column_title = "Number of ighe-transcribing cells\nColor : % of ighe-transcribing cells in cluster")


