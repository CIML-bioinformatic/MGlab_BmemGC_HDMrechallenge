# ##############################################################################
# Analyze IgE expression 
# ##############################################################################


# ..............................................................................
## @knitr plot_umap_ige_positive_cells_prepare_data
# ..............................................................................

# Create a copy of cellsData for this specific chunks 
cellsData_IgE = cellsData

# Add the IgE expression values to 'cellsData' dataframe
cellsData_IgE[["IgE"]] = GetAssayData( sc10x, assay = "RNA", layer = "data")["Ighe",rownames(cellsData)]
cellsData_IgE[["IgE_pos_HTO"]] = ifelse(cellsData_IgE[["IgE"]]>0, as.character(cellsData_IgE[["HTO"]]), NA )

# Reorder dataset to get IgE positive cells plotted on top
cellsData_IgE = cellsData_IgE[order(cellsData_IgE[["IgE"]]>0, decreasing = FALSE), ]



# ..............................................................................
## @knitr plot_umap_ige_positive_cells
# ..............................................................................

ggFigure = ggplot(  cellsData_IgE, 
                    aes( x = UMAP_1, 
                         y = UMAP_2,
                         color = gsub(".*_", "", IgE_pos_HTO),
                         alpha = IgE>0,
                         size = IgE>0,
                         shape = gsub("_.*", "", HTO))) +
  geom_point() +
  scale_size_manual(values = c("TRUE" = 3.5, "FALSE"= 1)) +
  scale_color_manual( name = "mouse", values = c("M1" = "#E90000", "M2" = "#A40D03", "M3" = "#7D0000", "M4" = "#4ACADE", "M5" = "#38A6D9", "M6" = "#4076D9")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE"= 0.15)) +
  scale_shape_manual( name = "tissue", values = c("LG" = 17, "LN"= 19))

print( ggFigure);



# ..............................................................................
## @knitr plot_umap_ige_positive_cells_split_tissue
# ..............................................................................

ggFigure = ggplot(  cellsData_IgE, 
                    aes( x = UMAP_1, 
                         y = UMAP_2,
                         color = gsub(".*_", "", IgE_pos_HTO),
                         alpha = IgE>0,
                         size = IgE>0,
                         shape = gsub("_.*", "", HTO))) +
  geom_point() +
  scale_size_manual(values = c("TRUE" = 3.5, "FALSE"= 1)) +
  scale_color_manual( name = "mouse", values = c("M1" = "#E90000", "M2" = "#A40D03", "M3" = "#7D0000", "M4" = "#4ACADE", "M5" = "#38A6D9", "M6" = "#4076D9")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE"= 0.15)) +
  scale_shape_manual( name = "tissue", values = c("LG" = 17, "LN"= 19)) +
  facet_wrap(vars( gsub("_.*", "", HTO)))


print( ggFigure);




# ..............................................................................
## @knitr plot_umap_ige_positive_cells_color_by_tissue
# ..............................................................................

ggFigure = ggplot(  cellsData_IgE, 
                    aes( x = UMAP_1, 
                         y = UMAP_2,
                         color = ifelse(IgE>0, gsub("_.*", "", HTO), NA),
                         alpha = IgE>0,
                         size = IgE>0)) +
  geom_point() +
  scale_size_manual(values = c("TRUE" = 3, "FALSE"= 1)) +
  scale_color_manual( name = "tissue", values = c("LG" = "#A40D03", "LN" = "#38A6D9")) +
  scale_alpha_manual(values = c("TRUE" = 0.75, "FALSE"= 0.15))


print( ggFigure);


# ..............................................................................
## @knitr plot_umap_ige_positive_cells_split_mouse_color_by_cluster
# ..............................................................................

ggFigure = ggplot(  cellsData, # Use the original cellsData that does not contain IgE additional information so fact_wrap puts it in all panels
                    aes(x = UMAP_1, 
                        y = UMAP_2)) +
  geom_point( color="grey",
              alpha = 0.1) +
  geom_point( data = cellsData_IgE[cellsData_IgE[["IgE"]]>0,], 
              #aes( color = ifelse(IgE>0, gsub("_.*", "", HTO), NA))) +
              aes( color = Cluster)) +
  #scale_color_manual( name = "tissue", values = c("LG" = "#A40D03", "LN" = "#38A6D9")) +
  scale_color_manual( name = "cluster", values = clustersColor) +
  facet_wrap( vars( IgE_pos_HTO))

print( ggFigure);


