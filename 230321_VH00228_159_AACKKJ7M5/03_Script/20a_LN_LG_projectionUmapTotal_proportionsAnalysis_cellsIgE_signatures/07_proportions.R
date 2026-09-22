# ##############################################################################
# Performing an analysis of cell idntity representation by condition using chi2
# and sccomp
# ##############################################################################



# ..............................................................................
## @knitr plot_cells_proportions
# ..............................................................................

# Compute percentage of cells by mice for each cluster
clustersCountByHTO_pctByMice = apply( clustersCountByHTO, 2, function(x){ 100*x/sum(x) })

longtable = reshape2::melt(clustersCountByHTO_pctByMice, value.name = "PctClusterByTissue")

longtable[["Tissue"]] = gsub( "_.*", "", longtable[["factorHTO"]])
longtable[["Mouse"]]  = gsub( ".*_", "", longtable[["factorHTO"]])

ggplot( longtable, aes( x=Tissue, y=PctClusterByTissue , col = Mouse)) + 
  geom_point() + 
  facet_wrap( vars( Identity), nrow = 1) +
  geom_line( aes( group = Mouse), 
             lty = 2,
             lwd = 0.25) +
  scale_color_manual(values = c("M1" = "#E90000", "M2" = "#A40D03", "M3" = "#7D0000", "M4" = "#4ACADE", "M5" = "#38A6D9", "M6" = "#4076D9")) +
theme_bw()


#######
# CHI 2
#######

# ..............................................................................
## @knitr chisquare_test_clusters_HTOgroupDay
# Test independance between rows and columns (clusters and sample condition/HTO)
# to detect counts over expected values under random condition.
# ..............................................................................

# Recompute counts for HTOs grouped by day instead of individual mice samples
clustersCountByHTO_groupDay = table( data.frame( "Identity" = as.character( Idents( sc10x)), 
                                                 "Category" = gsub( "M4|M5|M6", "D7", gsub( "M1|M2|M3", "D0", sc10x[["factorHTO", drop = TRUE]])) ))



chisq = chisq.test(clustersCountByHTO)

# Show the expected counts considering margins statistics, to contrast with
# observed counts shown above
print( knitr::kable( addmargins(round(chisq[["expected"]],2)),
                     align   = "c",
                     caption = "Expected number of cells by population (Cluster) and condition (HTO) from X².\nTo compare with observed counts.")) #%>% kable_styling(bootstrap_options = c("condensed"), full_width = FALSE);

# Show result as simple table
pander(chisq)


# ..............................................................................
## @knitr chisquare_test_clusters_HTO_plotResidualsAndContributions
# ..............................................................................
# Set a layout to combine pearson residuals and contribution plots
def.par <- par(no.readonly = TRUE) # save default, for resetting...
layout(matrix(c(1,2), nrow=2, byrow = TRUE)) # Two next figures on same row

# Plot Pearson residuals
corrplot(chisq$residuals, 
         is.cor=FALSE, 
         addCoef.col="darkgrey", 
         tl.srt=60, 
         cl.ratio=0.35, 
         col = rev(colorRampPalette(RColorBrewer::brewer.pal(name="RdBu", 11))(100)),
         #col = colorRampPalette(c("red", "white", "blue"))(100),
         mar=c(0,0,2,0), 
         title = "Pearson residuals")

# Plot contributions (see: http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r)
contrib = round(100*(chisq$residuals^2)/chisq$statistic,3)
corrplot(contrib, 
         is.cor=FALSE, 
         addCoef.col="darkgrey", 
         tl.srt=60, 
         cl.ratio=0.35, 
         col = colorRampPalette(c("white", "yellow", "darkorange", "darkred"))(100),
         mar=c(0,0,2,0), 
         title = "Contribution to X² statistic (%)")

# Reset layout to previously saved default
par(def.par)



# ..............................................................................
## @knitr chisquare_test_clusters_HTO_plotMosaic
# ..............................................................................
mosaicplot(clustersCountByHTO, shade = TRUE, las = 2, main = "X² residuals", xlab="", ylab="")



########
# SCCOMP
########

# ..............................................................................
## @knitr sccomp_computation
# ..............................................................................

# ##############################################################
# This script aims to analyse the composition of clusters
# across condition and time
# ##############################################################

# .....................................................
## @knitr composition_analysis_by_tissue
# .....................................................

# See https://github.com/stemangiola/sccomp
  
cat("\n \n")
cat("#### Barplot")
cat("\n \n")

# Look at the distribution of condition across clusters (group previous)
tissue_vs_cluster_table = table(cluster = Idents(sc10x), tissue = gsub("_.*", "", as.character(sc10x[["factorHTO", drop = TRUE]])))

print(
  ggplot( reshape2::melt(tissue_vs_cluster_table), aes(fill=cluster, y=value, x=tissue)) + 
    geom_bar(position = 'fill', stat="identity") + 
    scale_fill_manual(values = clustersColor)
)

cat("\n \n")
cat("#### Composition analysis")
cat("\n \n")

# Make sure seurat cluster identity is available in "seurat_clusters"
sc10x[["seurat_clusters"]] = Idents(sc10x)
sc10x[["tissue"]] = gsub("_.*", "", as.character(sc10x[["factorHTO", drop = TRUE]]))

# Analyze the data composition by clusters against conditions
sccomp_result_var = sc10x |>
  sccomp_estimate( 
    formula_composition = ~ tissue, 
    formula_variability = ~ tissue,
    .sample =  factorHTO, 
    .cell_group = seurat_clusters,
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers() |> 
  sccomp_test(test_composition_above_logit_fold_change = 0.2)

# Print the result table
print(
  sccomp_result_var %>%kable() %>% 
    kable_styling( full_width = FALSE) %>%
    kable_paper("hover", full_width = T) %>%
    scroll_box( height = "500px")
)


# Get the plots from the results
plots = plot( sccomp_result_var) 

# Plot the posterior predictive check
print( plots$boxplot)

# Plot the 1D significance plots
print( plots$credible_intervals_1D)

# Save result to RDS and plot to SVG
saveRDS( object = sccomp_result_var, file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "CompositionAnalysis_ByCondition_ScCompObject.RDS")))

ggsave( filename = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "CompositionAnalysis_ByCondition_IntervalPlot.svg")),
        plot = plots$credible_intervals_1D)


