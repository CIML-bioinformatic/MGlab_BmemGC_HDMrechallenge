

# SCATTER_GENES 
##################

## @knitr scatter_genes
# Create a scatterplot of cells for two selected features
# Can define parentLevel as the number of '#' in parent rmarkdown title so the 
# loop adds paragraph title in the correct (parentLevel+1) level.

# Function creating a scatterplot of selected features scores colored by identity.
# Returns a list with "scatter" ggplot object, and a "tablePos" dataframe giving count and pct of positive cells for each identity
# DEBUG: seuratObject = sc10x;featureNames = SCATTER_GENES[[1]]; assay = "RNA"; slot = "data"; main = ""
genesTableScatter = function(seuratObject, featureNames, assay = "RNA", slot = "data", main = "")
{
  expScores = as.data.frame( t( GetAssayData( object = seuratObject,
                                              slot = slot,
                                              assay = assay)[featureNames[1:2],]))
  if(!(dim(expScores)[2]==2)) stop(paste("Could not retrieve expression values to create scatterplot for selected features:", paste(featureNames, collapse = " - ") ))
  
  #expPos = expScores>0
  #colnames(expPos) = paste0(colnames(expScores), "_pos")
  
  dataPlot = cbind(expScores, ident = Idents(sc10x))
  
  scatter = ggplot( data = dataPlot,
                    aes( x = !!sym(colnames(expScores)[1]),
                         y = !!sym(colnames(expScores)[2]),
                         color = ident,
                         fill = ident)) +
            geom_point() +
            ggtitle(main)

  # Compute counts of positive cells by identity
  identSorted = factor(Idents(sc10x), levels = sort(levels(Idents(sc10x))))  
  posCount = do.call(rbind, by(expScores, identSorted, function(x){colSums(x>0)}))
  
  # Count nb of cells by identity
  ident_cells = as.numeric(table(identSorted)[rownames(posCount)])
  
  # Use previous results to compute percentage of positive cell of each identity
  posPct = round((posCount/ident_cells)*100, digits = 2)
  colnames(posPct) = paste0(colnames(posCount), "_pct_ident")
  
  # Combine stats in a single table
  tablePos = cbind(posCount, posPct, ident_cells)
  
  #ggplot( data = reshape2::melt(posCount),
  #        aes( x = "", y = value, fill = as.character(Var1))) +
  #  geom_bar(stat = "identity") +
  #  coord_polar("y", start = 0) +
  #  facet_grid(rows = Var2, 
  
  #ggMarginal( p,
  #            color = NA,
  #            groupColour = FALSE,
  #            groupFill = TRUE,
  #            type = "histogram",
  #            alpha = 1)
    
  return(list("scatter" = scatter,
              "tablePos" = tablePos))
}

for(scatterName in names(SCATTER_GENES))
{
  cat("\n\n", rep('#', if(exists("parentLevel")) parentLevel+1 else 2), " ", scatterName, " {.tabset .tabset-pills .tabset-fade}\n", sep = "")
  message( paste("Scatter comparison:", scatterName))
  
  ts = genesTableScatter( seuratObject = sc10x,
                          featureNames = SCATTER_GENES[[scatterName]],
                          main = scatterName)
  
  cat("\n\n", rep('#', if(exists("parentLevel")) parentLevel+2 else 3), " Scatterplot\n", sep = "")
  print(ts[["scatter"]])
  
  cat("\n\n", rep('#', if(exists("parentLevel")) parentLevel+2 else 3), " By cluster\n", sep = "")
  print(ts[["scatter"]] + 
          facet_wrap(vars(factor(ident, levels = sort(levels(ident))))) + 
          theme(legend.position = "none"))
  
  cat("\n\n", rep('#', if(exists("parentLevel")) parentLevel+2 else 3), " Table\n", sep = "")
  print( kbl(ts[["tablePos"]], caption = paste0(scatterName, " - positive cells")) %>%
           kable_paper("hover", full_width = F))
  cat("\n")
}

