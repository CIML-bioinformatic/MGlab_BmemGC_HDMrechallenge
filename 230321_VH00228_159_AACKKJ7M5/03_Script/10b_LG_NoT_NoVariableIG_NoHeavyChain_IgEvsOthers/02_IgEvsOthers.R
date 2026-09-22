# ##################################################
# Compare expression of Ighe against other immunoglobulins
# ##################################################



# ..............................................................................
## @knitr ighe_vs_others_total
# ..............................................................................

# Get expression values of genes
expmat = as.matrix(GetAssayData(sc10x, assay = "RNA", slot = "data")[c("Ighe", IgGenesSelection),])


for(currentGene in rownames(expmat)[-1])
{
  message(currentGene)
  
  df = as.data.frame(t(expmat[c("Ighe", currentGene),]))
  df[["Cluster"]] = Idents(sc10x)[rownames(df)]
  df[["HTO"]] = sc10x[["factorHTO", drop = TRUE]][rownames(df)]
  print( ggplot(df, aes_string(x = "Ighe", y = currentGene, col= "Cluster")) +
           geom_point(alpha = 0.4) +
           ggtitle(currentGene))
}


# ..............................................................................
## @knitr ighe_vs_others_byHTO
# ..............................................................................

# Get expression values of genes
expmat = as.matrix(GetAssayData(sc10x, assay = "RNA", slot = "data")[c("Ighe", IgGenesSelection),])


for(currentGene in rownames(expmat)[-1])
{
  cat(paste0("\n## ", currentGene, "\n"))
  message(currentGene)
  
  df = as.data.frame(t(expmat[c("Ighe", currentGene),]))
  df[["Cluster"]] = Idents(sc10x)[rownames(df)]
  df[["HTO"]] = sc10x[["factorHTO", drop = TRUE]][rownames(df)]
  print( ggplot(df, aes_string(x = "Ighe", y = currentGene, col= "Cluster")) +
           geom_point(alpha = 0.4) +
           facet_wrap(vars(HTO)) +
           ggtitle(currentGene))
  cat("\n \n")
}

