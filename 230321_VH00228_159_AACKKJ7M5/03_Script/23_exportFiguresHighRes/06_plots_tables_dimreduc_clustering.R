# ##############################################################################
# Plotting chunks show cells on a 2D plot (e.g. dimensionality reduction). Some
# use javascript link between report objects to select cells in plots and show
# their characteristics on attached graphs or tables. 
# Desired cells coordinates must be extracted in a  dataframe/matrix named
# 'cellsCoordinates' (see 01_prepareData for creation).
# Individual cell properties are extracted from Seurat object ("numID",
# "nCount_RNA", "nFeature_RNA", and "percent.mito", "percent.ribo" if available)
# Cluster identity is extracted from current Idents of Seurat object, and HTO
# category from metadata "factorHTO" which must be added previously.
# The plotting chunks also use the named vectors 'clustersColor' and 'HTOsColor'
# to configure consistent colors.
# It is possible to set variables 'showDimReducClusterLabels' and 'showLegend'
# to control optional settings but a fallback behavior is planned if they don't
# exist.
# ##############################################################################

# This is a stripped-down version imported from anlysis step 14, to focus on 
# plotting umaps with updated cluster nomenclature and colors


# ..............................................................................
## @knitr dimreduc_plotsPreparation
## Prepare SharedData object and variables for plots 
# ..............................................................................

# For each cell, gather all data to be used for plots generation
cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                    sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                    "percent.mito" = if("percent.mito" %in% colnames(sc10x[[]])) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                    "percent.ribo" = if("percent.ribo" %in% colnames(sc10x[[]])) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                    "Cluster" = Idents( sc10x),
                    "HTO" = sc10x[["factorHTO", drop = TRUE]],
                    "Batch" = sc10x[["orig.ident", drop = TRUE]],
                    cellsCoordinates,
                    jitterCoords = jitter(rep(0, ncol( sc10x)), amount = 0.15) + 0.17); # Precompute jitter for eventual dotplots so coordinates are not recomputed at each refresh

# Create text to show under cursor for each cell for plotly hover
hoverText = do.call(paste, c(Map( paste,
                                  c( "",
                                      "Cell ID: ",
                                      "# UMIs: ",
                                      "# Genes: ",
                                      if("percent.mito" %in% colnames(sc10x[[]])) "% Mito: ",
                                      if("percent.ribo" %in% colnames(sc10x[[]])) "% Ribo: ",
                                      "Cluster: ",
                                      "HTO: "),
                                  cellsData[1:(ncol(cellsData)-3)], # Do not include coordinates in hover text
                                  sep = ""),
                    sep = "\n"));

# Create SharedData object (using plotly accessor) for use by all plots/widgets
dataShare = highlight_key(cellsData);

# Compute the median coordinates for each cluster to plot corresponding name 
dimReducLabels = as.data.frame( do.call( rbind, 
                                         by( cellsData, 
                                             cellsData[["Cluster"]], 
                                             function(x)
                   { 
                     return( data.frame( lapply( x[colnames( cellsCoordinates)], median), 
                                         "Cluster" = x[1, "Cluster"]));
                   })));






##########################
# DIMENSIONALITY REDUCTION
##########################

# ..............................................................................
## @knitr plotDimReducRaster_thumbnail
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster")) +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA,
              size = if(PLOT_DIMREDUC_GROUPS_POINTSIZE==0) Seurat:::AutoPointSize(cellsData) else PLOT_DIMREDUC_GROUPS_POINTSIZE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text_repel( data = dimReducLabels,
                     aes( label = Cluster),
                     color = "black",
                     size = 5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else FALSE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));
print( ggFigure);




# ..............................................................................
## @knitr plotDimReducInteractive_colorClusters
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

#### Dimensionality reduction 2D (plotly)
# Parameters defined in chunk
# dimReducWidth  = 800;
# dimReducHeight = 800;
# initialPointSize = 5;
plotDimReduc = plot_ly( dataShare, # 'SharedData' object for plots interactions
                        x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
                        y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
                        #                       name = ~paste("Cluster", Cluster),
                        color = ~Cluster,
                        colors = clustersColor,
                        #symbol = ~HTO,
                        width = dimReducWidth,
                        height = dimReducHeight) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = initialPointSize,
                            opacity = 0.6),
             text = hoverText,
             hoverinfo = "text+name", 
             showlegend = if(exists( "showLegend")) showLegend else TRUE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  plotDimReduc = plotDimReduc %>% 
    # Use precomputed median of coordinates for each group as new data for Cluster text annotations (type='scatter', mode='text')
    add_trace( data = dimReducLabels,
               x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
               y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
               type = "scattergl",
               mode = "text",
               symbol = NULL,
               text = ~Cluster,
               textfont = list( color = '#000000', size = 20),
               hoverinfo = 'skip',
               showlegend = FALSE)
}

plotDimReduc = plotDimReduc %>% 
  #  hide_colorbar() %>%
  #  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

# Set a unique html ID to access plot with html/javascript for external slider
figID = if(exists("figID")) figID+1 else 0;
plotDimReduc[["elementId"]] = paste0("plotDimReduc_withSlider", figID);

# Render plot
(div(htmltools::tagList(plotDimReduc)));

# Create a html slider using javascript for controlling point size
cat( paste0( '\n',
             '<p>Point size: <span id="valueTag_', figID, '">', initialPointSize, '</span></p>',
             '<input type="range" min="2" max="50" value="', initialPointSize, '" id="slider_', figID, '" style="width: ', dimReducWidth, 'px;">',
             '<script> document.getElementById("slider_', figID, '").oninput = function() { Plotly.restyle("', plotDimReduc[["elementId"]], '", {"marker.size": this.value}); document.getElementById("valueTag_', figID, '").innerHTML = this.value; } </script>',
             '\n'));

# Create a html/js block to change color of clusters and get selected values
cat( paste0( '\n', 
             'Color:&nbsp;',
             # Create a selector, each option corresponding to a plotly trace (defined here by color)
             '<select id="colorChange_traceSelector_plotly', figID, '" style="width: 200px;">', 
             paste('<option value="', c( 'undefined', 0:(nlevels( cellsData[["Cluster"]])-1)), '">', c('All', levels( cellsData[["Cluster"]])), '</option>', sep ="", collapse = " " ), 
             '</select>', 
             '&nbsp;&nbsp;',
             # Create HTML anchor for color picker
             '<input type="text" id="colorChange_colorPicker_plotly', figID, '" />', 
             # Call javascript that transforms it in an actual color picker (updates plotly color of selected trace on change, and resets its own value to 'clear')
             '<script> $("#colorChange_colorPicker_plotly', figID, '").spectrum({ allowEmpty: true, showInput: true, preferredFormat: "rgb", clickoutFiresChange: false, showPalette: true, showSelectionPalette: true, palette: [ ', paste('"', clustersColor, '"', collapse = ", ", sep = ""), ' ], localStorageKey: "myLocalPalette", change: function(color) { if(color) { Plotly.restyle("', plotDimReduc[["elementId"]], '", {"marker.color": color.toHexString()}, $("#colorChange_traceSelector_plotly', figID, '").val()); $("#colorChange_colorPicker_plotly', figID, '").spectrum("set", ""); } } }); </script>', 
             '&nbsp;&nbsp;Get current colors:&nbsp;',
             # Add button to extract current colors from plotly figure (hexadecimal format)
             '<button id="colorGetter_alertButton_plotly', figID, '" onclick="getColorTab(', plotDimReduc[["elementId"]], ', \'scattergl\', \'hex\')">Hexadecimal</button>',
             # Add button to extract current colors from plotly figure (RGB format)
             '<button id="colorGetter_alertButton_plotly', figID, '" onclick="getColorTab(', plotDimReduc[["elementId"]], ', \'scattergl\', \'rgb\')">RGB</button>',
             '\n'));
cat('\n \n')




# ..............................................................................
## @knitr plotDimReducRaster_colorClusters
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster")) +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels,
               aes( label = Cluster),
               color = "black",
               size = 5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));
print( ggFigure);
cat('\n \n')







# ..............................................................................
## @knitr plotDimReducInteractive_colorClusters_cellselection
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

# Subset data used for plotting and update all related variables
#cellselection = sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M4", "M5", "M6")

# For each cell, gather all data to be used for plots generation
cellsData_cellselection = cellsData[ cellselection, ]

# Create text to show under cursor for each cell for plotly hover
hoverText_cellselection = hoverText[cellselection]

# Create SharedData object (using plotly accessor) for use by all plots/widgets
dataShare_cellselection = highlight_key(cellsData_cellselection);



#### Dimensionality reduction 2D (plotly)
# Parameters defined in chunk
# dimReducWidth  = 800;
# dimReducHeight = 800;
# initialPointSize = 5;
plotDimReduc = plot_ly( dataShare_cellselection, # 'SharedData' object for plots interactions
                        x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
                        y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
                        #                       name = ~paste("Cluster", Cluster),
                        color = ~Cluster,
                        colors = clustersColor,
                        #symbol = ~HTO,
                        width = dimReducWidth,
                        height = dimReducHeight) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = initialPointSize,
                            opacity = 0.6),
             text = hoverText_cellselection,
             hoverinfo = "text+name", 
             showlegend = if(exists( "showLegend")) showLegend else TRUE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  plotDimReduc = plotDimReduc %>% 
    # Use precomputed median of coordinates for each group as new data for Cluster text annotations (type='scatter', mode='text')
    add_trace( data = dimReducLabels,
               x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
               y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
               type = "scattergl",
               mode = "text",
               symbol = NULL,
               text = ~Cluster,
               textfont = list( color = '#000000', size = 20),
               hoverinfo = 'skip',
               showlegend = FALSE)
}

plotDimReduc = plotDimReduc %>% 
  #  hide_colorbar() %>%
  #  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

# Set a unique html ID to access plot with html/javascript for external slider
figID = if(exists("figID")) figID+1 else 0;
plotDimReduc[["elementId"]] = paste0("plotDimReduc_withSlider", figID);

# Render plot
(div(htmltools::tagList(plotDimReduc)));

# Create a html slider using javascript for controlling point size
cat( paste0( '\n',
             '<p>Point size: <span id="valueTag_', figID, '">', initialPointSize, '</span></p>',
             '<input type="range" min="2" max="50" value="', initialPointSize, '" id="slider_', figID, '" style="width: ', dimReducWidth, 'px;">',
             '<script> document.getElementById("slider_', figID, '").oninput = function() { Plotly.restyle("', plotDimReduc[["elementId"]], '", {"marker.size": this.value}); document.getElementById("valueTag_', figID, '").innerHTML = this.value; } </script>',
             '\n'));

# Create a html/js block to change color of clusters and get selected values
cat( paste0( '\n', 
             'Color:&nbsp;',
             # Create a selector, each option corresponding to a plotly trace (defined here by color)
             '<select id="colorChange_traceSelector_plotly', figID, '" style="width: 200px;">', 
             paste('<option value="', c( 'undefined', 0:(nlevels( cellsData[["Cluster"]])-1)), '">', c('All', levels( cellsData_cellselection[["Cluster"]])), '</option>', sep ="", collapse = " " ), 
             '</select>', 
             '&nbsp;&nbsp;',
             # Create HTML anchor for color picker
             '<input type="text" id="colorChange_colorPicker_plotly', figID, '" />', 
             # Call javascript that transforms it in an actual color picker (updates plotly color of selected trace on change, and resets its own value to 'clear')
             '<script> $("#colorChange_colorPicker_plotly', figID, '").spectrum({ allowEmpty: true, showInput: true, preferredFormat: "rgb", clickoutFiresChange: false, showPalette: true, showSelectionPalette: true, palette: [ ', paste('"', clustersColor, '"', collapse = ", ", sep = ""), ' ], localStorageKey: "myLocalPalette", change: function(color) { if(color) { Plotly.restyle("', plotDimReduc[["elementId"]], '", {"marker.color": color.toHexString()}, $("#colorChange_traceSelector_plotly', figID, '").val()); $("#colorChange_colorPicker_plotly', figID, '").spectrum("set", ""); } } }); </script>', 
             '&nbsp;&nbsp;Get current colors:&nbsp;',
             # Add button to extract current colors from plotly figure (hexadecimal format)
             '<button id="colorGetter_alertButton_plotly', figID, '" onclick="getColorTab(', plotDimReduc[["elementId"]], ', \'scattergl\', \'hex\')">Hexadecimal</button>',
             # Add button to extract current colors from plotly figure (RGB format)
             '<button id="colorGetter_alertButton_plotly', figID, '" onclick="getColorTab(', plotDimReduc[["elementId"]], ', \'scattergl\', \'rgb\')">RGB</button>',
             '\n'));
cat('\n \n')




# ..............................................................................
## @knitr plotDimReducRaster_colorClusters_cellselection
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

# Subset data used for plotting and update all related variables
#cellselection = sc10x[["HTOMICE_classification", drop = TRUE]] %in% c("M4", "M5", "M6")

# For each cell, gather all data to be used for plots generation
cellsData_cellselection = cellsData[ cellselection, ]

# Create text to show under cursor for each cell for plotly hover
hoverText_cellselection = hoverText[cellselection]

# Create SharedData object (using plotly accessor) for use by all plots/widgets
dataShare_cellselection = highlight_key(cellsData_cellselection);


# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData_cellselection, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster")) +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels,
               aes( label = Cluster),
               color = "black",
               size = 5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));
print( ggFigure);
cat('\n \n')







# ..............................................................................
## @knitr plotDimReducInteractive_colorClusters_symbolHTOs
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

#### Dimensionality reduction 2D (plotly)
# Parameters defined in chunk
# dimReducWidth  = 800;
# dimReducHeight = 800;
# initialPointSize = 5;

plotDimReduc = plot_ly( dataShare, # 'SharedData' object for plots interactions
                        x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
                        y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
                        #                       name = ~paste("Cluster", Cluster),
                        color = ~Cluster,
                        colors = clustersColor,
                        symbol = ~HTO,
                        width = dimReducWidth,
                        height = dimReducHeight) %>%
  add_trace( type = "scattergl",
             mode = "markers",
             marker = list( size = initialPointSize,
                            opacity = 0.6),
             text = hoverText,
             hoverinfo = "text+name", 
             showlegend = if(exists( "showLegend")) showLegend else TRUE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  plotDimReduc = plotDimReduc %>% 
    # Use precomputed median of coordinates for each group as new data for Cluster text annotations (type='scatter', mode='text')
    add_trace( data = dimReducLabels,
               x = as.formula( paste( "~", colnames( cellsCoordinates)[1])),
               y = as.formula( paste( "~", colnames( cellsCoordinates)[2])),
               type = "scattergl",
               mode = "text",
               symbol = NULL,
               text = ~Cluster,
               textfont = list( color = '#000000', size = 20),
               hoverinfo = 'skip',
               showlegend = FALSE)
}

plotDimReduc = plotDimReduc %>% 
  #  hide_colorbar() %>%
  #  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

# Set a unique html ID to access plot with html/javascript for external slider
figID = if(exists("figID")) figID+1 else 0;
plotDimReduc[["elementId"]] = paste0("plotDimReduc_withSlider", figID);

# Render plot
(div(htmltools::tagList(plotDimReduc)));

# Create a html slider using javascript for controlling point size
cat( paste0( '\n',
             '<p>Point size: <span id="valueTag_', figID, '">', initialPointSize, '</span></p>',
             '<input type="range" min="2" max="50" value="', initialPointSize, '" id="slider_', figID, '" style="width: ', dimReducWidth, 'px;">',
             '<script> document.getElementById("slider_', figID, '").oninput = function() { Plotly.restyle("', plotDimReduc[["elementId"]], '", {"marker.size": this.value}); document.getElementById("valueTag_', figID, '").innerHTML = this.value; } </script>',
             '\n'));
cat('\n \n')




# ..............................................................................
## @knitr plotDimReducRaster_colorClusters_symbolHTOs
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster",
                                pch = "HTO")) +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE)

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels,
               aes( label = Cluster, pch = NULL),
               color = "black",
               size = 3.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));
print( ggFigure);
cat('\n \n')




# ..............................................................................
## @knitr plotDimReducInteractiveGG_colorClusters_facetHTOs
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

# Using 'ggplotly' as facetting is not trivial with plotly 

#### Facetted dimensionality reduction 2D (ggplot + plotly conversion)
# Parameters defined in chunk
# dimReducWidth  = 800;
# dimReducHeight = 800;
# initialPointSize = 5;

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster")) +
  geom_point( stroke = 0) + # Just ignore point size and alpha here as it will be overriden after plotly conversion
  facet_wrap(facets = vars(HTO))

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels, # 'repel' not implemented by ggplotly yet
               aes( label = Cluster),
               color = "black",
               size = 3.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));

plotDimReducGG =  ggplotly( ggFigure, 
                            width = dimReducWidth, 
                            height = dimReducHeight) %>% 
  style( type = "scattergl",   # trick to use webGL for rendering the scatterplot since toWebGL does not have width/height arguments, and use of width/height in layout is deprecated
         marker.size = initialPointSize, # Set point size and opacity with plotly so we can synchronize them with slider (ggplot and plotly units don't match)
         marker.opacity = 0.6) %>% 
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

# Set a unique html ID to access plot with html/javascript for external slider
figID = if(exists("figID")) figID+1 else 0;
plotDimReducGG[["elementId"]] = paste0("plotDimReducGG_withSlider", figID);

# Render plot. ggplot conversion generates warnings about incompatibility between 'scatter' and 'hoveron'
div(htmltools::tagList(plotDimReducGG));

# Create a html slider using javascript for controlling point size
cat( paste0( '\n',
             '<p>Point size: <span id="valueTag_', figID, '">', initialPointSize, '</span></p>',
             '<input type="range" min="2" max="50" value="', initialPointSize, '" id="slider_', figID, '" style="width: ', dimReducWidth, 'px;">',
             '<script> document.getElementById("slider_', figID, '").oninput = function() { Plotly.restyle("', plotDimReducGG[["elementId"]], '", {"marker.size": this.value}); document.getElementById("valueTag_', figID, '").innerHTML = this.value; } </script>',
             '\n'));
cat('\n \n')




# ..............................................................................
## @knitr plotDimReducRaster_colorClusters_facetHTOs
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster"))  +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE) +
  facet_wrap(facets = vars(HTO))

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels,
               aes( label = Cluster),
               color = "black",
               size = 2.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));
print( ggFigure);
cat('\n \n')





# ..............................................................................
## @knitr plotDimReducInteractiveGG_colorHTOs_facetClusters
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................
# Using 'ggplotly' as facetting is not trivial with plotly 

#### Facetted dimensionality reduction 2D (ggplot + plotly conversion)
# Parameters defined in chunk
# dimReducWidth  = 800;
# dimReducHeight = 800;
# initialPointSize = 5;

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "HTO")) +
  geom_point( stroke = 0) + # Just ignore point size and alpha here as it will be overriden after plotly conversion
  facet_wrap(facets = vars(Cluster), scales = "free") # Wrap and rescale around cluster coordinates

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels, # 'repel' not implemented by ggplotly yet
               aes( label = Cluster),
               color = "black",
               size = 3.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "HTO", values = HTOsColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));

plotDimReducGG =  ggplotly( ggFigure, 
                            width = dimReducWidth, 
                            height = dimReducHeight) %>% 
  style( type = "scattergl",   # trick to use webGL for rendering the scatterplot since toWebGL does not have width/height arguments, and use of width/height in layout is deprecated
         marker.size = initialPointSize, # Set point size and opacity with plotly so we can synchronize them with slider (ggplot and plotly units don't match)
         marker.opacity = 0.6) %>% 
  config( displaylogo = FALSE,
          toImageButtonOptions = list( format='svg'),
          modeBarButtons = list(list( 'toImage'),
                                list( 'zoom2d', 
                                      'zoomIn2d', 
                                      'zoomOut2d', 
                                      'pan2d', 
                                      'resetScale2d')))

# Set a unique html ID to access plot with html/javascript for external slider
figID = if(exists("figID")) figID+1 else 0;
plotDimReducGG[["elementId"]] = paste0("plotDimReducGG_withSlider", figID);

# Render plot
div(htmltools::tagList(plotDimReducGG));

# Create a html slider using javascript for controlling point size
cat( paste0( '\n',
             '<p>Point size: <span id="valueTag_', figID, '">', initialPointSize, '</span></p>',
             '<input type="range" min="2" max="50" value="', initialPointSize, '" id="slider_', figID, '" style="width: ', dimReducWidth, 'px;">',
             '<script> document.getElementById("slider_', figID, '").oninput = function() { Plotly.restyle("', plotDimReducGG[["elementId"]], '", {"marker.size": this.value}); document.getElementById("valueTag_', figID, '").innerHTML = this.value; } </script>',
             '\n'));
cat('\n \n')




# ..............................................................................
## @knitr plotDimReducRaster_colorHTOs_facetClusters
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "HTO"))  +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE) +
  facet_wrap(facets = vars(Cluster), scales = "free") # Wrap and rescale around cluster coordinates

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels,
               aes( label = Cluster),
               color = "black",
               size = 2.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "HTO", values = HTOsColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"),
    plot.title = element_text( face = "bold",
                               size = rel( 1.8), #rel( 16/14),
                               hjust = 0.5,
                               vjust = 1,
                               margin = margin( b = 7)));
print( ggFigure);
cat('\n \n')




# ..............................................................................
## @knitr plotDimReducRaster_colorHTOs_facetClusters_allInBackground
# Can set logical for 'showDimReducClusterLabels' and 'showLegend' variables
# ..............................................................................

# Prepare dimreduc figure as ggplot
ggFigure = ggplot(  cellsData[c(colnames( cellsCoordinates)[1], colnames( cellsCoordinates)[2], "HTO")], # Omit the column used for facetting to get all points repeated in all facets
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2])) +
  geom_point( alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE/2,
              color = "grey") +
  geom_point( data = cellsData, # Now provide data including column for facetting
              aes(color = HTO),
              alpha = PLOT_DIMREDUC_GROUPS_ALPHA, 
              size  = PLOT_DIMREDUC_GROUPS_POINTSIZE) +
  facet_wrap(facets = vars(Cluster)) 

if(if(exists( "showDimReducClusterLabels")) showDimReducClusterLabels else TRUE)
{
  ggFigure = ggFigure +
    geom_text( data = dimReducLabels,
               aes( label = Cluster),
               color = "black",
               size = 2.5)
}

ggFigure = ggFigure +
  scale_color_manual( name = "HTO", values = HTOsColor) +
  theme_classic() +
  theme( #axis.text  = element_blank(),
    axis.title = element_blank(),
    #axis.ticks  = element_blank(),
    #axis.line = element_blank(),
    legend.position = if(if(exists( "showLegend")) showLegend else TRUE) NULL else "none",
    plot.margin = margin( 0, 0, 0, 0, "cm"));
print( ggFigure);
cat('\n \n')




###############
# DENSITY PLOTS
###############

# ..............................................................................
## @knitr plotDimReducRaster_colorDensity_circleDensity_prepareData
# ..............................................................................

# Prepare data for plotting density plots (cluster groups are highlighted by 2D
# density estimated on global population).

# ggFigure = ggplot(  cellsData, 
#                     aes_string( x = colnames( cellsCoordinates)[1], 
#                                 y = colnames( cellsCoordinates)[2]))  +
#   geom_density_2d_filled( aes( alpha = stat(level)),
#                           contour_var = "ndensity",
#                           n = 100 ,
#                           bins = nbBins,
#                           adjust = 1/2,
#                           show.legend = FALSE) +
#   scale_alpha_manual(values = c(rep(0, nbHideBins), seq(0.75,1,length.out = nbBins-nbHideBins))) +
#   scale_fill_manual(values = rainbow( nbBins, end = 0.75, rev = TRUE)) +
#   #geom_point() +
#   facet_grid(rows = vars( ifelse( grepl( "Panth",HTO), "Panth", "PBS")), # vars( HTO)
#              cols = vars( Cluster),
#              scales = "free",
#              space = "fixed")
# 
# ggFigure

nbBins = 100
nbHideBins = 10 # Number of 'low' bins to exclude from representation

# Prepare a copy of actual data to be reshaped for graphics (extend range for 
# density, then repeat coords for all HTOs to show on all facets)
cellsData_subset = cellsData[!grepl("negative", cellsData[["HTO"]], ignore.case = TRUE) ,c(colnames( cellsCoordinates), "Cluster", "HTO")]

# Add factors of interest for facetting
cellsData_subset[["Day"]] =  ifelse( gsub( ".*_", "", cellsData_subset[["HTO"]]) %in% c("M1", "M2", "M3"), "D0", "D7")
cellsData_subset[["Tissue"]] =  gsub( "_.*", "", cellsData_subset[["HTO"]])
cellsData_subset[["Mouse"]] = gsub(".*_", "", cellsData_subset[["HTO"]]) 

# Add flag whether each point is a bounding one (for debug, should not affect
# density computation as long as it is far enough and breakpoint is high enough)
cellsData_subset[["isBoundingPoint"]] = FALSE

# Artificially extend range around each combination of points with bounding 
# points around coordinates so 2D density covers the full plot (mul)
extendCellsCoordRange = function(x, clustNA = FALSE, extendRangeMultiplier = 0.45)
                        {
                          # Compute coordinates of bounding points
                          boundingPoints = data.frame( expand_range( range(x[[1]]), mul = extendRangeMultiplier), # Coord x
                                                       expand_range( range(x[[2]]), mul = extendRangeMultiplier), # Coord y
                                                       if(clustNA) NA else x[[3]][1], # Cluster (clustNA useful when applying function on all cells, we don't want to attribute a random cluster to bounding points, which interferes with plotting)
                                                       x[[4]][1], # HTO
                                                       x[[5]][1], # Day
                                                       x[[6]][1], # Tissue
                                                       x[[7]][1], # Mouse
                                                       TRUE)      # isBoundingPoint
                          
                          colnames(boundingPoints) = colnames(x)
                          # Add to existing data
                          rbind(x, boundingPoints)
}

# Extend clusters range (adds pair of bounding points for each cluster)
cellsData_extendRangeByCluster = do.call( rbind, 
                                          by( cellsData_subset, 
                                              cellsData_subset[c("Cluster")],
                                              extendCellsCoordRange))



# ..............................................................................
## @knitr plotDimReducRaster_colorDensity_circleDensity_globalRefGroups
## Plot 1 - Global figure with reference groups circles
# ..............................................................................
ggFigure = ggplot(  cellsData_subset, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2],
                                color = "Cluster"))  +
  geom_point(alpha = 0.35)  +
  # Highlight single-level density estimation (breaks), based on extended coords
  # set of points (large white line behind colored line).
  # Extend coordinates based on all points as we don't focus facets on cluster.
  geom_density_2d( data = extendCellsCoordRange(cellsData_subset, clustNA = TRUE, extendRangeMultiplier = 0.15), # Setting "clustNA = TRUE" prevents assigning an identity to bounding points as we extend independently of clustering
                   aes( group = Cluster, # Mimic the 'color' grouping 
                        linetype = isBoundingPoint), # Make sure bounding points are considered as separate group for computing density and plotting, so never mixed with actual data points
                   color = "white",
                   contour_var = "ndensity",
                   breaks = c( .05), # Unique break to get a simple circling
                   adjust = 1.5,
                   lwd = 2,
                   n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                   show.legend = FALSE) +
  geom_density_2d( data = extendCellsCoordRange(cellsData_subset, clustNA = TRUE, extendRangeMultiplier = 0.15), # "clustNA = TRUE" prevents assigning an identity to bounding points as we extend independently of clustering.
                   aes( color = Cluster,
                        linetype = isBoundingPoint), # Make sure bounding points are considered as separate group for computing density and plotting, so never mixed with actual data points
                   contour_var = "ndensity",
                   breaks = c( .05), # Unique break to get a simple circling
                   adjust = 1.5,
                   lwd = .75,
                   n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                   show.legend = FALSE) +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  scale_linetype_manual(values = c("TRUE" = "blank", "FALSE" = "solid")) # Make sure that even if density in bounding points group is enough to plot a line (should not), we would hide it

print( ggFigure);
cat('\n \n')




# ..............................................................................
## @knitr plotDimReducRaster_colorDensity_circleDensity_facetByHTO
### Plot 2 - Density facetted by HTO (Mice + Tissue)
# ..............................................................................
# ggFigure = ggplot(  cellsData_subset, 
#                     aes_string( x = colnames( cellsCoordinates)[1], 
#                                 y = colnames( cellsCoordinates)[2]))  +
#   geom_density_2d_filled( aes( alpha = stat(level)),
#                           contour_var = "ndensity",
#                           bins = nbBins,
#                           adjust = .25,
#                           n = 500, # Nb grid points to draw lines, prevents angles/aliasing
#                           show.legend = FALSE) +
#   # Define a manual color scale for density (+ make lower bins transparent)
#   scale_alpha_manual(values = c(rep(0, nbHideBins), seq(0.75,1,length.out = nbBins-nbHideBins))) +
#   scale_fill_manual(values = rainbow( nbBins, end = 0.75, rev = TRUE)) +
#   #geom_point() +
#   #facet_wrap( vars( factor(HTO, levels = HTO_FACTOR_LEVELS))) +
#   facet_wrap( vars( factor(HTO))) +
#   # Highlight single-level density estimation (breaks), based on extended coords
#   # set of points.
#   # Extend coordinates based on all points as we don't focus facets on cluster.
#   geom_density_2d( data = extendCellsCoordRange(cellsData_subset, clustNA = TRUE, extendRangeMultiplier = 0.15)[,-4], # Remove column HTO to repeat all data in all facets. Setting "clustNA = TRUE" prevents assigning an identity to bounding points as we extend independently of clustering.
#                    aes( color = Cluster,
#                         linetype = isBoundingPoint), # Make sure bounding points are considered as separate group for computing density and plotting, so never mixed with actual data points
#                    contour_var = "ndensity",
#                    breaks = c( .05), # Unique break to get a simple circling
#                    adjust = 1.5,
#                    lwd = .5,
#                    n = 500, # Nb grid points to draw lines, prevents angles/aliasing
#                    show.legend = FALSE) +
#   scale_linetype_manual(values = c("TRUE" = "blank", "FALSE" = "dashed")) # Make sure that even if density in bounding points group is enough to plot a line (should not), we would hide it
# 
# print( ggFigure);
# cat('\n \n')




# ..............................................................................
## @knitr plotDimReducRaster_colorDensity_circleDensity_facetByTissue
## Plot 2 - Density facetted by Tissue
# ..............................................................................
ggFigure = ggplot(  cellsData_subset, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2]))  +
  geom_density_2d_filled( aes( alpha = stat(level)),
                          contour_var = "ndensity",
                          bins = nbBins,
                          adjust = .2,
                          n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                          show.legend = FALSE) +
  # Define a manual color scale for density (+ make lower bins transparent)
  scale_alpha_manual(values = c(rep(0, nbHideBins), seq(0.75,1,length.out = nbBins-nbHideBins))) +
  scale_fill_manual(values = rainbow( nbBins, end = 0.75, rev = TRUE)) +
  #geom_point() +
  #facet_wrap( vars( factor(HTO, levels = HTO_FACTOR_LEVELS))) +
  facet_wrap( vars( factor(Tissue))) +
  # Highlight single-level density estimation (breaks), based on extended coords
  # set of points.
  # Extend coordinates based on all points as we don't focus facets on cluster.
  geom_density_2d( data = extendCellsCoordRange(cellsData_subset, clustNA = TRUE, extendRangeMultiplier = 0.1)[,-c(4:7)], # Remove column Tissue to repeat all data in all facets. Setting "clustNA = TRUE" prevents assigning an identity to bounding points as we extend independently of clustering.
                   aes( color = Cluster,
                        linetype = isBoundingPoint), # Make sure bounding points are considered as separate group for computing density and plotting, so never mixed with actual data points
                   contour_var = "ndensity",
                   breaks = c( .05), # Unique break to get a simple circling
                   adjust = 1.5,
                   lwd = .5,
                   n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                   show.legend = FALSE) +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  scale_linetype_manual(values = c("TRUE" = "blank", "FALSE" = "dashed")) # Make sure that even if density in bounding points group is enough to plot a line (should not), we would hide it

print( ggFigure);
cat('\n \n')


# 
# # ..............................................................................
# ## @knitr plotDimReducRaster_colorDensity_circleDensity_facetByHTOAndCluster
# ## Plot 3 - Density facetted by HTO & Cluster
# # ..............................................................................
# ggFigure = ggplot(  cellsData_subset, 
#                     aes_string( x = colnames( cellsCoordinates)[1], 
#                                 y = colnames( cellsCoordinates)[2]))  +
#   geom_density_2d_filled( aes( alpha = stat(level)),
#                           contour_var = "ndensity",
#                           bins = nbBins,
#                           adjust = 1,
#                           n = 500, # Nb grid points to draw lines, prevents angles/aliasing
#                           show.legend = FALSE) +
#   # Define a manual color scale for density (+ make lower bins transparent)
#   scale_alpha_manual(values = c(rep(0, nbHideBins), seq(0.75,1,length.out = nbBins-nbHideBins))) +
#   scale_fill_manual(values = rainbow( nbBins, end = 0.75, rev = TRUE)) +
#   #geom_point() +
#   # Use 'facet_grid2' to allow independent y scale on grid
#   #facet_grid2( rows = vars( factor(HTO, levels = HTO_FACTOR_LEVELS)), #vars( ifelse( grepl( "Panth", HTO), "Panth", "PBS")),
#   facet_grid2( rows = vars( factor(HTO)), #vars( ifelse( grepl( "Panth", HTO), "Panth", "PBS")),
#                cols = vars( Cluster),
#                scales = "free",
#                space = "fixed",
#                independent = "y") + # Specific to facet_grid2 (package ggh4x)
#   # Highlight single-level density estimation (breaks), based on extended coords
#   # set of points for each cluster (we show facets on each cluster coords here).
#   geom_density_2d( data = cellsData_extendRangeByCluster[-4], # Remove column HTO to repeat 'Cluster' data in all HTOs (rows)
#                    aes( color = Cluster,
#                         linetype = isBoundingPoint), # Make sure bounding points are considered as separate group for computing density and plotting, so never mixed with actual data points
#                    contour_var = "ndensity",
#                    breaks = c( .05), # Unique break to get a simple circling
#                    adjust = 1.5,
#                    alpha = 1,
#                    lwd = .5,
#                    n = 500, # Nb grid points to draw lines, prevents angles/aliasing
#                    show.legend = FALSE) +
#   scale_linetype_manual(values = c("TRUE" = "blank", "FALSE" = "dashed")) # Make sure that even if density in bounding points group is enough to plot a line (should not), we would hide it
# 
# print( ggFigure);
# cat('\n \n')
# 




# ..............................................................................
## @knitr plotDimReducRaster_colorDensity_circleDensity_facetByTissueAndDay
## Plot 3 - Density facetted by Tissue and Day
# ..............................................................................

ggFigure = ggplot(  cellsData_subset, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2]))  +
  geom_density_2d_filled( aes( alpha = stat(level)),
                          contour_var = "ndensity",
                          bins = nbBins,
                          adjust = 0.15,
                          n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                          show.legend = FALSE) +
  # Define a manual color scale for density (+ make lower bins transparent)
  scale_alpha_manual(values = c(rep(0, nbHideBins), seq(0.75,1,length.out = nbBins-nbHideBins))) +
  scale_fill_manual(values = rainbow( nbBins, end = 0.75, rev = TRUE)) +
  #geom_point() +
  # Use 'facet_grid2' to allow independent y scale on grid
  facet_grid2( rows = vars( Day), # Faceting on computed variable (merged timepoints)
               cols = vars( Tissue),
               scales = "fixed",
               space = "fixed") +
  #independent = "y") + # Specific to facet_grid2 (package ggh4x)
  # Highlight single-level density estimation (breaks), based on extended coords
  # set of points for each cluster (we show facets on each cluster coords here).
  geom_density_2d( data = extendCellsCoordRange(cellsData_subset, clustNA = TRUE, extendRangeMultiplier = 0.1)[,-c(4:7)], # Remove any facetting columns to repeat cpoints between conditions and compute global proper cluster limits  
                   aes( color = Cluster,
                        linetype = isBoundingPoint), # Make sure bounding points are considered as separate group for computing density and plotting, so never mixed with actual data points
                   contour_var = "ndensity",
                   breaks = c( .05),
                   adjust = 1.5,
                   alpha = 1,
                   lwd = .5,
                   n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                   show.legend = FALSE) +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  scale_linetype_manual(values = c("TRUE" = "blank", "FALSE" = "dashed")) # Make sure that even if density in bounding points group is enough to plot a line (should not), we would hide it

print( ggFigure);
cat('\n \n')




# ..............................................................................
## @knitr plotDimReducRaster_colorDensity_circleDensity_facetByTissue_M123
## Plot 4 - Density facetted by Tissue (Mice 123 only = D0)
# ..............................................................................

# Subselect Mice 123 (D0) and compute correspoding extended coordinates
cellsData_subset_M123 = cellsData_subset[cellsData_subset[["Day"]] %in% "D0", ]

# Extend clusters range (adds pair of bounding points for each cluster)
cellsData_extendRangeByCluster_M123 = do.call( rbind, 
                                               by( cellsData_subset_M123, 
                                                   cellsData_subset_M123[c("Cluster")],
                                                   extendCellsCoordRange,
                                                   extendRangeMultiplier = 0.3)) # Must increase to accomodate for less cells


ggFigure = ggplot(  cellsData_subset_M123, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2]))  +
  geom_density_2d_filled( aes( alpha = stat(level)),
                          contour_var = "ndensity",
                          bins = nbBins,
                          adjust = 0.1,
                          n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                          show.legend = FALSE) +
  # Define a manual color scale for density (+ make lower bins transparent)
  scale_alpha_manual(values = c(rep(0, nbHideBins), seq(0.75,1,length.out = nbBins-nbHideBins))) +
  scale_fill_manual(values = rainbow( nbBins, end = 0.75, rev = TRUE)) +
  #geom_point() +
  # Use 'facet_grid2' to allow independent y scale on grid
  facet_grid2( rows = vars( Day), # Faceting on computed variable (merged timepoints)
               cols = vars( Tissue),
               scales = "fixed",
               space = "fixed") +
  #independent = "y") + # Specific to facet_grid2 (package ggh4x)
  # Highlight single-level density estimation (breaks), based on extended coords
  # set of points for each cluster (we show facets on each cluster coords here).
  geom_density_2d( data = cellsData_extendRangeByCluster_M123[-c(4:7)], # Remove any facetting columns to repeat cpoints between conditions and compute global proper cluster limits  
                   aes( color = Cluster,
                        linetype = isBoundingPoint), # Make sure bounding points are considered as separate group for computing density and plotting, so never mixed with actual data points
                   contour_var = "ndensity",
                   breaks = c( .05),
                   adjust = 1.5,
                   alpha = 1,
                   lwd = .5,
                   n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                   show.legend = FALSE) +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  scale_linetype_manual(values = c("TRUE" = "blank", "FALSE" = "dashed")) # Make sure that even if density in bounding points group is enough to plot a line (should not), we would hide it

print( ggFigure);
cat('\n \n')







# ..............................................................................
## @knitr plotDimReducRaster_colorDensity_circleDensity_facetByTissue_M456
## Plot 5 - Density facetted by Tissue (Mice 456 only = D7)
# ..............................................................................

# Subselect Mice 456 (D7) and compute correspoding extended coordinates
cellsData_subset_M456 = cellsData_subset[cellsData_subset[["Day"]] %in% "D7", ]

# Extend clusters range (adds pair of bounding points for each cluster)
cellsData_extendRangeByCluster_M456 = do.call( rbind, 
                                               by( cellsData_subset_M456, 
                                                   cellsData_subset_M456[c("Cluster")],
                                                   extendCellsCoordRange,
                                                   extendRangeMultiplier = 0.1))


ggFigure = ggplot(  cellsData_subset_M456, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2]))  +
  geom_density_2d_filled( aes( alpha = stat(level)),
                          contour_var = "ndensity",
                          bins = nbBins,
                          adjust = 0.2,
                          n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                          show.legend = FALSE) +
  # Define a manual color scale for density (+ make lower bins transparent)
  scale_alpha_manual(values = c(rep(0, nbHideBins), seq(0.75,1,length.out = nbBins-nbHideBins))) +
  scale_fill_manual(values = rainbow( nbBins, end = 0.75, rev = TRUE)) +
  #geom_point() +
  # Use 'facet_grid2' to allow independent y scale on grid
  facet_grid2( rows = vars( Day), # Faceting on computed variable (merged timepoints)
               cols = vars( Tissue),
               scales = "fixed",
               space = "fixed") +
  #independent = "y") + # Specific to facet_grid2 (package ggh4x)
  # Highlight single-level density estimation (breaks), based on extended coords
  # set of points for each cluster (we show facets on each cluster coords here).
  geom_density_2d( data = cellsData_extendRangeByCluster_M456[-c(4:7)], # Remove any facetting columns to repeat cpoints between conditions and compute global proper cluster limits  
                   aes( color = Cluster,
                        linetype = isBoundingPoint), # Make sure bounding points are considered as separate group for computing density and plotting, so never mixed with actual data points
                   contour_var = "ndensity",
                   breaks = c( .05),
                   adjust = 1.5,
                   alpha = 1,
                   lwd = .5,
                   n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                   show.legend = FALSE) +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  scale_linetype_manual(values = c("TRUE" = "blank", "FALSE" = "dashed")) # Make sure that even if density in bounding points group is enough to plot a line (should not), we would hide it

print( ggFigure);
cat('\n \n')





# ..............................................................................
## @knitr plotDimReducRaster_colorDensity_circleDensity_facetByMiceAndTissue
### Plot 6 - Density facetted by Mice and Tissue
# ..............................................................................


ggFigure = ggplot(  cellsData_subset, 
                    aes_string( x = colnames( cellsCoordinates)[1], 
                                y = colnames( cellsCoordinates)[2]))  +
  geom_density_2d_filled( aes( alpha = stat(level)),
                          contour_var = "ndensity",
                          bins = nbBins,
                          adjust = 0.2,
                          n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                          show.legend = FALSE) +
  # Define a manual color scale for density (+ make lower bins transparent)
  scale_alpha_manual(values = c(rep(0, nbHideBins), seq(0.75,1,length.out = nbBins-nbHideBins))) +
  scale_fill_manual(values = rainbow( nbBins, end = 0.75, rev = TRUE)) +
  #geom_point() +
  # Use 'facet_grid2' to allow independent y scale on grid
  facet_grid2( rows = vars( Mouse),
               cols = vars( Tissue),
               scales = "fixed",
               space = "fixed") +
  #independent = "y") + # Specific to facet_grid2 (package ggh4x)
  # Highlight single-level density estimation (breaks), based on extended coords
  # set of points for each cluster (we show facets on each cluster coords here).
  geom_density_2d( data = extendCellsCoordRange(cellsData_subset, clustNA = TRUE, extendRangeMultiplier = 0.1)[,-c(4:7)], # Remove any facetting columns to repeat cpoints between conditions and compute global proper cluster limits  
                   aes( color = Cluster,
                        linetype = isBoundingPoint), # Make sure bounding points are considered as separate group for computing density and plotting, so never mixed with actual data points
                   contour_var = "ndensity",
                   breaks = c( .05),
                   adjust = 1.5,
                   alpha = 1,
                   lwd = .5,
                   n = 500, # Nb grid points to draw lines, prevents angles/aliasing
                   show.legend = FALSE) +
  scale_color_manual( name = "Cluster", values = clustersColor) +
  scale_linetype_manual(values = c("TRUE" = "blank", "FALSE" = "dashed")) # Make sure that even if density in bounding points group is enough to plot a line (should not), we would hide it

print( ggFigure);
cat('\n \n')


