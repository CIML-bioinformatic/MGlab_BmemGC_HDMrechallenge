###############################################################################
# This file defines SAMPLE parameters as global variables that will be loaded
# before analysis starts. 
#


SAMPLE_NAME = "RUN1"

SAMPLE_COLOR = c( "#00BA38")
names( SAMPLE_COLOR) = SAMPLE_NAME

MOUSE_COLORS = c("M1" = "#E90000", "M2" = "#A40D03", "M3" = "#7D0000", "M4" = "#4ACADE", "M5" = "#38A6D9", "M6" = "#4076D9")


library(RColorBrewer)
TISSUE_COLORS = brewer.pal(name = "Pastel2", n=3)
names(TISSUE_COLORS) = c("LN", "LG", "LN+LG")
