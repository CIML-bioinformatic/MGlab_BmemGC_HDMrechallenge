#!/bin/bash

# This script execute the CellRanger count analysis on MGlab data
# with mode multi to also assign HTO counts and BCR counts.
# It requires following informations:
# CONFIG_CSV_FILE, OUTPUT_FOLDER
# Path to reference and cellranger command (--id=) must fit correct genome (currently mm10)

# Path to the file describing libraries and other description files (paths to fastq files and sample/HTO names)
CONFIG_CSV_FILE="/mnt/DOSI/MGLAB/BIOINFO/Project/AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge/230321_VH00228_159_AACKKJ7M5/01_Reference/01_CellRanger/config_multi.csv"

# Path to the CellRanger singularity image
SINGULARITY_IMAGE="/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/02_CONTAINER/02_SINGULARITY/CellRanger/7_0_1/cellranger_701.sif"

# Output folder for the CellRanger output files
OUTPUT_FOLDER="/mnt/DOSI/MGLAB/BIOINFO/Project/AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge/230321_VH00228_159_AACKKJ7M5/05_Output/01_CellRanger_FeatureBarcoding"

mkdir -p $OUTPUT_FOLDER

#### Execution of the analysis
echo CONFIG_CSV_FILE = $CONFIG_CSV_FILE 
echo OUTPUT_FOLDER = $OUTPUT_FOLDER

cd $OUTPUT_FOLDER

echo -e "\n\nCommand:\nsingularity exec -B /mnt/DOSI:/mnt/DOSI $SINGULARITY_IMAGE cellranger multi --id=mm10 --csv=$CONFIG_CSV_FILE"

singularity exec -B /mnt/DOSI:/mnt/DOSI $SINGULARITY_IMAGE cellranger multi --id=mm10 --csv=$CONFIG_CSV_FILE





