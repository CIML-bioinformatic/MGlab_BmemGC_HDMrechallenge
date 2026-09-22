#!/bin/bash

# Execute 'velocyto' command line program to compute spliced/unspliced matrices


# Set experiment base path
export EXPERIMENT_PATH="/mnt/DOSI/MGLAB/BIOINFO/Project/AIDYFP_BmemGC_PC_Day0and7_HDMrechallenge/230321_VH00228_159_AACKKJ7M5"

# Provide the path to the genome annotation file
export GENOME_ANNOTATION_PATH="/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/01_REFERENCE/01_GENOME/CellRanger/mouse/2020-A/refdata-gex-mm10-2020-A/genes/genes.gtf"

# Path to the BAM file
export BAM_PATH="${EXPERIMENT_PATH}/05_Output/01_CellRanger_FeatureBarcoding/mm10/outs/per_sample_outs/mm10/count/sample_alignments.bam"

# Path to barcode file corresponding to BAM
export BARCODE_PATH="${EXPERIMENT_PATH}/05_Output/01_CellRanger_FeatureBarcoding/mm10/outs/per_sample_outs/mm10/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"


# Create the folder where to put the results (loom file)
export OUTPUT_PATH="${EXPERIMENT_PATH}/05_Output/08_velocyto"
mkdir -p $OUTPUT_PATH 


# Execute Velocyto on 10x data
cd $OUTPUT_PATH

# Execute velocyto for each
velocyto run -o "$OUTPUT_PATH" -b "${BARCODE_PATH}" --sampleid "230321_VH00228_159_AACKKJ7M5" "${BAM_PATH}" "$GENOME_ANNOTATION_PATH"

 
