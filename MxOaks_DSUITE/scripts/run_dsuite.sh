#!/bin/bash

# Activate environment if not already
# conda activate vcftools

# Directory containing the SETS files
SETS_DIR="data/dsuite_sets/"

# Directory for output files
OUTPUT_DIR="out/20240913_dsuite_output/"

# Path to your tree file (nwk)
TREE_FILE="data/dsuite_white_redout.nwk"

# Path to your VCF file
VCF_FILE="data/20240915_whiteoak_redout_renamed.vcf"

# Array of prefixes
PREFIXES=("SETS_white")

# Create output directory if it doesn't exist
# mkdir -p "$OUTPUT_DIR"

# Loop through each prefix
for PREFIX in "${PREFIXES[@]}"; do
    echo "Processing $PREFIX files..."
    
    # Loop through each file with the current prefix
    for SETS_FILE in ${SETS_DIR}/${PREFIX}*; do
        # Get the base name of the SETS file (without path and extension)
        BASE_NAME=$(basename "$SETS_FILE" | sed 's/\.[^.]*$//')
        # BASE_NAME=$(basename "$SETS_FILE")
        
        echo "Processing file: $BASE_NAME"
        
        # Run Dsuite Dtrios
        ../Dsuite/Build/Dsuite Dtrios -c -o "$OUTPUT_DIR" -n "$BASE_NAME" -t "$TREE_FILE" "$VCF_FILE" "$SETS_FILE"
        
        echo "Completed processing $BASE_NAME"
        echo
    done
    
    echo "Completed processing all $PREFIX files."
    echo
done

echo "All processing complete."