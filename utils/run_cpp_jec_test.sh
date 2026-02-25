#!/bin/bash
# run_cpp_jec_test.sh
# 
# Run C++ JEC map creator on all QCD HT-binned datasets
# Processes all files at once to generate combined correction maps

set -e  # Exit on error

echo "=============================================================="
echo "Scouting JEC Map Creator - Production Run (All HT bins)"
echo "=============================================================="

# Configuration
BASE_PATH="/store/user/mgaisdor/SVJScouting_ntuples/MC_offline/2017"
MAX_FILES_PER_BIN=500  # Set to -1 to use all available files
OUTPUT_DIR="test"
SERVER="cmsdcache-kit-disk.gridka.de:1094"

# QCD HT-binned datasets
DATASETS=(
    "QCD_HT300to500"
    "QCD_HT500to700"
    "QCD_HT700to1000"
    "QCD_HT1000to1500"
    "QCD_HT1500to2000"
    "QCD_HT2000toInf"
)

# Check if compiled
if [ ! -f "./create_scouting_jec_map" ]; then
    echo "Error: Executable not found. Please compile first:"
    echo "  make"
    exit 1
fi

echo ""
echo "Configuration:"
echo "  Base path: $BASE_PATH"
if [ $MAX_FILES_PER_BIN -eq -1 ]; then
    echo "  Max files per HT bin: ALL (no limit)"
else
    echo "  Max files per HT bin: $MAX_FILES_PER_BIN"
fi
echo "  Output: $OUTPUT_DIR/"
echo "  Server: $SERVER"
echo "  Datasets: ${#DATASETS[@]} QCD HT bins"
echo ""

# Create output directory
mkdir -p $OUTPUT_DIR

# Collect all files from all HT bins
echo "=============================================================="
echo "Finding ROOT files from all HT bins..."
echo "=============================================================="

ALL_FILES=""
TOTAL_FILES=0

for DATASET in "${DATASETS[@]}"; do
    DATASET_PATH="$BASE_PATH/$DATASET"
    
    # Find files from this HT bin
    if [ $MAX_FILES_PER_BIN -eq -1 ]; then
        FILES=$(xrdfs $SERVER ls $DATASET_PATH 2>/dev/null | grep ".root$" | sed "s|^|root://$SERVER/|")
    else
        FILES=$(xrdfs $SERVER ls $DATASET_PATH 2>/dev/null | grep ".root$" | head -n $MAX_FILES_PER_BIN | sed "s|^|root://$SERVER/|")
    fi
    
    if [ -z "$FILES" ]; then
        echo "  Warning: No files found in $DATASET, skipping..."
        continue
    fi
    
    N_FILES=$(echo "$FILES" | wc -l)
    TOTAL_FILES=$((TOTAL_FILES + N_FILES))
    echo "  $DATASET: $N_FILES files"
    
    # Append to master list
    ALL_FILES="$ALL_FILES $FILES"
done

echo ""
echo "Total files to process: $TOTAL_FILES"
echo ""

if [ -z "$ALL_FILES" ]; then
    echo "Error: No input files found!"
    exit 1
fi

# Process all files at once - processes BOTH Jet and FatJet simultaneously
echo "=============================================================="
echo "Processing both AK4 (Jet) and AK8 (FatJet) from all HT bins..."
echo "=============================================================="
TOTAL_START=$(date +%s)
time ./create_scouting_jec_map $ALL_FILES
TOTAL_END=$(date +%s)

# Create plots
echo ""
echo "=============================================================="
echo "Creating plots..."
echo "=============================================================="

if [ -f "$OUTPUT_DIR/scouting_jec_corrections_Jet.json" ]; then
    echo "Creating AK4 plots..."
    python3 plot_jec_corrections.py \
        $OUTPUT_DIR/scouting_jec_corrections_Jet.json \
        --output-dir $OUTPUT_DIR
fi

if [ -f "$OUTPUT_DIR/scouting_jec_corrections_FatJet.json" ]; then
    echo "Creating AK8 plots..."
    python3 plot_jec_corrections.py \
        $OUTPUT_DIR/scouting_jec_corrections_FatJet.json \
        --output-dir $OUTPUT_DIR
fi

TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))

echo ""
echo "=============================================================="
echo "Production run complete!"
echo "=============================================================="
echo "Processed $TOTAL_FILES files from ${#DATASETS[@]} HT bins"
echo "  Total time: ${TOTAL_TIME} seconds ($((TOTAL_TIME / 60)) minutes)"
echo ""
echo "Output files in $OUTPUT_DIR/:"
ls -lh $OUTPUT_DIR/scouting_jec_* 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
