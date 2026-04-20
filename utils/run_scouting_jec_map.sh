#!/bin/bash
#
# Helper script to run the scouting JEC map creation
# This script sets up the environment and runs the correction map creator
#

# Source the main setup if it exists
if [ -f setup.sh ]; then
    source setup.sh
fi

# Default settings
YEAR="2017"
JET_TYPE="both"
DR_THRESHOLD="0.15"
N_FILES=""
MAX_JETS="--max-jets 2"
PLOT_FLAG=""

# Parse arguments
MAX_JETS_CHANGED=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --test)
            N_FILES="--n-files 1"
            echo "Test mode: Processing only 1 file per dataset"
            shift
            ;;
        --n-files)
            N_FILES="--n-files $2"
            echo "Processing only $2 file(s) per dataset"
            shift 2
            ;;
        --max-jets)
            MAX_JETS="--max-jets $2"
            MAX_JETS_CHANGED=true
            echo "Using only $2 leading FatJet(s) per event (Jet uses all jets)"
            shift 2
            ;;
        --no-plot)
            PLOT_FLAG="--no-plot"
            shift
            ;;
        --year)
            YEAR="$2"
            shift 2
            ;;
        --jet-type)
            JET_TYPE="$2"
            shift 2
            ;;
        --dr-threshold)
            DR_THRESHOLD="$2"
            shift 2
            ;;
        --help)
            python utils/create_scouting_jec_map.py --help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--test] [--n-files N] [--max-jets N] [--no-plot] [--year YEAR] [--jet-type TYPE] [--dr-threshold DR] [--help]"
            exit 1
            ;;
    esac
done

# Print default max-jets if not changed
if [ "$MAX_JETS_CHANGED" = false ]; then
    echo "Using only 2 leading FatJet(s) per event (Jet uses all jets) [default]"
fi

# Run the script
echo "Running scouting JEC map creation..."
echo "  Year: $YEAR"
echo "  Jet type: $JET_TYPE"
echo "  deltaR threshold: $DR_THRESHOLD"

python utils/create_scouting_jec_map.py \
    --year $YEAR \
    --jet-type $JET_TYPE \
    --dr-threshold $DR_THRESHOLD \
    $N_FILES \
    $MAX_JETS \
    $PLOT_FLAG

echo ""
echo "Done! Check the following directories:"
echo "  Corrections: data/scouting_jec_residuals_${YEAR}_*.coffea"
echo "  Plots: jec_plots/scouting_correction_map_*.png"
