#!/bin/bash

N_WORKERS=6

dataset_directory=/ceph/mgais/Run2ScoutingSkims_v1

selection_name=s_channel_scouting_pre_selection_validation_region

year=2018

# Output directory for nominal samples - no variation of the uncertainties
output_directory=root://cmsdcache-kit-disk.gridka.de:1094//store/user/mgaisdor/SVJScouting_skims


dataset_names=(
    scouting_data
)

check_number_of_events() {

    local dataset_directory=$1
    local selection_name=$2
    local year=$3
    local dataset_name=$4
    local output_directory=$5

    # Path automatically built when preparing input files lists
    local files_list=${dataset_directory}/files_list/${year}/${dataset_name}.csv
    local output_directory=${output_directory}/${year}/${selection_name}/${dataset_name}

    python check_number_of_events.py -i ${files_list} -o ${output_directory} -n ${N_WORKERS} -nano
}


for dataset_name in ${dataset_names[@]}; do

    check_number_of_events ${dataset_directory} ${selection_name} ${year} ${dataset_name} ${output_directory}

done
