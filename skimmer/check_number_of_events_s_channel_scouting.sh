#!/bin/bash

N_WORKERS=6

dataset_directory=/ceph/mgais/Run2ScoutingSkims_v0

selection_name=s_channel_scouting_pre_selection
#selection_name=t_channel_wnae_qcd_training_region
#selection_name=t_channel_wnae_top_training_region
#selection_name=t_channel_lost_lepton_control_region

year=2018

# Output directory for nominal samples - no variation of the uncertainties
output_directory=root://cmsdcache-kit-disk.gridka.de:1094//store/user/mgaisdor/skimmed_background/


dataset_names=(
    #
    # Signals
    #
    # s-channel_mMed-700_mDark-20_rinv-0.3
    # s-channel_mMed-700_mDark-20_rinv-0.5
    # s-channel_mMed-700_mDark-20_rinv-0.7

    # s-channel_mMed-800_mDark-20_rinv-0.3
    # s-channel_mMed-800_mDark-20_rinv-0.5
    # s-channel_mMed-800_mDark-20_rinv-0.7

    # s-channel_mMed-900_mDark-20_rinv-0.3
    # s-channel_mMed-900_mDark-20_rinv-0.5
    # s-channel_mMed-900_mDark-20_rinv-0.7

    # s-channel_mMed-1000_mDark-20_rinv-0.3
    # s-channel_mMed-1000_mDark-20_rinv-0.5
    # s-channel_mMed-1000_mDark-20_rinv-0.7

    # s-channel_mMed-1100_mDark-20_rinv-0.3
    # s-channel_mMed-1100_mDark-20_rinv-0.5
    # s-channel_mMed-1100_mDark-20_rinv-0.7

    # s-channel_mMed-1200_mDark-20_rinv-0.3
    # s-channel_mMed-1200_mDark-20_rinv-0.5
    # s-channel_mMed-1200_mDark-20_rinv-0.7

    # s-channel_mMed-1300_mDark-20_rinv-0.3
    # s-channel_mMed-1300_mDark-20_rinv-0.5
    # s-channel_mMed-1300_mDark-20_rinv-0.7

    # s-channel_mMed-1400_mDark-20_rinv-0.3
    # s-channel_mMed-1400_mDark-20_rinv-0.5
    # s-channel_mMed-1400_mDark-20_rinv-0.7

    # s-channel_mMed-1500_mDark-20_rinv-0.3
    # s-channel_mMed-1500_mDark-20_rinv-0.5
    # s-channel_mMed-1500_mDark-20_rinv-0.7

    # s-channel_mMed-3000_mDark-20_rinv-0.3
    # s-channel_mMed-3000_mDark-20_rinv-0.5
    # s-channel_mMed-3000_mDark-20_rinv-0.7

    #
    # Backgrounds
    #
    # QCD
    QCD_HT300to500
    QCD_HT500to700
    QCD_HT700to1000
    QCD_HT1000to1500
    QCD_HT1500to2000
    QCD_HT2000toInf

    TTJets
    TTJets_2430000

    WJets
    WJets_inclusive_260000
    WJets_inclusive_270000
    WJets_inclusive_280000
    WJets_inclusive_40000
    WJets_inclusive_50000
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

