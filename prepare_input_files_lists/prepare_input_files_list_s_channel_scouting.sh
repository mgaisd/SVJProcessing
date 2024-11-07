#!/bin/bash

dataset_directory=/ceph/mgais/Run2ScoutingSkims_v0
dataset_config=dataset_configs.s_channel_scouting_new_signal_paths

module=analysis_configs.s_channel_scouting_pre_selection
selection_name=s_channel_scouting_pre_selection

#module=analysis_configs.t_channel_wnae_qcd_training_region
#selection_name=t_channel_wnae_qcd_training_region

#module=analysis_configs.t_channel_wnae_top_training_region
#selection_name=t_channel_wnae_top_training_region

#module=analysis_configs.t_channel_lost_lepton_control_region
#selection_name=t_channel_lost_lepton_control_region

year=2018

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
    #
    QCD_HT200to300
    QCD_HT300to500
    QCD_HT500to700
    QCD_HT700to1000
    QCD_HT1000to1500
    QCD_HT1500to2000
    QCD_HT2000toInf
    QCD_HT100to200

    #
    # TTJets
    #
    TTJets_2430000

    #
    # WJets
    #
    WJets_inclusive_260000
    WJets_inclusive_270000
    WJets_inclusive_280000
    WJets_inclusive_40000
    WJets_inclusive_50000
)

prepare_input_files_list() {

    local dataset_config=$1
    local dataset_directory=$2
    local module=$3
    local selection_name=$4
    local year=$5
    local dataset_name=$6

    echo ""
    echo "Preparing input files for dataset ${dataset_name} year ${year} and selection ${selection_name}"

    python list_dataset_files.py -d ${dataset_name} -y ${year} -c ${dataset_config} -o ${dataset_directory} -nano_scout
    python compute_unweighted_selection_efficiency.py -d ${dataset_name} -y ${year} -p ${module} -s ${selection_name} -i ${dataset_directory} -o ${dataset_directory} -n 6 -e futures -c 10000 -nano_scout
    python prepare_input_files_list.py -d ${dataset_name} -y ${year} -s ${selection_name} -i ${dataset_directory} -o ${dataset_directory} -m 50000
}


for dataset_name in ${dataset_names[@]}; do

    prepare_input_files_list ${dataset_config} ${dataset_directory} ${module} ${selection_name} ${year} ${dataset_name}

done

