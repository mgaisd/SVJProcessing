#!/bin/bash

MEMORY=4GB
CORES=1
CHUNK_SIZE=500 #10000
N_WORKERS=30
#EXECUTOR=dask/etpcondor   # HTCondor at KIT ETP
PORT=3719 # port for dask scheduler, needs to be opened by admins
#N_WORKERS=6
EXECUTOR=futures     # local job
FORCE_RECREATE=0   # 1 to recreate output file if it exists, 0 else
FIRST_FILE=0
LAST_FILE=-1  # Use -1 to skim all input files

dataset_directory=/ceph/mgais/Run2ScoutingSkims_v0

module=analysis_configs.s_channel_scouting_pre_selection
selection_name=s_channel_scouting_pre_selection

#module=analysis_configs.t_channel_wnae_qcd_training_region
#selection_name=t_channel_wnae_qcd_training_region

#module=analysis_configs.t_channel_wnae_top_training_region
#selection_name=t_channel_wnae_top_training_region

#module=analysis_configs.t_channel_lost_lepton_control_region
#selection_name=t_channel_lost_lepton_control_region

year=2018

variation=nominal  # nominal jec_up jec_down jer_up jer_down

# Output directory for nominal samples - no variation of the uncertainties
#output_directory=root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/${year}/Full/PrivateSkims/${variation}
output_directory=root://cmsdcache-kit-disk.gridka.de:1094//store/user/mgaisdor/skimmed_background/


dataset_names=(
    # #
    # # Signals
    # #
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
    #QCD_HT700to1000
    #QCD_HT1000to1500
    QCD_HT1500to2000
    QCD_HT2000toInf
    QCD_HT300to500
    QCD_HT500to700

    # TTJets
    TTJets_2430000

    # # WJets
    # WJets_inclusive_260000
    # WJets_inclusive_270000
    # WJets_inclusive_280000
    # WJets_inclusive_40000
    # WJets_inclusive_50000
)

cross_sections=(
    # # see https://github.com/CMS-SVJ-scouting/SVJScouting_ntuplizer/wiki/Cross-sections
    # #
    # # Signals
    # #
    # 34.55
    # 34.55
    # 34.55

    # 23.28
    # 23.28
    # 23.28

    # 15.69
    # 15.69
    # 15.69

    # 10.57
    # 10.57
    # 10.57

    # 7.122
    # 7.122
    # 7.122

    # 4.924
    # 4.924
    # 4.924

    # 3.476
    # 3.476
    # 3.476

    # 2.498
    # 2.498
    # 2.498

    # 1.822
    # 1.822
    # 1.822

    # 0.0412
    # 0.0412
    # 0.0412

    # #
    # # Backgrounds
    # #
    # # QCD
    #6310
    #1094
    99.38
    20.20
    323400
    30140

    # TTJets
    471.7

    # # WJets
    # 52940
    # 52940
    # 52940
    # 52940
    # 52940
)


make_skims() {

    local dataset_directory=$1
    local module=$2
    local selection_name=$3
    local year=$4
    local dataset_name=$5
    local output_directory=$6
    local xsec=$7

    # Path automatically built when preparing input files lists
    local suffix_dir=${year}/${selection_name}/${dataset_name}
    local files_list_directory=${dataset_directory}/skim_input_files_list/${suffix_dir}
    local output_directory=${output_directory}/${suffix_dir}

    local output_redirector=$(echo ${output_directory} | cut -d/ -f 1-4)
    local output_dir=$(echo ${output_directory} | cut -d/ -f 4-)
    xrdfs ${output_redirector} ls ${output_dir} > /dev/null 2>&1
    if [ "$?" != "0" ]; then
        xrdfs ${output_redirector} mkdir -p ${output_dir}
    fi

    i_file=-1
    for files_list in $(ls ${files_list_directory} | sort -V); do
        ((i_file++))
        if [ ${i_file} -le ${LAST_FILE} ] || [ "${LAST_FILE}" == "-1" ]; then
            if [ ${i_file} -ge ${FIRST_FILE} ]; then

                local input_files=${files_list_directory}/${files_list}
                local output_file=${output_directory}/${files_list/.txt/.root}
                local output_file_name_tmp=$(echo ${ouput_file}_$(date +"%Y%m%d-%H%M%S") | shasum | cut -d " " -f1).root
                local output_file_tmp=/tmp/${USER}/${output_file_name_tmp}

                echo ""
                echo "Making skim file ${output_file}"

                local output_redirector=$(echo ${output_file} | cut -d/ -f 1-4)
                local output_file_path=$(echo ${output_file} | cut -d/ -f 4-)
                xrdfs ${output_redirector} ls ${output_file_path} > /dev/null 2>&1
                if [ "$?" != "0" ] || [ "${FORCE_RECREATE}" == "1" ]; then
	            if [ "${variation}" == "nominal" ]; then
		        variation_flag=''
		    else
		        variation_flag="--variation ${variation}"
		    fi
                    python skim.py -i ${input_files} -o ${output_file_tmp} -p ${module} -pd ${dataset_name} -y ${year} -nano_scout -xsec ${xsec} -e ${EXECUTOR} -port ${PORT} -n ${N_WORKERS} -c ${CHUNK_SIZE} --memory ${MEMORY} --cores ${CORES} -pn_tagger ${variation_flag[@]}
                    xrdcp -f ${output_file_tmp} ${output_file}
                    echo ${output_file} has been saved.
                    rm ${output_file_tmp}
                else
                    echo ${output_file} already exists and FORCE_RECREATE is 0. Skipping.
                fi
            fi
        fi
    done
}


#for dataset_name in ${dataset_names[@]}; do
#    make_skims ${dataset_directory} ${module} ${selection_name} ${year} ${dataset_name} ${output_directory}
#done

n_datasets=${#dataset_names[@]}

for ((i=0; i<$n_datasets; i++)); do

    dataset_name=${dataset_names[i]}
    cross_section=${cross_sections[i]}
    make_skims ${dataset_directory} ${module} ${selection_name} ${year} ${dataset_name} ${output_directory} ${cross_section}
done
