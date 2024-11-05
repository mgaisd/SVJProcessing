###################################  README  ###################################
#
# Is called "dataset" a set of files corresponding to the same physics process.
# The object `datasets_info` describes the location of the different datasets.
# Its structure is the following:
#    * Keys are year
#    * Values are the "datasets_per_year_info"
# The structure of the `datasets_per_year_info` is the following:
#    * Keys are dataset names
#    * Values are the "dataset_info" defining which files belong to the dataset
# 
# The "dataset_info" has the following structure. It is a list of dict, which
# has 2 keys:
#    * "redirector": The XRootD redirector to the remote storage element
#    * "path": The path to the directory at which the files are located
#    * "regex": The regex to apply to select some files from that directory. 
#               The regex must be "" if no regex is applied.
#
################################################################################


years = ["2018"]

datasets_info = {
    year: {} for year in years
}

# SVJ scouting signal models
signal_models = [
    "s-channel_mMed-700_mDark-20_rinv-0.3",
    "s-channel_mMed-700_mDark-20_rinv-0.5",
    "s-channel_mMed-700_mDark-20_rinv-0.7",

    "s-channel_mMed-800_mDark-20_rinv-0.3",
    "s-channel_mMed-800_mDark-20_rinv-0.5",
    "s-channel_mMed-800_mDark-20_rinv-0.7",

    "s-channel_mMed-900_mDark-20_rinv-0.3",
    "s-channel_mMed-900_mDark-20_rinv-0.5",
    "s-channel_mMed-900_mDark-20_rinv-0.7",

    "s-channel_mMed-1000_mDark-20_rinv-0.3",
    "s-channel_mMed-1000_mDark-20_rinv-0.5",
    "s-channel_mMed-1000_mDark-20_rinv-0.7",

    "s-channel_mMed-1100_mDark-20_rinv-0.3",
    "s-channel_mMed-1100_mDark-20_rinv-0.5",
    "s-channel_mMed-1100_mDark-20_rinv-0.7",

    "s-channel_mMed-1200_mDark-20_rinv-0.3",
    "s-channel_mMed-1200_mDark-20_rinv-0.5",
    "s-channel_mMed-1200_mDark-20_rinv-0.7",

    "s-channel_mMed-1300_mDark-20_rinv-0.3",
    "s-channel_mMed-1300_mDark-20_rinv-0.5",
    "s-channel_mMed-1300_mDark-20_rinv-0.7",

    "s-channel_mMed-1400_mDark-20_rinv-0.3",
    "s-channel_mMed-1400_mDark-20_rinv-0.5",
    "s-channel_mMed-1400_mDark-20_rinv-0.7",

    "s-channel_mMed-1500_mDark-20_rinv-0.3",
    "s-channel_mMed-1500_mDark-20_rinv-0.5",
    "s-channel_mMed-1500_mDark-20_rinv-0.7",

    "s-channel_mMed-3000_mDark-20_rinv-0.3",
    "s-channel_mMed-3000_mDark-20_rinv-0.5",
    "s-channel_mMed-3000_mDark-20_rinv-0.7",
]

# signal_models2 = [
#     "SVJ_hadronic_std2_s-channel_mMed-700_mDark-20_rinv-0.3",
#     "SVJ_hadronic_std2_s-channel_mMed-700_mDark-20_rinv-0.5",
#     "SVJ_hadronic_std2_s-channel_mMed-700_mDark-20_rinv-0.7",

#     "SVJ_hadronic_std2_s-channel_mMed-800_mDark-20_rinv-0.3",
#     "SVJ_hadronic_std2_s-channel_mMed-800_mDark-20_rinv-0.5",
#     "SVJ_hadronic_std2_s-channel_mMed-800_mDark-20_rinv-0.7",

#     "SVJ_hadronic_std2_s-channel_mMed-900_mDark-20_rinv-0.3",
#     "SVJ_hadronic_std2_s-channel_mMed-900_mDark-20_rinv-0.5",
#     "SVJ_hadronic_std2_s-channel_mMed-900_mDark-20_rinv-0.7",

#     "SVJ_hadronic_std2_s-channel_mMed-1000_mDark-20_rinv-0.3",
#     "SVJ_hadronic_std2_s-channel_mMed-1000_mDark-20_rinv-0.5",
#     "SVJ_hadronic_std2_s-channel_mMed-1000_mDark-20_rinv-0.7",

#     "SVJ_hadronic_std2_s-channel_mMed-1500_mDark-20_rinv-0.3",
#     "SVJ_hadronic_std2_s-channel_mMed-1500_mDark-20_rinv-0.5",
#     "SVJ_hadronic_std2_s-channel_mMed-1500_mDark-20_rinv-0.7",
# ]

# signal_models3 = [
#     "s-channel_mMed-700GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-700GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-700GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV-pythia8_n-2000",

#     "s-channel_mMed-800GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-800GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-800GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV-pythia8_n-2000",

#     "s-channel_mMed-900GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-900GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-900GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV-pythia8_n-2000",

#     "s-channel_mMed-1100GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-1100GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-1100GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV-pythia8_n-2000",

#     "s-channel_mMed-1200GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-1200GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-1200GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV-pythia8_n-2000",

#     "s-channel_mMed-1300GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-1300GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-1300GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV-pythia8_n-2000",

#     "s-channel_mMed-1400GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-1400GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-1400GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV-pythia8_n-2000",

#     "s-channel_mMed-1500GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-1500GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV-pythia8_n-2000",
#     "s-channel_mMed-1500GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV-pythia8_n-2000",
# ]

qcd_bins = [
    "QCD_HT100to200",
    "QCD_HT200to300",
    "QCD_HT300to500",
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf",
]

ttjets_bins = [
    "TTJets_2430000",
]

wjets_bins = [
    "WJets_inclusive_260000",
    "WJets_inclusive_270000",
    "WJets_inclusive_280000",
    "WJets_inclusive_40000",
    "WJets_inclusive_50000",
]

for year in years:
    for signal_model in signal_models:
        if "3000" in signal_model:
            datasets_info[year].update({
                signal_model: [
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096//",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/signals/SVJ_hadronic_std/{signal_model}/",
                        "regex": f"",
                    },
                ]
            })
        elif "1000" in signal_model:
            datasets_info[year].update({
                signal_model: [
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096//",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/signals/SVJ_hadronic_std/{signal_model}/",
                        "regex": f"",
                    },
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096//",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/signals/SVJ_hadronic_std2/{signal_model}/",
                        "regex": f"",
                    
                    },
                ]
            })
        elif ("1100" in signal_model) or ("1200" in signal_model) or ("1300" in signal_model) or ("1400" in signal_model):
            datasets_info[year].update({
                signal_model: [
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096//",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/signals/SVJ_hadronic_std3/{signal_model}/",
                        "regex": f"",
                    },
                ]
            })
        else:
            datasets_info[year].update({
                signal_model: [
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096//",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/signals/SVJ_hadronic_std/{signal_model}/",
                        "regex": f"",
                    },
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096//",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/signals/SVJ_hadronic_std2/{signal_model}/",
                        "regex": f"",
                    
                    },
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096//",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/signals/SVJ_hadronic_std3/{signal_model}_alpha-peak_13TeV-pythia8_n-2000/",
                        "regex": f"",
                    },
                ]
            })         

# for year in years:
#     datasets_info[year].update({
#         signal_model: [
#             {
#                 "redirector": "root://storage01.lcg.cscs.ch:1096//",
#                 "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/signals/SVJ_hadronic_std3/{signal_model}/",
#                 "regex": f"",
#             },
#         ]
#         for signal_model in signal_models3
#     })


for year in years:
    datasets_info[year].update({
        bin: [
            {
                "redirector": "root://storage01.lcg.cscs.ch:1096//",
                "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/QCD_HT_binned_2018_v0/{bin}/",
                "regex": f"",
                
            },
        ]
        for bin in qcd_bins
    })

for year in years:
    datasets_info[year].update({
        bin: [
            {
                "redirector": "root://storage01.lcg.cscs.ch:1096//",
                "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/TTJets_inclusive_2018_v0/{bin}/",
                "regex": f"",
                
            },
        ]
        for bin in ttjets_bins
    })

for year in years:
    datasets_info[year].update({
        bin: [
            {
                "redirector": "root://storage01.lcg.cscs.ch:1096//",
                "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/PFNano/WJets_inclusive_2018_v0/{bin}/",
                "regex": f"",
                
            },
        ]
        for bin in wjets_bins
    })