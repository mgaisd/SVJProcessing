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


for year in years:
    datasets_info[year].update({
        "scouting_data": [
            {
                "redirector": "root://cmsdcache-kit-disk.gridka.de:1094/",
                "path": f"/store/user/mgaisdor/SVJScouting_ntuples/Data/{year}/",
                "regex": f"",
            },
        ]
    })
