from multiprocessing.pool import ThreadPool
from tqdm import tqdm
import subprocess
import uproot

redirector = "root://cmsdcache-kit-disk.gridka.de:1094/"

# Get list of root files
files = subprocess.check_output(
    ["xrdfs", redirector, "ls", "-R", "/store/user/mgaisdor/SVJScouting_ntuples/MC/2017_v2/"],
    text=True,
).splitlines()
files = [f for f in files if f.endswith(".root")]

def check_file(file):
    fullpath = f"{redirector}{file}"
    try:
        with uproot.open(fullpath) as f:
            _ = f.keys()
        return None, None
    except Exception as e:
        return file, str(e)

with ThreadPool(64) as p:
    results = list(tqdm(p.imap(check_file, files), total=len(files), desc="Checking"))

broken_files = [(f, e) for f, e in results if f is not None]
print(f"\nFound {len(broken_files)} broken files out of {len(files)} total files:")
for f, e in broken_files:
    print(f"Broken: {f}\n  Error: {e}\n")

if not broken_files:
    print("No broken files found!")
