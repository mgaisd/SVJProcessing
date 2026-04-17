from multiprocessing import Pool
from subprocess import run, DEVNULL
from tqdm import tqdm
import subprocess

redirector = "root://cmsdcache-kit-disk.gridka.de:1094"

# Get list of root files
files = subprocess.check_output(
    ["xrdfs", redirector, "ls", "-R", "/store/user/mgaisdor/SVJScouting_ntuples/MC/2017_v2/QCD_HT700to1000/"],
    text=True,
).splitlines()
files = [f for f in files if f.endswith(".root")]

def check_file(file):
    fullpath = f"{redirector}{file}"
    result = run(["rootls", fullpath], stdout=DEVNULL, stderr=DEVNULL)
    # Return the file if it's broken (non-zero exit code), None if it's good
    if result.returncode != 0:
        return file
    return None

with Pool(16) as p:
    results = list(tqdm(p.imap(check_file, files), total=len(files), desc="Checking"))

# Print broken files and summary
broken_files = [f for f in results if f is not None]
print(f"\nFound {len(broken_files)} broken files out of {len(files)} total files:")
for f in broken_files:
    print(f"Broken: {f}")

if not broken_files:
    print("No broken files found!")
