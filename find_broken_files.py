from multiprocessing import Pool
from subprocess import run, DEVNULL
from tqdm import tqdm
import subprocess

redirector = "root://cmsdcache-kit-disk.gridka.de:1094"

# Get list of root files
files = subprocess.check_output(
    ["xrdfs", redirector, "ls", "-R", "/store/user/mgaisdor/SVJScouting_ntuples/MC/2018_v2/"],
    text=True,
).splitlines()
files = [f for f in files if f.endswith(".root")]

def check_file(file):
    fullpath = f"{redirector}{file}"
    result = run(["rootls", fullpath], stdout=DEVNULL, stderr=DEVNULL)
    return file if result.returncode != 0 else None

with Pool(16) as p:
    results = list(tqdm(p.imap(check_file, files), total=len(files), desc="Checking"))

# Print broken files
for f in results:
    if f:
        print(f"Broken: {f}")
