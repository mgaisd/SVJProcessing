import ROOT

filelist = "input_files.txt"  # list of xrootd paths (one per line)
broken_count = 0  # Counter for broken files

with open(filelist) as f:
    for line in f:
        path = line.strip()
        print(f"Checking: {path}")
        try:
            f = ROOT.TFile.Open(path)
            if not f or f.IsZombie():
                print(f"[BROKEN] {path}")
                broken_count += 1  # Increment broken count
            else:
                print(f"[OK] {path}")
            f.Close()
        except Exception as e:
            print(f"[ERROR] {path} - {e}")
            broken_count += 1  # Increment broken count for errors

print(f"\nTotal broken files: {broken_count}/{len(open(filelist).readlines())}")
