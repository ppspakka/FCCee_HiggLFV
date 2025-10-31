import os
import re
import json
import subprocess
from pathlib import Path

# Defaults
LUMI = os.environ.get("LUMI", "1")
SEED = os.environ.get("SEED", "1")
QUANTILES = ["0.5", "0.84", "0.16", "0.975", "0.025"]
OUTPUT_JSON = "limits.json"

limit_regex = re.compile(r"Limit:\s*r\s*<\s*([0-9]+(?:\.[0-9]+)?)")

def parse_limit(log: str):
    m = limit_regex.search(log)
    return float(m.group(1)) if m else None

def extract_name_and_mass(filename: str):
    # Expect "datacard_HMuTauE_LFV_145.txt" -> ("HMuTauE_LFV_145", 145)
    m = re.match(r"datacard_(.+)_(\d+)\.txt$", filename)
    if m:
        name = f"{m.group(1)}_{m.group(2)}"
        return name, int(m.group(2))
    # Fallback if pattern doesn’t match
    stem = Path(filename).stem
    return stem, None

def run_combine(datacard: Path, quantile: str):
    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"
    cmd = [
        "combine",
        "-M", "HybridNew",
        "--LHCmode", "LHC-limits",
        datacard.name,
        "-s", str(SEED),
        "--expectedFromGrid", str(quantile),
        "--saveToys",
        "-T", "5000",
        "--setParameters", f"lumiscale={LUMI}",
    ]
    try:
        res = subprocess.run(cmd, cwd=str(datacard.parent), env=env,
                             capture_output=True, text=True)
    except FileNotFoundError:
        raise SystemExit("combine not found in PATH. Activate your CMSSW env (cmsenv).")
    output = (res.stdout or "") + "\n" + (res.stderr or "")
    return parse_limit(output), output, res.returncode

def main():
    here = Path(".").resolve()
    txt_files = sorted(p for p in here.glob("*.txt"))
    if not txt_files:
        print("No .txt files found in current directory.")
        return

    all_files = {}
    for f in txt_files:
        key, mass = extract_name_and_mass(f.name)
        all_files.setdefault(key, {})
        if mass is not None:
            all_files[key]["mass"] = mass

        for q in QUANTILES:
            print(f"Running combine: file={f.name}, quantile={q}, LUMI={LUMI}")
            limit, log, code = run_combine(f, q)
            if limit is None:
                print(f"Warning: no limit parsed for {f.name} @ {q}. Return code={code}")
            all_files[key][q] = limit

    with open(OUTPUT_JSON, "w") as fp:
        json.dump(all_files, fp, indent=2, sort_keys=True)
    print(f"Wrote {OUTPUT_JSON}")

if __name__ == "__main__":
    main()