import os
import re
import json
import subprocess
import argparse
import concurrent.futures
import logging
import shutil
from pathlib import Path

# Defaults
DEFAULT_LUMI = "1"
DEFAULT_SEED = "1"
# DEFAULT_QUANTILES = ["0.5"]  # For testing (only 0.5)
DEFAULT_QUANTILES = ["0.5", "0.84", "0.16", "0.975", "0.025"]
DEFAULT_OUTPUT_JSON = "limits.json"
DEFAULT_THREADS = 1

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

limit_regex = re.compile(r"Limit:\s*r\s*<\s*([0-9]+(?:\.[0-9]+)?)")

def get_prefix():
    status_dir = Path("status")
    queue = len(list(status_dir.glob("*.Queue")))
    running = len(list(status_dir.glob("*.Running")))
    complete = len(list(status_dir.glob("*.Complete")))
    return f"[Queue {queue}][Running {running}][Complete {complete}] "

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

def run_combine(datacard: Path, quantile: str, lumi: str, seed: str, logger):
    status_file = Path("status") / f"{datacard.stem}_q{quantile}.Queue"
    status_file.rename(status_file.with_suffix(".Running"))
    logger.info(get_prefix() + f"Running {datacard.name} @ {quantile}")
    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"
    cmd = [
        "combine",
        "-M", "HybridNew",
        "--LHCmode", "LHC-limits",
        datacard.name,
        "-s", str(seed),
        "--expectedFromGrid", str(quantile),
        # "--saveToys",
        "-T", "50000",
        # "--fork", "25",
        "--setParameters", f"lumiscale={lumi}",
    ]
    try:
        # save output to a log file ./logs/{datacard.stem}_q{quantile}.log
        log_dir = Path("logs")
        log_dir.mkdir(exist_ok=True)
        log_file = log_dir / f"{datacard.stem}_q{quantile}.log"
        with open(log_file, "w") as lf:
            res = subprocess.run(cmd, cwd=str(datacard.parent), env=env,
                                 stdout=lf, stderr=subprocess.STDOUT, text=True)
    except FileNotFoundError:
        raise SystemExit("combine not found in PATH. Activate your CMSSW env (cmsenv).")
    with open(log_file, "r") as lf:
        output = lf.read()
    limit = parse_limit(output)
    if limit is not None:
        logger.info(get_prefix() + f"Complete {datacard.name} @ {quantile} : r < {limit}")
    else:
        logger.info(get_prefix() + f"Failed {datacard.name} @ {quantile} return code {res.returncode}")
    status_file.with_suffix(".Running").rename(status_file.with_suffix(".Complete"))
    return limit, output, res.returncode, datacard.name, quantile

def main():
    parser = argparse.ArgumentParser(description="Run combine on datacards with multithreading.")
    parser.add_argument("--lumi", default=DEFAULT_LUMI, help="Luminosity scale (default: 1)")
    parser.add_argument("--seed", default=DEFAULT_SEED, help="Random seed (default: 1)")
    parser.add_argument("--quantiles", nargs="+", default=DEFAULT_QUANTILES, help="List of quantiles (default: ['0.5'])")
    parser.add_argument("--output", default=DEFAULT_OUTPUT_JSON, help="Output JSON file (default: limits.json)")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="Number of threads (default: 4)")
    args = parser.parse_args()

    status_dir = Path("status")
    status_dir.mkdir(exist_ok=True)

    here = Path(".").resolve()
    txt_files = sorted(p for p in here.glob("*.txt"))
    if not txt_files:
        logger.info("No .txt files found in current directory.")
        return

    all_files = {}
    for f in txt_files:
        key, mass = extract_name_and_mass(f.name)
        all_files.setdefault(key, {})
        if mass is not None:
            all_files[key]["mass"] = mass

    # Prepare tasks
    tasks = []
    for f in txt_files:
        for q in args.quantiles:
            tasks.append((f, q, args.lumi, args.seed))
            # Touch .Queue file
            status_file = status_dir / f"{f.stem}_q{q}.Queue"
            status_file.touch()

    # Run in parallel
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_task = {executor.submit(run_combine, datacard, quantile, lumi, seed, logger): (datacard, quantile) for datacard, quantile, lumi, seed in tasks}
        for future in concurrent.futures.as_completed(future_to_task):
            datacard, quantile = future_to_task[future]
            try:
                limit, log, code, fname, q = future.result()
                key, _ = extract_name_and_mass(fname)
                all_files[key][q] = limit
            except Exception as exc:
                logger.info(get_prefix() + f"Task for {datacard.name} @ {quantile} generated an exception: {exc}")

    with open(args.output, "w") as fp:
        json.dump(all_files, fp, indent=2, sort_keys=True)
    logger.info(f"Wrote {args.output}")

    shutil.rmtree(status_dir)

if __name__ == "__main__":
    main()