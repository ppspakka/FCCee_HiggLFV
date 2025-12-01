#!/usr/bin/env python3
"""
Parallel runner converted from runall.sh

- Builds a queue of jobs (signals & backgrounds) which can be enabled/disabled via CONFIG.
- Runs jobs in parallel (default workers = 4).
- Saves per-job log file next to the output, and also prints logs to the terminal (prefixed).
- After successful completion, runs makecard.py for each processed channel.
"""
from pathlib import Path
import subprocess
import sys
import concurrent.futures
import argparse
from dataclasses import dataclass
from typing import List, Dict, Any

# -------------------------
# User-configurable section
# -------------------------
PARALLEL = 4  # default parallel workers

# Pipeline files
MUTAUE_PIPELINE = "./pipeline_mutaue.json"
ETAUMU_PIPELINE = "./pipeline_etaumu.json"
MUTAUE_OFFSHELL_PIPELINE = "./pipeline_mutaue_Zoffshell.json"
ETAUMU_OFFSHELL_PIPELINE = "./pipeline_etaumu_Zoffshell.json"

# Save path
PARENT_DIR = Path("./normal_cuts")
MUTAUE_DIR = PARENT_DIR / "mutaue_hist"
ETAUMU_DIR = PARENT_DIR / "etaumu_hist"
MUTAUE_OFFSHELL_DIR = PARENT_DIR / "mutaue_Zoffshell_hist"
ETAUMU_OFFSHELL_DIR = PARENT_DIR / "etaumu_Zoffshell_hist"

# Backgrounds
BACKGROUNDS = [
    "/work/project/physics/psriling/FCC/FCCee/ZWW/ZWW.root",
    # "/work/project/physics/psriling/FCC/FCCee/HZFourLepton/HZFourLep.root", # Unused
    "/work/project/physics/psriling/FCC/FCCee/zh_ll_ww/zh_ll_ww.root",
    "/work/project/physics/psriling/FCC/FCCee/zh_ll_tautau/zh_ll_tautau.root",
    "/work/project/physics/psriling/FCC/FCCee/zz_ll_tautau/zz_ll_tautau.root",
]

# Toggle which channels / types to include (set True/False)
RUN_SIGNALS = True
RUN_BACKGROUNDS = True
CONFIG = {
    "mutaue_signal": RUN_SIGNALS,
    "etaumu_signal": RUN_SIGNALS,
    
    "mutaue_offshell_signal": RUN_SIGNALS,
    "etaumu_offshell_signal": RUN_SIGNALS,
    
    "mutaue_background": RUN_BACKGROUNDS,
    "etaumu_background": RUN_BACKGROUNDS,
    
    "mutaue_offshell_background": RUN_BACKGROUNDS,
    "etaumu_offshell_background": RUN_BACKGROUNDS,
}
# -------------------------
# End user-configurable
# -------------------------

@dataclass
class Job:
    input_path: str
    output_dir: Path
    pipeline: str
    sample_type: str  # "signal" or "background"
    channel: str      # e.g. "mutaue", "etaumu", ...

    def out_root(self) -> Path:
        filename = Path(self.input_path).name
        return self.output_dir / f"{self.sample_type}_{filename}"

    def log_path(self) -> Path:
        return self.output_dir / f"{self.sample_type}_{Path(self.input_path).name}.txt"

    def command(self) -> List[str]:
        # root expects the analyze_pipeline.cpp(...) as a single argument
        cpp_arg = f'analyze_pipeline.cpp("{self.input_path}","{self.out_root()}","{self.pipeline}")'
        return ["root", "-l", "-b", "-q", cpp_arg]


def ensure_dirs():
    for d in (MUTAUE_DIR, ETAUMU_DIR, MUTAUE_OFFSHELL_DIR, ETAUMU_OFFSHELL_DIR):
        d.mkdir(parents=True, exist_ok=True)


def build_jobs() -> List[Job]:
    jobs: List[Job] = []

    # MuTauE signals
    if CONFIG["mutaue_signal"]:
        for mass in range(110, 161, 5):
            inp = f"/work/project/physics/psriling/FCC/FCCee/HMuTauE_LFV/Hmass{mass}/HMuTauE_LFV_{mass}.root"
            jobs.append(Job(input_path=inp, output_dir=MUTAUE_DIR, pipeline=MUTAUE_PIPELINE,
                            sample_type="signal", channel="mutaue"))
    if CONFIG["mutaue_background"]:
        for bg in BACKGROUNDS:
            jobs.append(Job(input_path=bg, output_dir=MUTAUE_DIR, pipeline=MUTAUE_PIPELINE,
                            sample_type="background", channel="mutaue"))

    # Etaumu signals
    if CONFIG["etaumu_signal"]:
        for mass in range(110, 161, 5):
            inp = f"/work/project/physics/psriling/FCC/FCCee/HETauMu_LFV/Hmass{mass}/HETauMu_LFV_{mass}.root"
            jobs.append(Job(input_path=inp, output_dir=ETAUMU_DIR, pipeline=ETAUMU_PIPELINE,
                            sample_type="signal", channel="etaumu"))
    if CONFIG["etaumu_background"]:
        for bg in BACKGROUNDS:
            jobs.append(Job(input_path=bg, output_dir=ETAUMU_DIR, pipeline=ETAUMU_PIPELINE,
                            sample_type="background", channel="etaumu"))

    # MuTauE Offshell
    if CONFIG["mutaue_offshell_signal"]:
        for mass in range(150, 161, 5):
            inp = f"/work/project/physics/psriling/FCC/FCCee/HMuTauE_LFV/Hmass{mass}/HMuTauE_LFV_{mass}.root"
            jobs.append(Job(input_path=inp, output_dir=MUTAUE_OFFSHELL_DIR, pipeline=MUTAUE_OFFSHELL_PIPELINE,
                            sample_type="signal", channel="mutaue_offshell"))
    if CONFIG["mutaue_offshell_background"]:
        for bg in BACKGROUNDS:
            jobs.append(Job(input_path=bg, output_dir=MUTAUE_OFFSHELL_DIR, pipeline=MUTAUE_OFFSHELL_PIPELINE,
                            sample_type="background", channel="mutaue_offshell"))

    # Etaumu Offshell
    if CONFIG["etaumu_offshell_signal"]:
        for mass in range(150, 161, 5):
            inp = f"/work/project/physics/psriling/FCC/FCCee/HETauMu_LFV/Hmass{mass}/HETauMu_LFV_{mass}.root"
            jobs.append(Job(input_path=inp, output_dir=ETAUMU_OFFSHELL_DIR, pipeline=ETAUMU_OFFSHELL_PIPELINE,
                            sample_type="signal", channel="etaumu_offshell"))
    if CONFIG["etaumu_offshell_background"]:
        for bg in BACKGROUNDS:
            jobs.append(Job(input_path=bg, output_dir=ETAUMU_OFFSHELL_DIR, pipeline=ETAUMU_OFFSHELL_PIPELINE,
                            sample_type="background", channel="etaumu_offshell"))

    return jobs


def run_job(job: Job, verbose_prefix: bool = True) -> int:
    job.output_dir.mkdir(parents=True, exist_ok=True)
    log_path = job.log_path()
    cmd = job.command()

    # Start subprocess; capture stdout+stderr merged
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)

    # Open log file
    with open(log_path, "w", encoding="utf-8", errors="replace") as logf:
        prefix = f"[{job.channel}/{job.sample_type}/{Path(job.input_path).name}] "
        # Stream output: write both to file and to terminal with prefix
        if proc.stdout is None:
            return proc.wait()
        for line in proc.stdout:
            # write to terminal
            try:
                sys.stdout.write(prefix + line)
                sys.stdout.flush()
            except BrokenPipeError:
                pass
            # write to log
            logf.write(line)
        return proc.wait()

def run_makecard_commands(dry_run: bool = False):
    """Builds and runs the makecard.py commands for channels that were processed."""
    print("\n" + "="*30)
    print("Starting makecard generation")
    print("="*30)

    makecard_jobs: Dict[str, Any] = {
        "mutaue": {
            "in_dir": MUTAUE_DIR,
            "out_dir": Path(f"{PARENT_DIR}/datacards_{MUTAUE_DIR.name.replace('_hist', '')}"),
            "enabled": CONFIG["mutaue_signal"] or CONFIG["mutaue_background"]
        },
        "etaumu": {
            "in_dir": ETAUMU_DIR,
            "out_dir": Path(f"{PARENT_DIR}/datacards_{ETAUMU_DIR.name.replace('_hist', '')}"),
            "enabled": CONFIG["etaumu_signal"] or CONFIG["etaumu_background"]
        },
        "mutaue_offshell": {
            "in_dir": MUTAUE_OFFSHELL_DIR,
            "out_dir": Path(f"{PARENT_DIR}/datacards_{MUTAUE_OFFSHELL_DIR.name.replace('_hist', '')}"),
            "enabled": CONFIG["mutaue_offshell_signal"] or CONFIG["mutaue_offshell_background"]
        },
        "etaumu_offshell": {
            "in_dir": ETAUMU_OFFSHELL_DIR,
            "out_dir": Path(f"{PARENT_DIR}/datacards_{ETAUMU_OFFSHELL_DIR.name.replace('_hist', '')}"),
            "enabled": CONFIG["etaumu_offshell_signal"] or CONFIG["etaumu_offshell_background"]
        },
    }

    lumi_pb = 1_000_000
    failures = 0

    for name, params in makecard_jobs.items():
        if not params["enabled"]:
            continue

        out_dir = params["out_dir"]
        out_dir.mkdir(exist_ok=True)
        
        cmd = [
            "python3", "makecard.py",
            "--in-dir", str(params["in_dir"]),
            "--lumi-pb", str(lumi_pb),
            "--out-root", str(out_dir / "merged.root"),
            "--out-card", str(out_dir)
        ]

        print(f"\nRunning makecard for '{name}':")
        print("CMD:", " ".join(cmd))

        if dry_run:
            print("DRY RUN: command not executed.")
            continue

        try:
            # Using subprocess.run for simpler commands
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(result.stdout)
            if result.stderr:
                print("--- STDERR ---")
                print(result.stderr)
            print(f"Makecard for '{name}' completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"ERROR: makecard for '{name}' failed with exit code {e.returncode}.")
            print("--- STDOUT ---")
            print(e.stdout)
            print("--- STDERR ---")
            print(e.stderr)
            failures += 1
        except FileNotFoundError:
            print("ERROR: 'makecard.py' not found. Make sure it is in the current directory.")
            failures += 1
            # No point in trying other jobs if script is missing
            break
    
    if failures > 0:
        print(f"\n{failures} makecard job(s) failed.")
        sys.exit(3)
    else:
        print("\nAll makecard jobs completed successfully.")


def main():
    parser = argparse.ArgumentParser(description="Parallel runner for analyze_pipeline jobs")
    parser.add_argument("--output-dir", "-o", type=str, default=str(PARENT_DIR), help="parent output directory")
    parser.add_argument("--parallel", "-p", type=int, default=PARALLEL, help="number of parallel jobs")
    parser.add_argument("--list", action="store_true", help="only list jobs (don't execute)")
    parser.add_argument("--dry-run", action="store_true", help="show commands without running")
    parser.add_argument("--skip-makecard", action="store_true", help="skip the final makecard step")
    args = parser.parse_args()

    ensure_dirs()
    jobs = build_jobs()

    if not jobs:
        print("No jobs to run (check CONFIG).")
        return

    print(f"Built {len(jobs)} jobs. Running with {args.parallel} parallel workers.")
    if args.list:
        for j in jobs:
            print("CMD:", " ".join(j.command()), "-> log:", j.log_path())
        return

    if args.dry_run:
        for j in jobs:
            print("DRY:", " ".join(j.command()), "->", j.out_root(), "log:", j.log_path())
        if not args.skip_makecard:
            run_makecard_commands(dry_run=True)
        return

    # Run in ThreadPoolExecutor: subprocesses are external so threads are fine for IO.
    failures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.parallel) as exe:
        future_to_job = {exe.submit(run_job, job): job for job in jobs}
        for fut in concurrent.futures.as_completed(future_to_job):
            job = future_to_job[fut]
            try:
                rc = fut.result()
                if rc != 0:
                    print(f"Job FAILED (rc={rc}): {job.input_path} -> see {job.log_path()}")
                    failures.append((job, rc))
                else:
                    print(f"Job DONE: {job.input_path} -> {job.out_root()}")
            except Exception as e:
                print(f"Job EXCEPTION for {job.input_path}: {e}")
                failures.append((job, -1))

    if failures:
        print(f"\n{len(failures)} jobs failed. Check logs.")
        print("Skipping makecard step due to failures.")
        sys.exit(2)
    else:
        print("\nAll processing completed successfully.")
        if not args.skip_makecard:
            run_makecard_commands()
        else:
            print("\nSkipping makecard step as requested.")


if __name__ == "__main__":
    main()