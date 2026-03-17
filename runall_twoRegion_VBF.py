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
import os
import concurrent.futures
import argparse
from dataclasses import dataclass
import shutil
from typing import List, Dict, Any

# REVISED SCRIPT FOR VBF ANALYSIS
# INPUT FILE FORMAT: /work/project/physics/vwachira/fcc-ee-higgs-lfv/delphes_outputs/XX
# XX = e.g. mh110_etau, ... (or mh110_mutau)

# -------------------------
# User-configurable section
# -------------------------
PARALLEL = 4  # default parallel workers

# Pipeline files
MUTAUE_81To101_PIPELINE = "./pipeline_mutaue_81To101.json"
ETAUMU_81To101_PIPELINE = "./pipeline_etaumu_81To101.json"
MUTAUE_21To81_PIPELINE = "./pipeline_mutaue_21To81.json"
ETAUMU_21To81_PIPELINE = "./pipeline_etaumu_21To81.json"

# Save path
PARENT_DIR = Path("./normal_cuts")

# Place dummy
MUTAUE_81To101_DIR = None
ETAUMU_81To101_DIR = None
MUTAUE_21To81_DIR = None
ETAUMU_21To81_DIR = None

# Backgrounds
BACKGROUNDS = [
    # ISR samples
    "/work/project/physics/psriling/FCC/FCCee/ISR_zh_ll_ww/ROOT/",
    "/work/project/physics/psriling/FCC/FCCee/ISR_zh_ll_tautau/ROOT/",
    "/work/project/physics/psriling/FCC/FCCee/ISR_zz_ll_tautau/ROOT/",
    "/work/project/physics/psriling/FCC/FCCee/ISR_zww/ROOT/",
    "/work/project/physics/psriling/FCC/FCCee/ISR_vbs/ROOT/"
]
BACKGROUNDS_NAMES = [
    # "HZFourLep",  # Unused
    "zh_ll_ww",
    "zh_ll_tautau",
    "zz_ll_tautau",
    "zww",
    "vbs"
]

# Toggle which channels / types to include (set True/False)
RUN_SIGNALS = True
RUN_BACKGROUNDS = True

MASS_RANGE = range(165, 241, 5)  # for signals

CONFIG = {
    "signal_mutaue_81To101": RUN_SIGNALS,
    "signal_etaumu_81To101": RUN_SIGNALS,
    "signal_mutaue_21To81": RUN_SIGNALS,
    "signal_etaumu_21To81": RUN_SIGNALS,
    "background_mutaue_81To101": RUN_BACKGROUNDS,
    "background_etaumu_81To101": RUN_BACKGROUNDS,
    "background_mutaue_21To81": RUN_BACKGROUNDS,
    "background_etaumu_21To81": RUN_BACKGROUNDS,
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
    name: str = None         # process name

    def out_root(self) -> Path:
        if self.sample_type == "signal":
            base_name = Path(self.input_path).stem if self.name is None else self.name
            return self.output_dir / f"{base_name}.root"
        elif self.sample_type == "background":
            idx = BACKGROUNDS.index(self.input_path)
            bg_name = BACKGROUNDS_NAMES[idx]
            return self.output_dir / f"{self.sample_type}_{bg_name}.root"
        else:
            raise ValueError(f"Unknown sample_type: {self.sample_type}")

    def log_path(self) -> Path:
        # return self.output_dir / f"{self.sample_type}_{Path(self.input_path).name}.txt"
        # Use name if provided
        if self.sample_type == "signal":
            base_name = Path(self.input_path).stem if self.name is None else self.name
            return self.output_dir / f"{base_name}.log"
        elif self.sample_type == "background":
            idx = BACKGROUNDS.index(self.input_path)
            bg_name = BACKGROUNDS_NAMES[idx]
            return self.output_dir / f"{self.sample_type}_{bg_name}.log"
        else:
            raise ValueError(f"Unknown sample_type: {self.sample_type}")

    def command(self) -> List[str]:
        # root expects the analyze_pipeline.cpp(...) as a single argument
        cpp_arg = f'analyze_pipeline.cpp("{self.input_path}","{self.out_root()}","{self.pipeline}")'
        return ["root", "-l", "-b", "-q", cpp_arg]


def ensure_dirs():
    for d in (MUTAUE_81To101_DIR, ETAUMU_81To101_DIR, MUTAUE_21To81_DIR, ETAUMU_21To81_DIR):
        d.mkdir(parents=True, exist_ok=True)


def build_jobs() -> List[Job]:
    jobs: List[Job] = []

    # MuTauE signals
    if CONFIG["signal_mutaue_81To101"]:
        for mass in MASS_RANGE:
            inp = f"/work/project/physics/vwachira/fcc-ee-higgs-lfv/delphes_outputs/mh{mass}_mutau/"
            jobs.append(Job(input_path=inp, output_dir=MUTAUE_81To101_DIR, pipeline=MUTAUE_81To101_PIPELINE,
                            sample_type="signal", channel="mutaue_81To101", name=f"signal_HMuTauE_LFV_{mass}"))
        
    if CONFIG["background_mutaue_81To101"]:
        for bg in BACKGROUNDS:
            jobs.append(Job(input_path=bg, output_dir=MUTAUE_81To101_DIR, pipeline=MUTAUE_81To101_PIPELINE,
                            sample_type="background", channel="mutaue_81To101", name=BACKGROUNDS_NAMES[BACKGROUNDS.index(bg)]))

    # Etaumu signals
    if CONFIG["signal_etaumu_81To101"]:
        for mass in MASS_RANGE:
            inp = f"/work/project/physics/vwachira/fcc-ee-higgs-lfv/delphes_outputs/mh{mass}_etau/"
            jobs.append(Job(input_path=inp, output_dir=ETAUMU_81To101_DIR, pipeline=ETAUMU_81To101_PIPELINE,
                            sample_type="signal", channel="etaumu_81To101", name=f"signal_HETauMu_LFV_{mass}"))
    if CONFIG["background_etaumu_81To101"]:
        for bg in BACKGROUNDS:
            jobs.append(Job(input_path=bg, output_dir=ETAUMU_81To101_DIR, pipeline=ETAUMU_81To101_PIPELINE,
                            sample_type="background", channel="etaumu_81To101", name=BACKGROUNDS_NAMES[BACKGROUNDS.index(bg)]))

    # MuTauE Offshell
    if CONFIG["signal_mutaue_21To81"]:
        for mass in MASS_RANGE:
            inp = f"/work/project/physics/vwachira/fcc-ee-higgs-lfv/delphes_outputs/mh{mass}_mutau/"
            jobs.append(Job(input_path=inp, output_dir=MUTAUE_21To81_DIR, pipeline=MUTAUE_21To81_PIPELINE,
                            sample_type="signal", channel="mutaue_21To81", name=f"signal_HMuTauE_LFV_{mass}"))
    if CONFIG["background_mutaue_21To81"]:
        for bg in BACKGROUNDS:
            jobs.append(Job(input_path=bg, output_dir=MUTAUE_21To81_DIR, pipeline=MUTAUE_21To81_PIPELINE,
                            sample_type="background", channel="mutaue_21To81", name=BACKGROUNDS_NAMES[BACKGROUNDS.index(bg)]))

    # Etaumu Offshell
    if CONFIG["signal_etaumu_21To81"]:
        for mass in MASS_RANGE:
            inp = f"/work/project/physics/vwachira/fcc-ee-higgs-lfv/delphes_outputs/mh{mass}_etau/"
            jobs.append(Job(input_path=inp, output_dir=ETAUMU_21To81_DIR, pipeline=ETAUMU_21To81_PIPELINE,
                            sample_type="signal", channel="etaumu_21To81", name=f"signal_HETauMu_LFV_{mass}"))
    if CONFIG["background_etaumu_21To81"]:
        for bg in BACKGROUNDS:
            jobs.append(Job(input_path=bg, output_dir=ETAUMU_21To81_DIR, pipeline=ETAUMU_21To81_PIPELINE,
                            sample_type="background", channel="etaumu_21To81", name=BACKGROUNDS_NAMES[BACKGROUNDS.index(bg)]))

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

def run_makecard_commands(args, dry_run: bool = False):
    """Builds and runs the makecard.py commands for channels that were processed."""
    print("\n" + "="*30)
    print("Starting makecard generation")
    print("="*30)

    makecard_jobs: Dict[str, Any] = {
        "mutaue_81To101": {
            "in_dir": MUTAUE_81To101_DIR,
            "out_dir": Path(f"{args.output_dir}/{MUTAUE_81To101_DIR.name.replace('hist_', 'datacards_')}"),
            "enabled": CONFIG["signal_mutaue_81To101"] or CONFIG["background_mutaue_81To101"]
        },
        "etaumu_81To101": {
            "in_dir": ETAUMU_81To101_DIR,
            "out_dir": Path(f"{args.output_dir}/{ETAUMU_81To101_DIR.name.replace('hist_', 'datacards_')}"),
            "enabled": CONFIG["signal_etaumu_81To101"] or CONFIG["background_etaumu_81To101"]
        },
        "mutaue_21To81": {
            "in_dir": MUTAUE_21To81_DIR,
            "out_dir": Path(f"{args.output_dir}/{MUTAUE_21To81_DIR.name.replace('hist_', 'datacards_')}"),
            "enabled": CONFIG["signal_mutaue_21To81"] or CONFIG["background_mutaue_21To81"]
        },
        "etaumu_21To81": {
            "in_dir": ETAUMU_21To81_DIR,
            "out_dir": Path(f"{args.output_dir}/{ETAUMU_21To81_DIR.name.replace('hist_', 'datacards_')}"),
            "enabled": CONFIG["signal_etaumu_21To81"] or CONFIG["background_etaumu_21To81"]
        },
    }

    lumi_pb = 5_000_000
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
        
def run_sbatch_commands(args):
    # Copy the datacards/run_limits.py and datacards/slurm_submit.slurm to each output directory
    script_files = ["datacards/run_limits.py", "datacards/slurm_submit.slurm"]
    makecard_jobs: Dict[str, Any] = {
        "mutaue_81To101": {
            "in_dir": MUTAUE_81To101_DIR,
            # "out_dir": Path(f"{args.output_dir}/datacards_{MUTAUE_81To101_DIR.name.replace('hist_', '')}"),
            "out_dir": Path(f"{args.output_dir}/{MUTAUE_81To101_DIR.name.replace('hist_', 'datacards_')}"),
            "enabled": CONFIG["signal_mutaue_81To101"] or CONFIG["background_mutaue_81To101"]
        },
        "etaumu_81To101": {
            "in_dir": ETAUMU_81To101_DIR,
            "out_dir": Path(f"{args.output_dir}/{ETAUMU_81To101_DIR.name.replace('hist_', 'datacards_')}"),
            "enabled": CONFIG["signal_etaumu_81To101"] or CONFIG["background_etaumu_81To101"]
        },
        "mutaue_21To81": {
            "in_dir": MUTAUE_21To81_DIR,
            "out_dir": Path(f"{args.output_dir}/{MUTAUE_21To81_DIR.name.replace('hist_', 'datacards_')}"),
            "enabled": CONFIG["signal_mutaue_21To81"] or CONFIG["background_mutaue_21To81"]
        },
        "etaumu_21To81": {
            "in_dir": ETAUMU_21To81_DIR,
            "out_dir": Path(f"{args.output_dir}/{ETAUMU_21To81_DIR.name.replace('hist_', 'datacards_')}"),
            "enabled": CONFIG["signal_etaumu_21To81"] or CONFIG["background_etaumu_21To81"]
        },
    }
    for name, params in makecard_jobs.items():
        if not params["enabled"]:
            continue
        out_dir = params["out_dir"]
        for script in script_files:
            dest = out_dir / Path(script).name
            try:
                shutil.copy(script, dest)
                print(f"Copied {script} to {dest}")
            except Exception as e:
                print(f"Failed to copy {script} to {dest}: {e}")

    # Copy THEN GO TO EACH DIRECTORY AND SUBMIT
    # THE SCRIPTS MUST BE USED IN THE SAME DIRECTORY
    # for name, params in makecard_jobs.items():
    #     if not params["enabled"]:
    #         continue
    #     out_dir = params["out_dir"]
    #     sbatch_script = out_dir / "slurm_submit.slurm"
    #     if not sbatch_script.exists():
    #         print(f"SBATCH script not found in {out_dir}, skipping sbatch submission.")
    #         continue
        
    #     # USE THE OS.CHDIR TO CHANGE DIRECTORY
    #     current_dir = os.getcwd()
    #     os.chdir(out_dir)
    #     try:
    #         cmd = ["sbatch", "slurm_submit.slurm"]
    #         print(f"Submitting sbatch job in {out_dir}: CMD: {' '.join(cmd)}")
    #         result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    #         print(result.stdout)
    #         if result.stderr:
    #             print("--- STDERR ---")
    #             print(result.stderr)
    #         print(f"SBATCH submission for '{name}' completed successfully.")
    #     except subprocess.CalledProcessError as e:
    #         print(f"ERROR: sbatch submission for '{name}' failed with exit code {e.returncode}.")
    #         print("--- STDOUT ---")
    #         print(e.stdout)
    #         print("--- STDERR ---")
    #         print(e.stderr)
    #     except FileNotFoundError:
    #         print("ERROR: 'sbatch' command not found. Make sure SLURM is installed and configured.")
    #     finally:
    #         os.chdir(current_dir)


def main():
    parser = argparse.ArgumentParser(description="Parallel runner for analyze_pipeline jobs")
    parser.add_argument("--output-dir", "-o", type=str, default=str(PARENT_DIR), help="parent output directory")
    parser.add_argument("--parallel", "-p", type=int, default=PARALLEL, help="number of parallel jobs")
    parser.add_argument("--list", action="store_true", help="only list jobs (don't execute)")
    parser.add_argument("--dry-run", action="store_true", help="show commands without running")
    parser.add_argument("--skip-makecard", action="store_true", help="skip the final makecard step")
    parser.add_argument("--skip-sbatch", action="store_true", help="skip the sbatch submission step after makecard")
    parser.add_argument("--skip-run", action="store_true", help="skip the job execution step")
    args = parser.parse_args()

    print("Output directory set to:", args.output_dir)
    # TOBE CLEAR: OVERRIDE THE GLOBAL PATHS HERE
    global MUTAUE_81To101_DIR, ETAUMU_81To101_DIR, MUTAUE_21To81_DIR, ETAUMU_21To81_DIR
    MUTAUE_81To101_DIR = Path(args.output_dir) / "hist_mutaue_81To101"
    ETAUMU_81To101_DIR = Path(args.output_dir) / "hist_etaumu_81To101"
    MUTAUE_21To81_DIR = Path(args.output_dir) / "hist_mutaue_21To81"
    ETAUMU_21To81_DIR = Path(args.output_dir) / "hist_etaumu_21To81"

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
    if not args.skip_run:
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
        
        
        pipeline_files = [MUTAUE_81To101_PIPELINE, ETAUMU_81To101_PIPELINE, MUTAUE_21To81_PIPELINE, ETAUMU_21To81_PIPELINE]
        for pf in pipeline_files:
            dest = Path(args.output_dir) / Path(pf).name
            try:
                shutil.copy(pf, dest)
                print(f"Copied {pf} to {dest}")
            except Exception as e:
                print(f"Failed to copy {pf} to {dest}: {e}")

    # Copy the datacards/merge_datacards.py to the parent output directory for later use
    merge_script = "datacards/merge_datacards.py"
    dest_merge = Path(args.output_dir) / Path(merge_script).name
    try:
        shutil.copy(merge_script, dest_merge)
        print(f"Copied {merge_script} to {dest_merge}")
    except Exception as e:
        print(f"Failed to copy {merge_script} to {dest_merge}: {e}")

    if failures:
        print(f"\n{len(failures)} jobs failed. Check logs.")
        print("Skipping makecard step due to failures.")
        sys.exit(2)
    else:
        print("\nAll processing completed successfully.")
        if not args.skip_makecard:
            run_makecard_commands(args)
        else:
            print("\nSkipping makecard step as requested.")
        
        run_sbatch_commands(args)
        
        # Run the merge_datacards.py with the option --submit if args.skip_sbatch is False
        cmd = ["python3", "merge_datacards.py"]
        if not args.skip_sbatch:
            cmd.append("--submit")
        print("\nRunning merge_datacards.py with command:", " ".join(cmd))
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=args.output_dir)
            print(result.stdout)
            if result.stderr:
                print("--- STDERR ---")
                print(result.stderr)
            print("merge_datacards.py completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"ERROR: merge_datacards.py failed with exit code {e.returncode}.")
            print("--- STDOUT ---")
            print(e.stdout)
            print("--- STDERR ---")
            print(e.stderr)
            sys.exit(4)
        except FileNotFoundError:
            print("ERROR: 'merge_datacards.py' not found. Make sure it is in the output directory.")
            sys.exit(4)

if __name__ == "__main__":
    main()