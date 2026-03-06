"""
template for slurm submission script

[psriling@frontend-03 FCCee_HiggLFV]$ more slurm_submit.slurm
#!/bin/bash

#SBATCH --qos=cu_htc
#SBATCH --partition=cpugpu
#SBATCH --job-name=Limits
#SBATCH --output=sbatch.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=32G

<python script here>


Content of pipeline json config file (etaumu, lowmass)
{
  "parameters": {
    "_comment": "Configuration parameters for the selection cuts",

    "_comment1": "==== Lepton cuts ====",
    "lepton_pt_min": 5.0,

    "_comment2": "==== Z -> ll cuts ====",
    "z_mass": 91,
    "zl_pt_min": 5.0,
    "z_mass_window_upper": 5.0,
    "z_mass_window_lower": 5.0,

    "_comment3": "==== H -> e tau_mu cut ====",
    "mu_pt_min": 0.0,
    "e_pt_min": 40.0,

    "_comment4": "==== MET cuts ====",
    "max_dphi_e_met": 0.1,
    "max_dphi_mu_met": 0.1
  },
  "selections": [
    { "name": "lepton_selection",       "enabled": true },

    { "name": "z_candidate_selection",  "enabled": true },
    { "name": "Z_to_ll",                "enabled": false },

    { "name": "H_to_etau_mu",           "enabled": true },

    { "name": "MET_mu_dphi",            "enabled": true },
    { "name": "finalstate_nocut",       "enabled": true }
  ]
}


2 steps:

1) optimization phase:
- edit the temporary json config file with new cut values
- write the temporary slrum submission script to call the analysis with the temp json config file,
execution: python3 -u runall.py -p <N_cpu> -o optimize_params/<output_dir>
output_dir: <group>_<variable>_<value>
group: mutaue_lowmass, mutaue_highmass, etaumu_lowmass, etaumu_highmass
variable: the cut variable being optimized
value: the cut value being optimized
NOTE: in runall.py, modify the target channel (mutaue, etaumu correspondingly) and mass region (lowmass, highmass)

Example of master switch in runall.py: (offshell group = highmass)
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
For lowmass optimization, set the offshell groups to False and only one channel (mutaue or etaumu) to True
For highmass optimization, set the onshell groups to False and only one channel (mutaue or etaumu) to True


- optimize seperatly for each channel
- there will be total 4 branches of optimization (mutaue_lowmass, mutaue_highmass, etaumu_lowmass, etaumu_highmass)
*lowmass = mH 110-145 GeV, highmass = mH >145 GeV

after submission, the runall.py will run the selections and submit the upper limit calculation jobs to slurm


2) collect results phase: *this can be implemented in the same script (using --mode option) or a separate script
- after all jobs are done, collect the limit results from the output directories
output path:
- etaumu_lowmass: optimize_params/<variable>_<value>/datacards_etaumu/limits.json
- etaumu_highmass: optimize_params/<variable>_<value>/datacards_etaumu_Zoffshell/limits.json
- mutaue_lowmass: optimize_params/<variable>_<value>/datacards_mutaue/limits.json
- mutaue_highmass: optimize_params/<variable>_<value>/datacards_mutaue_Zoffshell/limits.json

# Structure of limits.json
{
  "HETauMu_LFV_110": {
    "0.025": 0.255903,
    "0.16": 0.264165,
    "0.5": 0.276002,
    "0.84": 0.374174,
    "0.975": 0.548557,
    "mass": 110
  },
  "HETauMu_LFV_115": {
    "0.025": 0.231562,
    "0.16": 0.235147,
    "0.5": 0.241444,
    "0.84": 0.337093,
    "0.975": 0.479649,
    "mass": 115
  },
  ... up to mass 145 GeV (lowmass), highmass 150-160 GeV (step size 5 GeV for all)
NOTE: to save the runtime, only the median expected limit (0.5) will be considered for optimization
"""

import argparse
import csv
import datetime as _dt
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# Example cut values to optimize over
# assume each cut is independent, only modify one cut at a tim
# Total iterations = N_cuts1 + N_cuts2 + ... + N_cutsN
scan_parameters = {
    "etaumu_lowmass": { # e hard, mu soft
        "lepton_pt_min": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        "mu_pt_min": [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
        "e_pt_min": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0],
        "max_dphi_mu_met": [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    },
    "etaumu_highmass": {
        "lepton_pt_min": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        "mu_pt_min": [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
        "e_pt_min": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0],
        "max_dphi_mu_met": [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    },
    "mutaue_lowmass": { # mu hard e soft
        "lepton_pt_min": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        "mu_pt_min": [0.0, 10, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0],
        "e_pt_min": [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
        "max_dphi_e_met": [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    },
    "mutaue_highmass": {
        "lepton_pt_min": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        "mu_pt_min": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0],
        "e_pt_min": [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
        "max_dphi_e_met": [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    },
}
REPO_ROOT = Path(__file__).resolve().parent


GROUP_TO_PIPELINE = {
    "etaumu_lowmass": "pipeline_etaumu.json",
    "etaumu_highmass": "pipeline_etaumu_Zoffshell.json",
    "mutaue_lowmass": "pipeline_mutaue.json",
    "mutaue_highmass": "pipeline_mutaue_Zoffshell.json",
}

GROUP_TO_DATACARDS_DIRNAME = {
    "etaumu_lowmass": "datacards_etaumu",
    "etaumu_highmass": "datacards_etaumu_Zoffshell",
    "mutaue_lowmass": "datacards_mutaue",
    "mutaue_highmass": "datacards_mutaue_Zoffshell",
}


ESSENTIAL_ITEMS = [
    "analyze_pipeline.cpp",
    "Delphes.C",
    "Delphes.h",
    "makecard.py",
    "include",
    "src",
    "datacards",
]


@dataclass
class ScanPoint:
    tag: str
    group: str
    param: str
    value: float
    workdir: str
    outdir: str
    pipeline_file: str
    jobid: Optional[str] = None
    submitted_at: Optional[str] = None


def eprint(*args: Any):
    print(*args, file=sys.stderr)


def slugify_value(value: Any) -> str:
    s = str(value)
    s = s.strip()
    s = s.replace("-", "m")
    s = s.replace("+", "p")
    s = s.replace(".", "p")
    s = re.sub(r"[^A-Za-z0-9_]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "val"


def safe_mkdir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def load_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def dump_json(path: Path, data: Dict[str, Any]):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=False)
        f.write("\n")


def ensure_symlink(dst: Path, src: Path):
    if dst.exists() or dst.is_symlink():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.symlink_to(src, target_is_directory=src.is_dir())


def write_text(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def group_to_runall_config(group: str) -> Dict[str, bool]:
    """Return CONFIG dict for runall.py that enables only one channel+mass-region."""
    enabled: Dict[str, bool] = {
        "mutaue_signal": False,
        "etaumu_signal": False,
        "mutaue_offshell_signal": False,
        "etaumu_offshell_signal": False,
        "mutaue_background": False,
        "etaumu_background": False,
        "mutaue_offshell_background": False,
        "etaumu_offshell_background": False,
    }
    if group == "etaumu_lowmass":
        enabled["etaumu_signal"] = True
        enabled["etaumu_background"] = True
    elif group == "etaumu_highmass":
        enabled["etaumu_offshell_signal"] = True
        enabled["etaumu_offshell_background"] = True
    elif group == "mutaue_lowmass":
        enabled["mutaue_signal"] = True
        enabled["mutaue_background"] = True
    elif group == "mutaue_highmass":
        enabled["mutaue_offshell_signal"] = True
        enabled["mutaue_offshell_background"] = True
    else:
        raise KeyError(f"Unknown group: {group}")
    return enabled


def replace_brace_block(lines: List[str], start_idx: int) -> Tuple[int, int]:
    """Given start line containing an opening '{', return (start_idx, end_idx_exclusive) for the matching block."""
    depth = 0
    started = False
    for i in range(start_idx, len(lines)):
        line = lines[i]
        if "{" in line:
            depth += line.count("{")
            started = True
        if "}" in line and started:
            depth -= line.count("}")
            if depth <= 0:
                return start_idx, i + 1
    raise ValueError("Unmatched brace block")


def patch_runall_config(workdir: Path, group: str):
    """Copy runall.py into workdir and patch CONFIG to only run the requested group."""
    src = REPO_ROOT / "runall.py"
    dst = workdir / "runall.py"
    shutil.copy2(src, dst)

    enabled = group_to_runall_config(group)
    lines = read_text(dst).splitlines(True)  # keep newlines
    start = None
    for i, line in enumerate(lines):
        if re.match(r"^\s*CONFIG\s*=\s*\{\s*$", line):
            start = i
            break
    if start is None:
        raise ValueError("Could not find 'CONFIG = {' block in runall.py")
    s, e = replace_brace_block(lines, start)
    indent = re.match(r"^(\s*)", lines[s]).group(1)
    inner = indent + "    "

    new_block = [
        f"{indent}CONFIG = {{\n",
        f'{inner}"mutaue_signal": {enabled["mutaue_signal"]},\n',
        f'{inner}"etaumu_signal": {enabled["etaumu_signal"]},\n',
        "\n",
        f'{inner}"mutaue_offshell_signal": {enabled["mutaue_offshell_signal"]},\n',
        f'{inner}"etaumu_offshell_signal": {enabled["etaumu_offshell_signal"]},\n',
        "\n",
        f'{inner}"mutaue_background": {enabled["mutaue_background"]},\n',
        f'{inner}"etaumu_background": {enabled["etaumu_background"]},\n',
        "\n",
        f'{inner}"mutaue_offshell_background": {enabled["mutaue_offshell_background"]},\n',
        f'{inner}"etaumu_offshell_background": {enabled["etaumu_offshell_background"]},\n',
        f"{indent}}}\n",
    ]

    patched = lines[:s] + new_block + lines[e:]
    write_text(dst, "".join(patched))


def prepare_workdir(workdir: Path, overwrite: bool = False):
    if workdir.exists() and overwrite:
        shutil.rmtree(workdir)
    safe_mkdir(workdir)

    # Always keep an isolated runall.py copy in the workdir (we patch CONFIG per scan point)
    # Note: patched later in submit_mode() after workdir creation.
    shutil.copy2(REPO_ROOT / "runall.py", workdir / "runall.py")

    # Symlink essential project files/dirs
    for rel in ESSENTIAL_ITEMS:
        src = REPO_ROOT / rel
        if not src.exists():
            raise FileNotFoundError(f"Missing required item in repo root: {src}")
        ensure_symlink(workdir / rel, src)

    # Copy all pipeline configs so each scan point is isolated
    for pf in [
        "pipeline_mutaue.json",
        "pipeline_etaumu.json",
        "pipeline_mutaue_Zoffshell.json",
        "pipeline_etaumu_Zoffshell.json",
    ]:
        src = REPO_ROOT / pf
        dst = workdir / pf
        if dst.exists():
            continue
        shutil.copy2(src, dst)


def patch_pipeline_parameter(
    workdir: Path, pipeline_file: str, param: str, value: float
):
    path = workdir / pipeline_file
    if not path.exists():
        raise FileNotFoundError(f"Pipeline file not found in workdir: {path}")
    data = load_json(path)
    if "parameters" not in data or not isinstance(data["parameters"], dict):
        raise ValueError(
            f"Invalid pipeline JSON structure (no 'parameters' dict): {path}"
        )

    if param not in data["parameters"]:
        raise KeyError(
            f"Parameter '{param}' not found in {pipeline_file}. "
            f"Available keys: {sorted(list(data['parameters'].keys()))}"
        )
    data["parameters"][param] = float(value)
    dump_json(path, data)


def read_sbatch_jobid(sbatch_stdout: str) -> Optional[str]:
    m = re.search(r"Submitted batch job\s+(\d+)", sbatch_stdout)
    return m.group(1) if m else None


def build_runall_slurm_script(
    *,
    workdir: Path,
    outdir: Path,
    parallel: int,
    job_name: str,
    cpus: int,
    mem: str,
    partition: Optional[str],
    qos: Optional[str],
    time_limit: Optional[str],
    extra_sbatch: List[str],
    python_exe: str,
) -> str:
    lines: List[str] = [
        "#!/bin/bash",
        "",
        f"#SBATCH --job-name={job_name}",
        "#SBATCH --output=sbatch_runall.log",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --mem={mem}",
    ]
    if partition:
        lines.append(f"#SBATCH --partition={partition}")
    if qos:
        lines.append(f"#SBATCH --qos={qos}")
    if time_limit:
        lines.append(f"#SBATCH --time={time_limit}")
    for item in extra_sbatch:
        item = item.strip()
        if not item:
            continue
        if not item.startswith("#SBATCH"):
            item = "#SBATCH " + item
        lines.append(item)

    cmd = f'{python_exe} -u runall.py -p {parallel} -o "{outdir}"'
    lines += [
        "",
        "set -euo pipefail",
        f'cd "{workdir}"',
        'echo "PWD=$(pwd)"',
        'echo "HOST=$(hostname)"',
        f'echo "CMD={cmd}"',
        cmd,
        "",
    ]
    return "\n".join(lines)


def load_scan_values(group: str, param: str) -> List[float]:
    if group not in scan_parameters:
        raise KeyError(
            f"Unknown group '{group}'. Available: {sorted(list(scan_parameters.keys()))}"
        )
    params = scan_parameters[group]
    if param not in params:
        raise KeyError(
            f"Unknown param '{param}' for group '{group}'. Available: {sorted(list(params.keys()))}"
        )
    return [float(v) for v in params[param]]


def write_manifest(manifest_path: Path, points: List[ScanPoint]):
    safe_mkdir(manifest_path.parent)
    # append-friendly JSONL (one json per line) to avoid read-modify-write races
    with open(manifest_path, "a", encoding="utf-8") as f:
        for p in points:
            f.write(json.dumps(asdict(p), sort_keys=False) + "\n")


def read_manifest(manifest_path: Path) -> List[ScanPoint]:
    if not manifest_path.exists():
        return []
    # manifest is JSONL append-only; allow multiple entries per tag.
    # Keep only the latest occurrence so reruns with --overwrite don't double count.
    latest_by_tag: Dict[str, ScanPoint] = {}
    with open(manifest_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            p = ScanPoint(**obj)
            if p.tag in latest_by_tag:
                # re-insert to preserve "latest" ordering
                del latest_by_tag[p.tag]
            latest_by_tag[p.tag] = p
    return list(latest_by_tag.values())


def median_expected_from_limits_json(path: Path) -> Tuple[Optional[float], int, int]:
    """Return (metric_value, n_ok, n_total) for quantile 0.5 across masses."""
    data = load_json(path)
    values: List[float] = []
    n_total = 0
    for _, entry in data.items():
        if not isinstance(entry, dict):
            continue
        n_total += 1
        v = entry.get("0.5")
        if v is None:
            continue
        try:
            values.append(float(v))
        except Exception:
            continue
    if not values:
        return None, 0, n_total
    # Default metric is mean; selection happens in caller
    return float(sum(values) / len(values)), len(values), n_total


def metric_reduce(values: List[float], metric: str) -> float:
    if not values:
        raise ValueError("No values")
    if metric == "mean":
        return float(sum(values) / len(values))
    if metric == "median":
        vs = sorted(values)
        mid = len(vs) // 2
        if len(vs) % 2 == 1:
            return float(vs[mid])
        return float((vs[mid - 1] + vs[mid]) / 2.0)
    if metric == "min":
        return float(min(values))
    if metric == "max":
        return float(max(values))
    raise ValueError(f"Unknown metric: {metric}")


def collect_point_limits(
    limits_path: Path, metric: str
) -> Tuple[Optional[float], int, int]:
    data = load_json(limits_path)
    values: List[float] = []
    n_total = 0
    for _, entry in data.items():
        if not isinstance(entry, dict):
            continue
        n_total += 1
        v = entry.get("0.5")
        if v is None:
            continue
        try:
            values.append(float(v))
        except Exception:
            continue
    if not values:
        return None, 0, n_total
    return metric_reduce(values, metric), len(values), n_total


def format_limit_lt(value: Optional[float], precision: int = 6) -> str:
    if value is None:
        return ""
    # Use general format by default for readability, but keep stable precision.
    fmt = f"{{:.{precision}g}}"
    return "<" + fmt.format(float(value))


def parse_limits_json_05(limits_path: Path) -> Tuple[List[Dict[str, Any]], int, int]:
    """Parse limits.json and return (entries, n_ok, n_total).

    entries items:
      { process, mass, limit_0p5, limit_str }
    """
    data = load_json(limits_path)
    entries: List[Dict[str, Any]] = []
    n_total = 0
    n_ok = 0

    for proc, entry in data.items():
        if not isinstance(entry, dict):
            continue
        n_total += 1

        mass_val = entry.get("mass")
        if mass_val is None:
            m = re.search(r"_(\d+)$", str(proc))
            mass_val = int(m.group(1)) if m else None

        limit_val = entry.get("0.5")
        if limit_val is None:
            limit_f = None
        else:
            try:
                limit_f = float(limit_val)
            except Exception:
                limit_f = None

        if limit_f is not None:
            n_ok += 1

        entries.append(
            {
                "process": str(proc),
                "mass": int(mass_val) if mass_val is not None else "",
                "limit_0p5": limit_f if limit_f is not None else "",
                "limit_str": format_limit_lt(limit_f),
            }
        )

    # Sort by mass then process name for stable output
    entries.sort(key=lambda x: (x["mass"] if x["mass"] != "" else 10**9, x["process"]))
    return entries, n_ok, n_total


def submit_mode(args: argparse.Namespace) -> int:
    group: str = args.group
    param: str = args.param
    pipeline_file = GROUP_TO_PIPELINE.get(group)
    if not pipeline_file:
        raise SystemExit(f"No pipeline mapping for group '{group}'")
    values = load_scan_values(group, param)

    base_out = Path(args.base_out).resolve()
    work_base = Path(args.work_base).resolve()
    manifest_path = base_out / args.manifest
    python_exe = args.python

    submitted: List[ScanPoint] = []
    skipped = 0

    for v in values:
        tag = f"{group}_{param}_{slugify_value(v)}"
        outdir = base_out / tag
        workdir = work_base / tag
        slurm_path = workdir / "slurm_runall.slurm"

        if outdir.exists() and not args.overwrite:
            skipped += 1
            continue

        prepare_workdir(workdir, overwrite=args.overwrite)
        patch_runall_config(workdir, group)
        patch_pipeline_parameter(workdir, pipeline_file, param, v)
        safe_mkdir(outdir)

        slurm_script = build_runall_slurm_script(
            workdir=workdir,
            outdir=outdir,
            parallel=args.parallel,
            job_name=tag[:128],
            cpus=args.cpus,
            mem=args.mem,
            partition=args.partition,
            qos=args.qos,
            time_limit=args.time,
            extra_sbatch=args.sbatch_extra,
            python_exe=python_exe,
        )
        slurm_path.write_text(slurm_script, encoding="utf-8")

        point = ScanPoint(
            tag=tag,
            group=group,
            param=param,
            value=float(v),
            workdir=str(workdir),
            outdir=str(outdir),
            pipeline_file=pipeline_file,
        )

        if args.dry_run:
            print(f"[DRY] would submit: {tag}")
            submitted.append(point)
            continue

        try:
            res = subprocess.run(
                ["sbatch", slurm_path.name],
                cwd=str(workdir),
                capture_output=True,
                text=True,
                check=True,
            )
        except FileNotFoundError:
            eprint("ERROR: 'sbatch' not found in PATH. Are you on a SLURM login node?")
            return 2
        except subprocess.CalledProcessError as e:
            eprint(f"ERROR: sbatch failed for {tag}")
            eprint(e.stdout)
            eprint(e.stderr)
            return 3

        jobid = read_sbatch_jobid(res.stdout)
        point.jobid = jobid
        point.submitted_at = _dt.datetime.now(tz=_dt.timezone.utc).isoformat()
        print(f"Submitted {tag} -> jobid={jobid or 'UNKNOWN'}")
        submitted.append(point)

    if submitted:
        write_manifest(manifest_path, submitted)

    print(
        f"Done. Submitted {len(submitted)} points. Skipped {skipped} (already exists; use --overwrite to redo)."
    )
    print(f"Manifest: {manifest_path}")
    return 0

def collect_mode(args: argparse.Namespace) -> int:
    base_out = Path(args.base_out).resolve()
    manifest_path = base_out / args.manifest
    # NOTE: collect mode no longer optimizes by averaging across masses.
    # The --metric arg is kept for backward CLI compatibility, but is not used.

    points = read_manifest(manifest_path)
    if not points:
        # fallback: infer from directories if manifest is missing
        points = []
        for d in sorted(base_out.glob("*")):
            if not d.is_dir():
                continue
            # best-effort parse: <group>_<param>_<value>
            parts = d.name.split("_")
            if len(parts) < 3:
                continue
            group = (
                "_".join(parts[0:2])
                if parts[1] in ("lowmass", "highmass")
                else parts[0]
            )
            # if group was joined above, param starts after 2nd token
            if group.endswith("lowmass") or group.endswith("highmass"):
                if len(parts) < 4:
                    continue
                param = parts[2]
            else:
                param = parts[1]
            points.append(
                ScanPoint(
                    tag=d.name,
                    group=group,
                    param=param,
                    value=float("nan"),
                    workdir="",
                    outdir=str(d),
                    pipeline_file=GROUP_TO_PIPELINE.get(group, ""),
                )
            )

    # filter
    if args.group:
        points = [p for p in points if p.group == args.group]
    if args.param:
        points = [p for p in points if p.param == args.param]
    if not points:
        eprint("No scan points to collect.")
        return 2

    rows: List[Dict[str, Any]] = []

    for p in points:
        outdir = Path(p.outdir)
        dc = GROUP_TO_DATACARDS_DIRNAME.get(p.group)
        if not dc:
            continue
        limits_path = outdir / dc / "limits.json"
        if not limits_path.exists():
            rows.append(
                {
                    "tag": p.tag,
                    "group": p.group,
                    "param": p.param,
                    "value": p.value,
                    "jobid": p.jobid or "",
                    "point_dir": str(outdir),
                    "limits_json": str(limits_path),
                    "done": False,
                    "n_ok": 0,
                    "n_total": 0,
                    "process": "",
                    "mass_point": "",
                    "limit_0p5": "",
                    "limit_str": "",
                }
            )
            continue

        entries, n_ok, n_total = parse_limits_json_05(limits_path)
        any_done = any(e["limit_0p5"] != "" for e in entries)
        if not entries:
            rows.append(
                {
                    "tag": p.tag,
                    "group": p.group,
                    "param": p.param,
                    "value": p.value,
                    "jobid": p.jobid or "",
                    "point_dir": str(outdir),
                    "limits_json": str(limits_path),
                    "done": False,
                    "n_ok": 0,
                    "n_total": 0,
                    "process": "",
                    "mass_point": "",
                    "limit_0p5": "",
                    "limit_str": "",
                }
            )
            continue

        for e in entries:
            rows.append(
                {
                    "tag": p.tag,
                    "group": p.group,
                    "param": p.param,
                    "value": p.value,
                    "jobid": p.jobid or "",
                    "point_dir": str(outdir),
                    "limits_json": str(limits_path),
                    "done": any_done,
                    "n_ok": n_ok,
                    "n_total": n_total,
                    "process": e["process"],
                    "mass_point": e["mass"],
                    "limit_0p5": e["limit_0p5"],
                    "limit_str": e["limit_str"],
                }
            )

    # Write CSV summary
    safe_mkdir(base_out)
    summary_csv = base_out / (args.summary or "summary.csv")
    with open(summary_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "tag",
                "group",
                "param",
                "value",
                "jobid",
                "done",
                "n_ok",
                "n_total",
                "point_dir",
                "limits_json",
                "mass_point",
                "process",
                "limit_0p5",
                "limit_str",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    n_done = sum(1 for r in rows if r["done"])
    print(f"Wrote: {summary_csv}")
    print(f"Collected {n_done}/{len(rows)} points with completed limits.")

    # -------- Pretty report (per-mass best + full per-mass tables) --------
    done_rows = [r for r in rows if r.get("limit_0p5") not in ("", None)]
    masses = sorted({int(r["mass_point"]) for r in done_rows if str(r["mass_point"]).isdigit()})
    if masses:
        print("\n=== Best per mass (expected 0.5) ===")
        for m in masses:
            candidates = [r for r in done_rows if r["mass_point"] == m]
            if not candidates:
                continue
            # pick smallest limit
            best_r = min(candidates, key=lambda x: float(x["limit_0p5"]))
            print(
                f"Mass {m}: best {best_r['param']}={best_r['value']}  limit {best_r['limit_str']}  ({best_r['tag']})"
            )

        print("\n=== Full details per mass (all scanned points) ===")
        for m in masses:
            candidates = [r for r in done_rows if r["mass_point"] == m]
            if not candidates:
                continue
            # sort by parameter value if numeric, else by limit
            def _sort_key(x: Dict[str, Any]):
                try:
                    return (0, float(x["value"]))
                except Exception:
                    return (1, float(x["limit_0p5"]))

            candidates.sort(key=_sort_key)
            print(f"\nM{m}")
            for r in candidates:
                print(f"  {r['param']}={r['value']}: limit {r['limit_str']}  [{r['tag']}]  ({r['process']})")

    if args.prune_workdirs and points:
        pruned = 0
        done_tags = {r["tag"] for r in rows if r["done"]}
        for p in points:
            if p.tag not in done_tags:
                continue
            if not p.workdir:
                continue
            wd = Path(p.workdir)
            if wd.exists():
                shutil.rmtree(wd)
                pruned += 1
        print(f"Pruned {pruned} workdir(s) under optimize_work.")

    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Brute-force optimization wrapper: modify pipeline JSON, submit runall.py via slurm, then collect limits.json",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--mode",
        choices=["submit", "collect"],
        required=True,
        help="submit scan jobs or collect limits",
    )
    parser.add_argument(
        "--base-out",
        default=str(REPO_ROOT / "optimize_params"),
        help="base output directory for scan points",
    )
    parser.add_argument(
        "--manifest",
        default="manifest.jsonl",
        help="manifest file (relative to --base-out)",
    )
    parser.add_argument(
        "--group", choices=sorted(list(scan_parameters.keys())), help="scan group"
    )
    parser.add_argument(
        "--param", help="parameter name to scan (key in pipeline 'parameters')"
    )

    # submit-specific
    parser.add_argument(
        "--work-base",
        default=str(REPO_ROOT / "optimize_work"),
        help="base working directory (isolated per point)",
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=8,
        help="runall.py -p (thread pool for ROOT jobs)",
    )
    parser.add_argument(
        "--python", default="python3", help="python executable used inside slurm job"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="overwrite existing point directories"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="generate workdirs/scripts but do not sbatch",
    )

    # slurm resources for runall job
    parser.add_argument("--cpus", type=int, default=8, help="slurm --cpus-per-task")
    parser.add_argument("--mem", default="8G", help="slurm --mem")
    parser.add_argument("--partition", default="cpugpu", help="slurm --partition")
    parser.add_argument("--qos", default="cu_htc", help="slurm --qos")
    parser.add_argument("--time", default=None, help="slurm --time (e.g. 02:00:00)")
    parser.add_argument(
        "--sbatch-extra",
        nargs="*",
        default=[],
        help="extra SBATCH options, e.g. --sbatch-extra --account=myacct --constraint=...",
    )

    # collect-specific
    parser.add_argument(
        "--metric",
        choices=["mean", "median", "min", "max"],
        default="mean",
        help="reduce expected(0.5) over masses",
    )
    parser.add_argument(
        "--summary", default=None, help="output CSV name (relative to --base-out)"
    )
    parser.add_argument(
        "--prune-workdirs",
        action="store_true",
        help="delete workdirs for points with completed limits.json",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.mode == "submit":
        if not args.group or not args.param:
            parser.error("--mode submit requires --group and --param")
        return submit_mode(args)

    if args.mode == "collect":
        # group/param optional filters
        return collect_mode(args)

    parser.error(f"Unknown mode: {args.mode}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

'''

example execution script.sh:

# lepton_pt_min
python3 optimize_results.py --mode submit --group mutaue_lowmass --param lepton_pt_min
python3 optimize_results.py --mode submit --group etaumu_lowmass --param lepton_pt_min
python3 optimize_results.py --mode submit --group mutaue_highmass --param lepton_pt_min
python3 optimize_results.py --mode submit --group etaumu_highmass --param lepton_pt_min

# deltaPhi e-met
python3 optimize_results.py --mode submit --group mutaue_lowmass --param max_dphi_e_met
python3 optimize_results.py --mode submit --group mutaue_highmass --param max_dphi_e_met

# deltaPhi mu-met
python3 optimize_results.py --mode submit --group etaumu_lowmass --param max_dphi_mu_met
python3 optimize_results.py --mode submit --group etaumu_highmass --param max_dphi_mu_met

# mu_pt_min
python3 optimize_results.py --mode submit --group mutaue_lowmass --param mu_pt_min
python3 optimize_results.py --mode submit --group etaumu_lowmass --param mu_pt_min
python3 optimize_results.py --mode submit --group mutaue_highmass --param mu_pt_min
python3 optimize_results.py --mode submit --group etaumu_highmass --param mu_pt_min

# e_pt_min
python3 optimize_results.py --mode submit --group mutaue_lowmass --param e_pt_min
python3 optimize_results.py --mode submit --group etaumu_lowmass --param e_pt_min
python3 optimize_results.py --mode submit --group mutaue_highmass --param e_pt_min
python3 optimize_results.py --mode submit --group etaumu_highmass --param e_pt_min


# Collect results
TARGET_DIR=opt_logs
# lepton_pt_min
python3 optimize_results.py --mode collect --group mutaue_lowmass --param lepton_pt_min | tee $TARGET_DIR/opt_mutaue_lowmass_lepton_pt_min.log
python3 optimize_results.py --mode collect --group etaumu_lowmass --param lepton_pt_min | tee $TARGET_DIR/opt_etaumu_lowmass_lepton_pt_min.log
python3 optimize_results.py --mode collect --group mutaue_highmass --param lepton_pt_min | tee $TARGET_DIR/opt_mutaue_highmass_lepton_pt_min.log
python3 optimize_results.py --mode collect --group etaumu_highmass --param lepton_pt_min | tee $TARGET_DIR/opt_etaumu_highmass_lepton_pt_min.log

# deltaPhi e-met
python3 optimize_results.py --mode collect --group mutaue_lowmass --param max_dphi_e_met | tee $TARGET_DIR/opt_mutaue_lowmass_max_dphi_e_met.log
python3 optimize_results.py --mode collect --group mutaue_highmass --param max_dphi_e_met | tee $TARGET_DIR/opt_mutaue_highmass_max_dphi_e_met.log
# deltaPhi mu-met
python3 optimize_results.py --mode collect --group etaumu_lowmass --param max_dphi_mu_met | tee $TARGET_DIR/opt_etaumu_lowmass_max_dphi_mu_met.log
python3 optimize_results.py --mode collect --group etaumu_highmass --param max_dphi_mu_met | tee $TARGET_DIR/opt_etaumu_highmass_max_dphi_mu_met.log

# mu_pt_min
python3 optimize_results.py --mode collect --group mutaue_lowmass --param mu_pt_min | tee $TARGET_DIR/opt_mutaue_lowmass_mu_pt_min.log
python3 optimize_results.py --mode collect --group etaumu_lowmass --param mu_pt_min | tee $TARGET_DIR/opt_etaumu_lowmass_mu_pt_min.log
python3 optimize_results.py --mode collect --group mutaue_highmass --param mu_pt_min | tee $TARGET_DIR/opt_mutaue_highmass_mu_pt_min.log
python3 optimize_results.py --mode collect --group etaumu_highmass --param mu_pt_min | tee $TARGET_DIR/opt_etaumu_highmass_mu_pt_min.log

# e_pt_min
python3 optimize_results.py --mode collect --group mutaue_lowmass --param e_pt_min | tee $TARGET_DIR/opt_mutaue_lowmass_e_pt_min.log
python3 optimize_results.py --mode collect --group etaumu_lowmass --param e_pt_min | tee $TARGET_DIR/opt_etaumu_lowmass_e_pt_min.log
python3 optimize_results.py --mode collect --group mutaue_highmass --param e_pt_min | tee $TARGET_DIR/opt_mutaue_highmass_e_pt_min.log
python3 optimize_results.py --mode collect --group etaumu_highmass --param e_pt_min | tee $TARGET_DIR/opt_etaumu_highmass_e_pt_min.log

'''