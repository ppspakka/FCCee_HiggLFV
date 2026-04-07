#!/usr/bin/env python3
#SBATCH --job-name=combined_sel
#SBATCH --qos=cu_hpc
#SBATCH --partition=cpugpu
#SBATCH --ntasks=1
#SBATCH --output=sbatch_combined_sel.log
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G           

import os
import subprocess
import concurrent.futures
import json
import sys

# 1. Define Universal Base Directories
PPor_Dir = "/work/project/physics/psriling/FCC/FCCee/"
Ptop_Dir = "/work/project/physics/vwachira/fcc-ee-higgs-lfv/delphes_outputs/"

BACKGROUND_DIRS = [
    f"{PPor_Dir}ISR_vbs/ROOT/", f"{PPor_Dir}ISR_zh_ll_tautau/ROOT/", f"{PPor_Dir}ISR_zh_ll_ww/ROOT/",
    f"{PPor_Dir}ISR_zz_ll_tautau/ROOT/", f"{PPor_Dir}ISR_zww/ROOT/"
]
BACKGROUND_NAMES = ["ISR_vbs", "ISR_zh_ll_tautau", "ISR_zh_ll_ww", "ISR_zz_ll_tautau", "ISR_zww"]

# 2. Define Channel-Specific Settings
CHANNELS = {
    # "MuE": {
    #     "output_dir": "MuE",
    #     "signal_dir": "/work/project/escience/ruttho/FCC-ee_SimpleDelphesAnalysis/EventSample/ISR_HEMu_LFV/",
    #     "hz_range": range(110, 225, 5),   
    #     "vbf_range": range(110, 240, 5),  
    #     "vbf_suffix": "_emu",
    #     "config_prefix": "pipeline_mue_"
    # },
    "ETauMu": {
        "output_dir": "ETauMu",
        "signal_dir": "/work/project/physics/psriling/FCC/FCCee/ISR_HETauMu_LFV/",
        "hz_range": range(110, 225, 5),   
        "vbf_range": range(110, 245, 5),  
        "vbf_suffix": "_etau",
        "config_prefix": "pipeline_etaumu_"
    }
}

# 3. Environment Extraction (Runs ONCE)
def get_cmssw_env():
    print("Pre-loading CMSSW and ROOT environment...")
    env_setup = """
    source /work/app/cms/cmsset_default.sh
    cd ~/binary/CMSSW_14_1_0_pre5/src && source /work/app/share_env/hepsw.sh && cmsenv && cd -
    export ROOT_INCLUDE_PATH="/work/app/delphes/src/Delphes-3.5.0/external/:/work/app/delphes/src/Delphes-3.5.0/:${ROOT_INCLUDE_PATH}"
    export LD_LIBRARY_PATH="/work/app/delphes/src/Delphes-3.5.0/:${LD_LIBRARY_PATH}"
    python3 -c "import os, json; print(json.dumps(dict(os.environ)))"
    """
    result = subprocess.run(env_setup, shell=True, stdout=subprocess.PIPE, text=True, executable="/bin/bash")
    for line in reversed(result.stdout.splitlines()):
        if line.startswith('{'):
            return json.loads(line)
    raise RuntimeError("Failed to parse environment variables.")

# 4. Pre-Compile the C++ Script (Runs ONCE)
def compile_macro(env_dict):
    print("Compiling analyze_pipeline.cpp into a shared library (.so)...")
    # The command uses ROOT's ACLiC (+ flag) to compile without executing
    compile_cmd = ["root", "-l", "-b", "-q", "-e", ".L analyze_pipeline.cpp+"]
    
    try:
        subprocess.run(compile_cmd, env=env_dict, check=True)
        print("Compilation successful!")
    except subprocess.CalledProcessError:
        print("CRITICAL ERROR: Failed to compile analyze_pipeline.cpp. Exiting.")
        sys.exit(1)

# 5. Define the worker function that processes a single job safely
def run_analysis(file_dir, res_name, phsp, suffix, target_dir, config_prefix, env_dict):
    print(f"--> Starting job: {res_name} in {target_dir}")
    
    config_json = f"../../../../{config_prefix}{phsp}{suffix}.json"
    HIST_DIR = os.path.join(target_dir, "HIST")
    LOG_DIR = os.path.join(target_dir, "LOG")
    os.makedirs(HIST_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    hist_path = os.path.join(HIST_DIR, f"{res_name}_{phsp}Hist.root")
    
    # Notice the '+' after .cpp. Because we pre-compiled it, ROOT will 
    # see the '+' and instantly load the .so binary instead of reading the text file.
    root_cmd = [
        "root", "-l", "-b", "-q", 
        f"../../../../analyze_pipeline.cpp+(\"{file_dir}\", \"{hist_path}\", \"{config_json}\")"
    ]
    
    log_path = os.path.join(LOG_DIR, f"{res_name}_{phsp}log.txt")
    with open(log_path, "w") as log_file:
        subprocess.run(root_cmd, stdout=log_file, stderr=subprocess.STDOUT, cwd=target_dir, env=env_dict)
        
    print(f"<-- Finished job: {res_name}")
    return res_name

# 6. Execution
if __name__ == "__main__":
    MAX_WORKERS = int(os.environ.get("SLURM_CPUS_PER_TASK", 20))
    print(f"Submitting jobs to a pool of {MAX_WORKERS} workers...")

    # Load the environment and compile the macro sequentially before threading begins
    cached_env = get_cmssw_env()
    compile_macro(cached_env)

    # Create the thread pool
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = []

        for channel_name, settings in CHANNELS.items():
            
            base_dir = os.path.join(os.getcwd(), "Results", "Selections")
            New_base_dir = os.path.join(base_dir, settings["output_dir"])
            path_dir1 = os.path.join(New_base_dir, "MasstoPT")
            path_dir2 = os.path.join(New_base_dir, "PTtoMass")
            path_dir3 = os.path.join(New_base_dir, "OnlyPT")
            path_dir4 = os.path.join(New_base_dir, "NoPTMin")
            
            os.makedirs(path_dir1, exist_ok=True)
            os.makedirs(path_dir2, exist_ok=True)
            os.makedirs(path_dir3, exist_ok=True)
            os.makedirs(path_dir4, exist_ok=True)


            datasets = []
            
            for f_dir, r_name in zip(BACKGROUND_DIRS, BACKGROUND_NAMES):
                datasets.append((f_dir, r_name))
                
            for i in settings["hz_range"]:
                datasets.append((f"{settings['signal_dir']}Hmass{i}/ROOT/", f"HLFV_{i}GeV"))
                
            for j in settings["vbf_range"]:
                datasets.append((f"{Ptop_Dir}mh{j}{settings['vbf_suffix']}/", f"VBF_HLFV_{j}GeV"))

            scenarios = [
                {"dir": path_dir1, "suffix": ""},
                {"dir": path_dir2, "suffix": "_PTthenMass"},
                {"dir": path_dir3, "suffix": "_OnlyPT"},
                {"dir": path_dir4, "suffix": "_NoPTMin"}
            ]

            for scenario in scenarios:
                for f_dir, res_name in datasets:
                    futures.append(executor.submit(
                        run_analysis, f_dir, res_name, "21To81", scenario["suffix"], scenario["dir"], settings["config_prefix"], cached_env
                    ))
                    futures.append(executor.submit(
                        run_analysis, f_dir, res_name, "81To101", scenario["suffix"], scenario["dir"], settings["config_prefix"], cached_env
                    ))

        for future in concurrent.futures.as_completed(futures):
            try:
                finished_job = future.result() 
            except Exception as e:
                print(f"A job failed with error: {e}")

    print("All combined jobs completed successfully.")