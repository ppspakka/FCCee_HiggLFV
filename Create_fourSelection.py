#!/usr/bin/env python3
#SBATCH --job-name=combined_sel
#SBATCH --qos=cu_hpc
#SBATCH --partition=cpugpu
#SBATCH --ntasks=1
#SBATCH --output=sbatch_combined_sel.log
#SBATCH --cpus-per-task=20  # <-- Doubled to 20 CPUs
#SBATCH --mem=128G           # <-- Doubled to 128 GB RAM

import os
import subprocess
import concurrent.futures

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
    "MuE": {
        "output_dir": "MuE_Selection",
        "signal_dir": "/work/project/escience/ruttho/FCC-ee_SimpleDelphesAnalysis/EventSample/ISR_HEMu_LFV/",
        "hz_range": range(110, 165, 5),   # 110 to 160
        "vbf_range": range(110, 240, 5),  # 110 to 235
        "vbf_suffix": "_emu",
        "config_prefix": "pipeline_mue_"
    },
    "ETauMu": {
        "output_dir": "ETauMu_Selection",
        "signal_dir": "/work/project/physics/psriling/FCC/FCCee/ISR_HETauMu_LFV/",
        "hz_range": range(110, 220, 5),   # 110 to 215
        "vbf_range": range(110, 245, 5),  # 110 to 240
        "vbf_suffix": "_etau",
        "config_prefix": "pipeline_etaumu_"
    }
}

# 3. Environment Setup String
env_setup = """
source /work/app/cms/cmsset_default.sh
cd ~/binary/CMSSW_14_1_0_pre5/src && source /work/app/share_env/hepsw.sh && eval `scramv1 runtime -sh` && cd -
export ROOT_INCLUDE_PATH="/work/app/delphes/src/Delphes-3.5.0/external/:/work/app/delphes/src/Delphes-3.5.0/:${ROOT_INCLUDE_PATH}"
export LD_LIBRARY_PATH="/work/app/delphes/src/Delphes-3.5.0/:${LD_LIBRARY_PATH}"
"""

# 4. Define the worker function that processes a single job safely
def run_analysis(file_dir, res_name, phsp, suffix, target_dir, config_prefix):
    print(f"--> Starting job: {res_name} in {target_dir}")
    
    # Pathing: Use the specific config_prefix for the channel being processed
    config_json = f"../../{config_prefix}{phsp}{suffix}.json"
    root_cmd = f'root -l -b -q "../../analyze_pipeline.cpp(\\"{file_dir}\\", \\"{res_name}_{phsp}Hist.root\\", \\"{config_json}\\")"'
    full_cmd = f"{env_setup}\n{root_cmd}"
    
    log_path = os.path.join(target_dir, f"{res_name}_{phsp}log.txt")
    with open(log_path, "w") as log_file:
        subprocess.run(full_cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT, executable="/bin/bash", cwd=target_dir)
        
    print(f"<-- Finished job: {res_name}")
    return res_name

# 5. Execution
if __name__ == "__main__":
    MAX_WORKERS = 40  # <-- Matches the 40 allocated CPUs
    print(f"Submitting jobs to a pool of {MAX_WORKERS} workers...")

    # Create the thread pool
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = []

        # Loop over each channel defined in our dictionary
        for channel_name, settings in CHANNELS.items():
            
            # Create directories for this specific channel
            base_dir = os.path.join(os.getcwd(), settings["output_dir"])
            path_dir1 = os.path.join(base_dir, "MasstoPT_Selection")
            path_dir2 = os.path.join(base_dir, "PTtoMass_Selection")
            
            os.makedirs(path_dir1, exist_ok=True)
            os.makedirs(path_dir2, exist_ok=True)

            # Compile the dataset list for this channel
            datasets = []
            
            # Backgrounds (Shared between both channels)
            for f_dir, r_name in zip(BACKGROUND_DIRS, BACKGROUND_NAMES):
                datasets.append((f_dir, r_name))
                
            # Higgstrahlung Signals (Uses channel-specific range and path)
            for i in settings["hz_range"]:
                datasets.append((f"{settings['signal_dir']}Hmass{i}/ROOT/", f"HLFV_{i}GeV"))
                
            # VBF Signals (Uses channel-specific range and suffix)
            for j in settings["vbf_range"]:
                datasets.append((f"{Ptop_Dir}mh{j}{settings['vbf_suffix']}/", f"VBF_HLFV_{j}GeV"))

            # Define configurations (Standard vs PTthenMass)
            scenarios = [
                {"dir": path_dir1, "suffix": ""},
                {"dir": path_dir2, "suffix": "_PTthenMass"}
            ]

            # Submit tasks to the executor
            for scenario in scenarios:
                for f_dir, res_name in datasets:
                    # Submit 21To81 Phase Space
                    futures.append(executor.submit(
                        run_analysis, f_dir, res_name, "21To81", scenario["suffix"], scenario["dir"], settings["config_prefix"]
                    ))
                    # Submit 81To101 Phase Space
                    futures.append(executor.submit(
                        run_analysis, f_dir, res_name, "81To101", scenario["suffix"], scenario["dir"], settings["config_prefix"]
                    ))

        # Monitor progress as jobs finish
        for future in concurrent.futures.as_completed(futures):
            try:
                finished_job = future.result() 
            except Exception as e:
                print(f"A job failed with error: {e}")

    print("All combined jobs completed successfully.")