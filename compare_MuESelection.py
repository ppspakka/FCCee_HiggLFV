#!/usr/bin/env python3
#SBATCH --job-name=mue_sel
#SBATCH --qos=cu_hpc
#SBATCH --partition=cpugpu
#SBATCH --ntasks=1
#SBATCH --output=sbatch_mue_sel.log
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G

import os
import subprocess
import concurrent.futures

# 1. Define Base Directories
ROOTFiles_DIR = "/work/project/escience/ruttho/FCC-ee_SimpleDelphesAnalysis/EventSample/ISR_HEMu_LFV/"
PPor_Dir = "/work/project/physics/psriling/FCC/FCCee/"
Ptop_Dir = "/work/project/physics/vwachira/fcc-ee-higgs-lfv/delphes_outputs/"

# 2. Define Signal and Background Arrays
# SIGNAL_DIRS = [
#     f"{ROOTFiles_DIR}Hmass110/ROOT/", f"{ROOTFiles_DIR}Hmass115/ROOT/", f"{ROOTFiles_DIR}Hmass120/ROOT/",
#     f"{ROOTFiles_DIR}Hmass125/ROOT/", f"{ROOTFiles_DIR}Hmass130/ROOT/", f"{ROOTFiles_DIR}Hmass135/ROOT/",
#     f"{ROOTFiles_DIR}Hmass140/ROOT/", f"{ROOTFiles_DIR}Hmass145/ROOT/", f"{ROOTFiles_DIR}Hmass150/ROOT/",
#     f"{ROOTFiles_DIR}Hmass155/ROOT/", f"{ROOTFiles_DIR}Hmass160/ROOT/"
# ]

BACKGROUND_DIRS = [
    f"{PPor_Dir}ISR_vbs/ROOT/", f"{PPor_Dir}ISR_zh_ll_tautau/ROOT/", f"{PPor_Dir}ISR_zh_ll_ww/ROOT/",
    f"{PPor_Dir}ISR_zz_ll_tautau/ROOT/", f"{PPor_Dir}ISR_zww/ROOT/"
]

# FILELIST = SIGNAL_DIRS + BACKGROUND_DIRS

BACKGROUND_NAMES = ["ISR_vbs", "ISR_zh_ll_tautau", "ISR_zh_ll_ww", "ISR_zz_ll_tautau", "ISR_zww"]

# 3. Environment Setup String
env_setup = """
source /work/app/cms/cmsset_default.sh
cd ~/binary/CMSSW_14_1_0_pre5/src && source /work/app/share_env/hepsw.sh && eval `scramv1 runtime -sh` && cd -
export ROOT_INCLUDE_PATH="/work/app/delphes/src/Delphes-3.5.0/external/:/work/app/delphes/src/Delphes-3.5.0/:${ROOT_INCLUDE_PATH}"
export LD_LIBRARY_PATH="/work/app/delphes/src/Delphes-3.5.0/:${LD_LIBRARY_PATH}"
"""

# 4. Define the worker function that processes a single job
# def run_analysis_21T81(file_dir, res_name):
#     print(f"--> Starting job: {res_name}")
    
#     # Prepare the ROOT C++ command
#     # root_cmd = f'root -l -b -q "../Analysis_Chain.cpp(\\"{file_dir}\\", \\"{res_name}_MassHeatMap.root\\")"'
#     # Assuming you have a variable for your config file
#     config_json = "../pipeline_mue_21To81.json"
#     root_cmd = f'root -l -b -q "../analyze_pipeline.cpp(\\"{file_dir}\\", \\"{res_name}_21To81Hist.root\\", \\"{config_json}\\")"'
#     full_cmd = f"{env_setup}\n{root_cmd}"
    
#     # Open the log file using a context manager (the 'with' block ensures it closes safely)
#     with open(f"{res_name}_21To81log.txt", "w") as log_file:
#         # subprocess.run blocks *this specific thread* until the ROOT script finishes
#         subprocess.run(full_cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT, executable="/bin/bash")
        
#     print(f"<-- Finished job: {res_name}")
#     return res_name # Return the name so we can track completion

def run_analysis(file_dir, res_name, phsp, suffix):
    print(f"--> Starting job: {res_name}")
    
    # Prepare the ROOT C++ command
    config_json = f"../../pipeline_mue_{phsp}{suffix}.json"
    root_cmd = f'root -l -b -q "../../analyze_pipeline.cpp(\\"{file_dir}\\", \\"{res_name}_{phsp}Hist.root\\", \\"{config_json}\\")"'
    full_cmd = f"{env_setup}\n{root_cmd}"
    
    # Open the log file using a context manager (the 'with' block ensures it closes safely)
    with open(f"{res_name}_{phsp}log.txt", "w") as log_file:
        # subprocess.run blocks *this specific thread* until the ROOT script finishes
        subprocess.run(full_cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT, executable="/bin/bash")
        
    print(f"<-- Finished job: {res_name}")
    return res_name # Return the name so we can track completion

# 5. Setup directories and execution
if __name__ == "__main__":
    output_dir1 = "MasstoPT_Selection"
    output_dir2 = "PTtoMass_Selection"
    os.makedirs("MuE_Selection", exist_ok=True)
    os.chdir("MuE_Selection")
    os.makedirs(output_dir1, exist_ok=True)
    os.makedirs(output_dir2, exist_ok=True)

    # Set to 10 to match your #SBATCH --cpus-per-task=10
    MAX_WORKERS = 10 

    print(f"Submitting jobs to a pool of {MAX_WORKERS} workers...")

    # 6. Create the thread pool and map the tasks
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # Submit all jobs to the executor
        # executor.submit takes the function, followed by its arguments
        futures = []

        # ------------ Mass to PT jobs ------------
        os.chdir(output_dir1) # Start in the first output directory for the 21To81 selection
        # Backgroud jobs
        for f_dir, r_name in zip(BACKGROUND_DIRS, BACKGROUND_NAMES):
            future = executor.submit(run_analysis, f_dir, r_name, "21To81", "")
            future = executor.submit(run_analysis, f_dir, r_name, "81To101", "")
            futures.append(future)
        
        # Higgstranlung jobs
        for i in range(110, 165, 5):
            res_name = f"HLFV_{i}GeV"
            f_dir = f"{ROOTFiles_DIR}Hmass{i}/ROOT/"
            future = executor.submit(run_analysis, f_dir, res_name, "21To81", "")
            future = executor.submit(run_analysis, f_dir, res_name, "81To101", "")
            futures.append(future)

        # VBF jobs
        for j in range(110, 240, 5):
            res_name = f"VBF_HLFV_{j}GeV"
            f_dir = f"{Ptop_Dir}mh{j}_emu/"
            future = executor.submit(run_analysis, f_dir, res_name, "21To81", "")
            future = executor.submit(run_analysis, f_dir, res_name, "81To101", "")
            futures.append(future)

        # ------------ PT to Mass jobs ------------
        os.chdir("../" + output_dir2) # Switch to the second output directory for the 81To101 selection
        # Backgroud jobs
        for f_dir, r_name in zip(BACKGROUND_DIRS, BACKGROUND_NAMES):
            future = executor.submit(run_analysis, f_dir, r_name, "21To81", "_PTthenMass")
            future = executor.submit(run_analysis, f_dir, r_name, "81To101", "_PTthenMass")
            futures.append(future)

        # Higgstranlung jobs
        for i in range(110, 165, 5):
            res_name = f"HLFV_{i}GeV"
            f_dir = f"{ROOTFiles_DIR}Hmass{i}/ROOT/"
            future = executor.submit(run_analysis, f_dir, res_name, "21To81", "_PTthenMass")
            future = executor.submit(run_analysis, f_dir, res_name, "81To101", "_PTthenMass")
            futures.append(future)

        # VBF jobs
        for j in range(110, 240, 5):
            res_name = f"VBF_HLFV_{j}GeV"
            f_dir = f"{Ptop_Dir}mh{j}_emu/"
            future = executor.submit(run_analysis, f_dir, res_name, "21To81", "_PTthenMass")
            future = executor.submit(run_analysis, f_dir, res_name, "81To101", "_PTthenMass")
            futures.append(future)

        # as_completed yields futures as soon as they finish, allowing us to monitor progress
        for future in concurrent.futures.as_completed(futures):
            try:
                finished_job = future.result() # This will raise an exception if the job failed
            except Exception as e:
                print(f"A job failed with error: {e}")

    # 

    print("All jobs completed successfully.")