#!/usr/bin/env python3
#SBATCH --job-name=combined_sel
#SBATCH --qos=cu_hpc
#SBATCH --partition=cpugpu
#SBATCH --ntasks=1
#SBATCH --output=sbatch_acceptance.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

import sys
import re
import os
import pandas as pd
import ROOT
import json # <-- Imported JSON
from concurrent.futures import ProcessPoolExecutor, as_completed

# --- Global Optimizations ---
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError 

PATTERN_INIT = re.compile(r"00_Initial_n_muons$")
PATTERN_INV = re.compile(r".*_finalstate_nocut_m_h_invariant_count$")

def extract_yields_and_count(file_path, masses, window=10):
    f = ROOT.TFile.Open(file_path, "READ")
    if not f or f.IsZombie():
        return 0.0, {m: 0.0 for m in masses}
    
    hist_init = None
    hist_inv = None
    
    dirs_to_process = [f]
    while dirs_to_process:
        current_dir = dirs_to_process.pop()
        for key in current_dir.GetListOfKeys():
            obj_class = key.GetClassName()
            name = key.GetName()
            
            if obj_class.startswith("TDirectory"):
                dirs_to_process.append(key.ReadObj())
            elif obj_class.startswith("TH1"):
                if not hist_init and PATTERN_INIT.search(name):
                    hist_init = key.ReadObj()
                elif not hist_inv and PATTERN_INV.search(name):
                    hist_inv = key.ReadObj()
                    
            if hist_init and hist_inv:
                break
        if hist_init and hist_inv:
            break

    initial_count = hist_init.Integral() if hist_init else 0.0
    
    yields = {}
    for m in masses:
        if hist_inv:
            axis = hist_inv.GetXaxis()
            bin_low = axis.FindBin(m - window)
            bin_high = axis.FindBin(m + window)
            yields[m] = hist_inv.Integral(bin_low, bin_high)
        else:
            yields[m] = 0.0
            
    f.Close()
    return initial_count, yields

def GetTotalEvt(filepath):
    if not os.path.exists(filepath): return None
    with open(filepath, 'r') as f:
        for line in f:
            if "Initial" in line:
                return int(line.split("(")[0].strip().split()[-1])
    return None

def GetFinalAcceptance(filepath):
    if not os.path.exists(filepath): return None
    with open(filepath, 'r') as f:
        for line in f:
            if "finalstate_nocut" in line:
                return int(line.split("(")[0].strip().split()[-1])
    return None

def GetThisCutTable(path):
    raw_table = {
        "TotalEvent": {},
        "21To81": {},
        "81To101": {}
    }

    masses = list(range(110, 245, 5))
    bg_list = ["vbs", "zh_ll_tautau", "zh_ll_ww", "zww", "zz_ll_tautau"]
    sig_list = ["HZ", "VBF"]

    for category in raw_table:
        for col in sig_list + bg_list:
            raw_table[category][col] = {}

    for bg in bg_list:
        f21 = os.path.join(path, f"ISR_{bg}_21To81Hist.root")
        count21, yields21 = extract_yields_and_count(f21, masses, 10) if os.path.exists(f21) else (0.0, {m: 0.0 for m in masses})
        
        f81 = os.path.join(path, f"ISR_{bg}_81To101Hist.root")
        _, yields81 = extract_yields_and_count(f81, masses, 10) if os.path.exists(f81) else (0.0, {m: 0.0 for m in masses})

        for m in masses:
            raw_table["TotalEvent"][bg][f"Mass_{m}"] = count21
            raw_table["21To81"][bg][f"Mass_{m}"] = yields21[m]
            raw_table["81To101"][bg][f"Mass_{m}"] = yields81[m]

    for m in masses:
        f_hz21 = os.path.join(path, f"HLFV_{m}GeV_21To81Hist.root")
        c_hz21, y_hz21 = extract_yields_and_count(f_hz21, [m], 10) if os.path.exists(f_hz21) else (0, {m: 0})
        f_hz81 = os.path.join(path, f"HLFV_{m}GeV_81To101Hist.root")
        _, y_hz81 = extract_yields_and_count(f_hz81, [m], 10) if os.path.exists(f_hz81) else (0, {m: 0})
        
        raw_table["TotalEvent"]["HZ"][f"Mass_{m}"] = c_hz21
        raw_table["21To81"]["HZ"][f"Mass_{m}"] = y_hz21[m]
        raw_table["81To101"]["HZ"][f"Mass_{m}"] = y_hz81[m]

        f_vbf21 = os.path.join(path, f"VBF_HLFV_{m}GeV_21To81Hist.root")
        c_vbf21, y_vbf21 = extract_yields_and_count(f_vbf21, [m], 10) if os.path.exists(f_vbf21) else (0, {m: 0})
        f_vbf81 = os.path.join(path, f"VBF_HLFV_{m}GeV_81To101Hist.root")
        _, y_vbf81 = extract_yields_and_count(f_vbf81, [m], 10) if os.path.exists(f_vbf81) else (0, {m: 0})
        
        raw_table["TotalEvent"]["VBF"][f"Mass_{m}"] = c_vbf21
        raw_table["21To81"]["VBF"][f"Mass_{m}"] = y_vbf21[m]
        raw_table["81To101"]["VBF"][f"Mass_{m}"] = y_vbf81[m]

    return raw_table


def process_worker(process, selection, path):
    return process, selection, GetThisCutTable(path)


if __name__ == "__main__":
    
    # --- NEW CACHE LOGIC ---
    CACHE_FILE = "all_final_acceptance.json"
    
    if os.path.exists(CACHE_FILE):
        print(f"Loading cached data from {CACHE_FILE}...")
        with open(CACHE_FILE, 'r') as f:
            All_Final_Acceptance = json.load(f)
            
    else:
        print("No cache found. Starting data extraction...")
        basedir = os.path.join(os.getcwd(), "Results", "Selections")
        if not os.path.exists(basedir):
            print(f"Directory {basedir} does not exist.")
            sys.exit(1)

        process_types = os.listdir(basedir)
        All_Final_Acceptance = {p: {} for p in process_types}
        
        tasks = []
        
        # Changed max_workers to 4 to prevent your out-of-memory error
        with ProcessPoolExecutor(max_workers=4) as executor:
            for process in process_types:
                process_path = os.path.join(basedir, process)
                if not os.path.isdir(process_path): continue
                
                selection_types = os.listdir(process_path)
                for selection in selection_types:
                    sel_path = os.path.join(process_path, selection, "HIST")
                    if not os.path.isdir(sel_path): continue
                    
                    tasks.append(executor.submit(process_worker, process, selection, sel_path))
                    
            for future in as_completed(tasks):
                proc, sel, table = future.result()
                All_Final_Acceptance[proc][sel] = table

        print("Data extraction complete.")
        
        # Save the dictionary to a JSON file
        print(f"Saving extracted data to {CACHE_FILE}...")
        with open(CACHE_FILE, 'w') as f:
            json.dump(All_Final_Acceptance, f, indent=4)
            
    print("Data is ready to use!")

    # df_dict = { ... }