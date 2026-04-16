#!/usr/bin/env python3
import os
import glob
import shutil
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Merge LFV Higgs datacards and submit limits.")
    parser.add_argument("--submit", action="store_true", help="Submit sbatch jobs after merging.")
    args = parser.parse_args()

    channels = ['etaumu', 'mutaue', 'mue'] 
    
    # Create an isolated environment for Combine scripts
    combine_env = os.environ.copy()
    combine_env["PYTHONNOUSERSITE"] = "1"
    
    for ch in channels:
        dir_21 = f"datacards_{ch}_21To81"
        dir_81 = f"datacards_{ch}_81To101"
        dir_comb = f"datacards_{ch}_combined"

        if not os.path.exists(dir_21) or not os.path.exists(dir_81):
            print(f"Skipping {ch}: Missing input directories.")
            continue

        os.makedirs(dir_comb, exist_ok=True)

        cards_21 = glob.glob(os.path.join(dir_21, "datacard_*.txt"))
        
        for card_21 in cards_21:
            filename = os.path.basename(card_21)
            card_81 = os.path.join(dir_81, filename)
            card_comb = os.path.join(dir_comb, filename)

            if not os.path.exists(card_81):
                print(f"Warning: {card_81} missing. Skipping {filename}.")
                continue

            # Execute combineCards.py with isolated environment and native Python file writing
            cmd = ["combineCards.py", f"bin21_81={card_21}", f"bin81_101={card_81}"]
            with open(card_comb, "w") as f_out:
                subprocess.run(cmd, env=combine_env, stdout=f_out, check=True)

        print(f"Merged datacards created in {dir_comb}/")

        # Copy auxiliary scripts
        for script in ["run_limits.py", "slurm_submit.slurm"]:
            src = os.path.join(dir_21, script)
            if os.path.exists(src):
                shutil.copy(src, dir_comb)

        # Execute submission within the combined directory
        if args.submit:
            print(f"Submitting job for {ch} combined...")
            cwd = os.getcwd()
            os.chdir(dir_comb)
            # Submit using the standard environment, or combine_env if your Slurm script also invokes combine directly
            subprocess.run(["sbatch", "slurm_submit.slurm"], check=True)
            os.chdir(cwd)

if __name__ == "__main__":
    main()