#! /usr/bin/env python3

import os

# Define the path to the directory containing the cards
merged_card_dir = "Merged_ZH_VBF_110To220_V1"
process_list = os.listdir(merged_card_dir)
for process in process_list:
    this_path = os.path.join(merged_card_dir, process)
    card_list = os.listdir(this_path)
    for card in card_list:
        mass = card.split("_")[-1].replace(".txt", "")
        try:
            int_mass = int(mass)
        except ValueError:
            print(f"Warning: Could not convert mass '{mass}' to an integer. Skipping card '{card}'.")
            continue
        os.system(f"sbatch Run_Upperlimit.sh -m {int_mass} -c {card} -p {this_path}")
        print(f"Submitted job for mass {int_mass} GeV using card {card} in process {process} with path {this_path}.")