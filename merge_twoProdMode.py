import os
import glob
import math

# --- Configuration ---
DATACARD_PATH = {
    'ZH': './ZH_110To220_integerBin',
    'VBF': './VBF_110To220_integerBin'
}
OUTPUT_BASE = './Merged_ZH_VBF_110To220'

# --- Core Functions ---
def align_columns(label, values, widths):
    row = f"{label:<20}"
    for val, width in zip(values, widths):
        row += f"{str(val):<{width + 3}}"
    return row.rstrip()

def parse_datacard(filepath):
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f.readlines()]

    data = {
        "header_comments": [], "imax": 0, "jmax": 0, "kmax": 0,
        "bins": [], "observation": [],
        "col_bins": [], "col_processes": [], "col_ids": [], "col_rates": [],
        "systematics": [], "parameters": []
    }

    for line in lines:
        if not line.strip(): continue
        parts = line.split()
        
        if line.startswith("Combination") or line.startswith("#"): data["header_comments"].append(line)
        elif parts[0] == "imax": data["imax"] = parts[1]
        elif parts[0] == "jmax": data["jmax"] = parts[1]
        elif parts[0] == "kmax": data["kmax"] = parts[1]
        elif parts[0] == "bin" and len(data["observation"]) == 0: data["bins"] = parts[1:]
        elif parts[0] == "observation": data["observation"] = parts[1:]
        elif parts[0] == "bin" and len(data["observation"]) > 0: data["col_bins"] = parts[1:]
        elif parts[0] == "process":
            if any(p.isalpha() for p in parts[1:]): data["col_processes"] = parts[1:]
            else: data["col_ids"] = parts[1:]
        elif parts[0] == "rate": data["col_rates"] = parts[1:]
        elif len(parts) > 1 and parts[1] in ["lnN", "shape", "gmN"]:
            data["systematics"].append({"name": parts[0], "type": parts[1], "values": parts[2:]})
        elif not parts[0].startswith("-"):
            data["parameters"].append(line)
    return data

def construct_datacard(data, output_path):
    all_table_rows = [data['col_bins'], data['col_processes'], data['col_ids'], data['col_rates']]
    for s in data['systematics']: all_table_rows.append(s['values'])
    
    col_widths = []
    for i in range(len(data['col_bins'])):
        widths_in_col = [len(str(row[i])) for row in all_table_rows if i < len(row)]
        col_widths.append(max(widths_in_col) if widths_in_col else 8)

    with open(output_path, 'w') as f:
        for comment in data["header_comments"]: f.write(f"{comment}\n")
        
        # Replaced parsed values with '*' for auto-counting
        f.write("imax * number of bins\n")
        f.write("jmax * number of processes minus 1\n")
        f.write("kmax * number of nuisance parameters\n")
        f.write("-" * 130 + "\n")
        
        f.write(f"{'bin':<20} {'   '.join(data['bins'])}\n")
        f.write(f"{'observation':<20} {'   '.join(data['observation'])}\n")
        f.write("-" * 130 + "\n")
        
        f.write(align_columns("bin", data['col_bins'], col_widths) + "\n")
        f.write(align_columns("process", data['col_processes'], col_widths) + "\n")
        f.write(align_columns("process", data['col_ids'], col_widths) + "\n")
        f.write(align_columns("rate", data['col_rates'], col_widths) + "\n")
        f.write("-" * 130 + "\n")
        
        for s in data["systematics"]:
            f.write(align_columns(f"{s['name']:<14} {s['type']}", s['values'], col_widths) + "\n")
        f.write("\n")
        
        for p in data["parameters"]: f.write(f"{p}\n")

# --- Processing Functions ---
def rename_template_signal(data, prefix):
    """Renames the ID 0 signal in the template to include the prefix, updating tail params."""
    old_sig_name = ""
    for i, pid in enumerate(data['col_ids']):
        if pid == '0':
            old_sig_name = data['col_processes'][i]
            data['col_processes'][i] = f"{prefix}_{old_sig_name}"
    
    new_params = []
    for param in data['parameters']:
        if len(param.split()) >= 4 and param.split()[3] == old_sig_name and not param.startswith("nuisance"):
            new_params.append(param.replace(old_sig_name, f"{prefix}_{old_sig_name}"))
        else:
            new_params.append(param)
    data['parameters'] = new_params
    return data, old_sig_name

def inject_signal(template_data, new_rates, prefix, base_sig_name, sig_id):
    """Injects a new signal into the template data directly after the last signal."""
    # Find insertion points (after existing signals, before backgrounds)
    insertion_indices = []
    for bin_name in template_data['bins']:
        # Find the last process in this bin that is a signal (ID <= 0)
        last_sig_idx = -1
        for i, (b, pid) in enumerate(zip(template_data['col_bins'], template_data['col_ids'])):
            if b == bin_name and int(pid) <= 0:
                last_sig_idx = i
        if last_sig_idx != -1:
            insertion_indices.append(last_sig_idx + 1)
    
    # Insert backwards to avoid index shifting issues
    new_sig_name = f"{prefix}_{base_sig_name}"
    for idx, bin_name, rate in zip(reversed(insertion_indices), reversed(template_data['bins']), reversed(new_rates)):
        template_data['col_bins'].insert(idx, bin_name)
        template_data['col_processes'].insert(idx, new_sig_name)
        template_data['col_ids'].insert(idx, str(sig_id))
        template_data['col_rates'].insert(idx, str(rate))
        
        # Mirror systematics from the base signal (assuming ID 0 is always at idx - abs(sig_id))
        ref_idx = idx - 1 
        for s in template_data['systematics']:
            s['values'].insert(idx, s['values'][ref_idx])

    # Mirror rate params
    new_params = list(template_data['parameters'])
    for param in template_data['parameters']:
        parts = param.split()
        if len(parts) >= 4 and parts[3].endswith(base_sig_name) and not param.startswith("nuisance"):
            new_params.append(param.replace(parts[3], new_sig_name))
    
    template_data['parameters'] = sorted(list(set(new_params))) # Clean duplicates just in case
    return template_data

def get_signal_rates(data):
    return [float(rate) for pid, rate in zip(data['col_ids'], data['col_rates']) if pid == '0']

def get_bkg_rates(data):
    return [float(rate) for pid, rate in zip(data['col_ids'], data['col_rates']) if int(pid) > 0]

# --- Main Execution ---
def main():
    if not os.path.exists(OUTPUT_BASE): os.makedirs(OUTPUT_BASE)

    signals = list(DATACARD_PATH.keys())
    template_key = signals[0]
    template_base = DATACARD_PATH[template_key]
    
    combined_dirs = [d for d in os.listdir(template_base) if "combined" in d]

    for c_dir in combined_dirs:
        out_dir_path = os.path.join(OUTPUT_BASE, c_dir)
        if not os.path.exists(out_dir_path): os.makedirs(out_dir_path)

        datacards = glob.glob(os.path.join(template_base, c_dir, "*.txt"))
        
        for template_card_path in datacards:
            filename = os.path.basename(template_card_path)
            
            # 1. Load Template
            template_data = parse_datacard(template_card_path)
            template_bkgs = get_bkg_rates(template_data)
            
            # 2. Rename Template Signal (e.g., HMuTauE -> ZH_HMuTauE)
            template_data, base_sig_name = rename_template_signal(template_data, template_key)

            # 3. Iterate over remaining signals dynamically
            for i, sig_key in enumerate(signals[1:], start=1):
                sig_card_path = os.path.join(DATACARD_PATH[sig_key], c_dir, filename)
                if not os.path.exists(sig_card_path):
                    print(f"Skipping {sig_key} for {filename} (not found)")
                    continue
                
                sig_data = parse_datacard(sig_card_path)
                
                # Sanity check backgrounds
                sig_bkgs = get_bkg_rates(sig_data)
                for tb, sb in zip(template_bkgs, sig_bkgs):
                    if not math.isclose(tb, sb, rel_tol=1e-4):
                        raise ValueError(f"Bkg mismatch in {filename} between {template_key} and {sig_key}")

                # Extract and inject
                sig_rates = get_signal_rates(sig_data)
                sig_id = -i # Assigns -1, -2, etc.
                template_data = inject_signal(template_data, sig_rates, sig_key, base_sig_name, sig_id)

            # 4. Save
            out_file = os.path.join(out_dir_path, filename)
            construct_datacard(template_data, out_file)
            print(f"Processed: {filename}")

if __name__ == "__main__":
    main()
