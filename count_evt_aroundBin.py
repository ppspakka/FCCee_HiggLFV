#!/usr/bin/env python3
import sys
import re
import ROOT

def get_yield_around_mass(hist, mass, window):
    """Sums bin contents in [mass - window, mass + window]."""
    if not hist:
        return 0.0
    
    yield_sum = 0.0
    axis = hist.GetXaxis()
    nbins = axis.GetNbins()
    
    low_edge = mass - window
    high_edge = mass + window

    for ibin in range(1, nbins + 1):
        bin_center = axis.GetBinCenter(ibin)
        if low_edge <= bin_center <= high_edge:
            yield_sum += hist.GetBinContent(ibin)
            
    return yield_sum

def find_histograms(tfile, patterns):
    """Finds histograms in a ROOT file matching a list of regex patterns."""
    found_hists = {p: None for p in patterns}
    
    def walk(directory):
        for key in directory.GetListOfKeys():
            obj = key.ReadObj()
            if obj.InheritsFrom("TDirectory"):
                walk(obj)
            elif obj.InheritsFrom("TH1"):
                name = obj.GetName()
                for p in patterns:
                    if re.search(p, name):
                        found_hists[p] = obj.Clone()
                        found_hists[p].SetDirectory(0)
    
    walk(tfile)
    return found_hists

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 count_events.py <input.root> <mass> <window_size>")
        sys.exit(1)

    file_path = sys.argv[1]
    target_mass = float(sys.argv[2])
    window = float(sys.argv[3])

    # 1. Added the initial count histogram to the patterns
    init_pattern = r"00_Initial_n_muons$"
    patterns = [
        init_pattern,
        r".*_finalstate_nocut_m_collinear$",
        r".*_finalstate_nocut_m_h_invariant$",
        r".*_finalstate_nocut_m_h_invariant_count$"
    ]

    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Error: Could not open {file_path}")
        sys.exit(1)

    hists = find_histograms(f, patterns)
    f.Close()

    # 2. Extract initial count
    initial_hist = hists.get(init_pattern)
    # Using Integral() is generally safer than GetEntries() for weighted events
    initial_count = initial_hist.Integral() if initial_hist else 0.0

    print(f"\n--- Analysis Results (Window: {target_mass} +/- {window}) ---")
    if initial_count > 0:
        print(f"{'Initial Events':15s} (Hist: {initial_hist.GetName()}): {initial_count:.2f}")
    else:
        print(f"{'Initial Events':15s}: NOT FOUND or EMPTY!")
    
    print("-" * 65)
    
    # 3. Print loop with efficiency calculation
    for pattern, hist in hists.items():
        if pattern == init_pattern: 
            continue  # Skip printing the initial count again in the loop
            
        label = "Collinear Mass" if "collinear" in pattern else "Invariant Mass"
        
        if hist:
            yield_val = get_yield_around_mass(hist, target_mass, window)
            
            # Calculate efficiency %
            efficiency = (yield_val / initial_count * 100) if initial_count > 0 else 0.0
            
            print(f"{label:15s} (Hist: {hist.GetName():25s}): {yield_val:10.4f} events | Eff: {efficiency:7.4f}%")
        else:
            print(f"{label:15s}: Histogram not found!")

if __name__ == "__main__":
    main()

