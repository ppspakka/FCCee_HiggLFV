# FCCee_HiggLFV

## Modular analysis pipeline (ROOT macro)

A modular analysis pipeline is provided in `analyze_pipeline.cpp`. It builds a sequence of selections and produces per-cut histograms with names like:

- 00_Initial_<var>
- 01_Z_to_ll_<var>
- 02_H_to_mue_<var>
- 03_MET_dphi_<var>

Current variables booked per cut:

- n_muons
- z_mass
- z_mass_diff
- dphi_e_met

### How to run

1) Activate your ROOT-enabled conda environment (as used in this project):

```bash
conda activate HEP_env
```

2) Run the macro non-interactively, pointing to your input file/pattern and desired output file. You can optionally pass a config file (defaults to `pipeline.json`):

```bash
# Using defaults baked in (or pipeline.json if present)
root -l -b -q 'analyze_pipeline.cpp("samples/HMuTauE_LFV_125.root","out_hist_modular.root")'

# Explicit config path
root -l -b -q 'analyze_pipeline.cpp("samples/HMuTauE_LFV_125.root","out_hist_cfg.root","pipeline.json")'
```

The macro will print a small cutflow summary and write all histograms plus a `cutflow` TH1F to the output ROOT file.

### Code layout

- analyze_pipeline.cpp: Thin ROOT macro entry that wires the pipeline and includes the implementation units to avoid a separate build; loads selections and parameters from JSON.
- include/
  - types.h: Parameters, Event, Meta structs
  - iselection.h: Selection interface
  - selections.h: Selection declarations (AtLeastOneMuon, Z_to_ll, H_to_mue, MET_dphi)
  - selection_factory.h: Factory to build selections by name
  - config_parser.h: Minimal JSON parser for pipeline config
  - histman.h: Histogram manager interface
  - utils.h: Small helpers (e.g., file listing)
- src/
  - selections.cpp: Selection implementations
  - selection_factory.cpp: Name → selection mapping
  - config_parser.cpp: Parses `pipeline.json`
  - histman.cpp: HistogramManager implementation
  - utils.cpp: getFileList implementation

- pipeline.json: Example config listing selection order and run parameters

Note: The macro includes `Delphes.C` once so all sources share the same Delphes class definition at compile time inside the ROOT interpreter.

### Adjusting parameters

Thresholds such as the Z-mass window and the MET dphi requirement live in the `Parameters` struct (see `include/types.h`) and can be set in `pipeline.json` under `parameters`.

### Outputs

The output ROOT file contains per-cut histograms named `NN_StepName_Var` and a `cutflow` histogram whose bin labels match the step names (Initial and each enabled selection in order).

## Add a new selection (cut)

You can add a new cut without touching the event loop by implementing a selection module and wiring it in the factory and config.

1) Define the selection class
   - File: `include/selections.h`
     - Declare a struct deriving from `ISelection` with `name() const` and `apply(const Event&, Meta&, const Parameters&)`.
   - File: `src/selections.cpp`
     - Implement `name()` (used for step naming) and `apply()` (return true/false, and update `Meta` as needed).

2) Add any state to pass between cuts
   - File: `include/types.h`
     - Extend `Meta` with fields your selection needs to read/write (e.g., chosen indices, computed angles/masses). Use NaN/-1 as “unset” for safe plotting.

3) Ensure needed branches are enabled
   - File: `analyze_pipeline.cpp`
     - In the binding section, add `chain->SetBranchStatus("BranchPrefix*", 1);` for any additional collections your selection reads.

4) Register the selection name in the factory
   - File: `src/selection_factory.cpp`
     - Map a user-facing name (e.g., `"MyNewCut"`) to `std::make_unique<MyNewCut>()`.
     - You can support a few aliases if useful.

5) Add the selection to the pipeline config
   - File: `pipeline.json`
     - Insert your selection in `selections` at the desired position: `{ "name": "MyNewCut", "enabled": true }`.
     - Reorder or toggle cuts by editing this list. No code changes required.

6) Optional: add new variables (plots)
   - File: `analyze_pipeline.cpp`
     - In the `variables` list, add a `VarSpec` with a compute lambda. Return NaN when the value isn’t valid yet; the histogram filler will skip it.

7) Optional: add new parameters
   - File: `include/types.h` (add to `Parameters`)
   - File: `pipeline.json` (set under `parameters`)
   - Use them inside your selection’s `apply()`.

Notes:
- Dependencies between cuts: implement prerequisite checks inside `apply()` (e.g., if your cut needs an electron index from a previous step, return false when it’s missing). We can add an optional startup validator if you want hard errors on mis-ordered configs.
- Keep `analyze_eeZH.cpp` untouched; it’s a reference implementation only.
