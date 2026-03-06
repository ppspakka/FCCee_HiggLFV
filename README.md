# FCCee_HiggLFV

### How to run

```bash
root -l -b -q 'analyze_pipeline.cpp("PATH_TO_ROOT_FILE_OR_PATTERN","OUTPUT_ROOT_FILE","CONFIG_JSON_FILE")'
```
- `PATH_TO_ROOT_FILE_OR_PATTERN`: e.g. `samples/HMuTauE_LFV_125.root` or `samples/` (to process all ROOT files in the directory)
- `OUTPUT_ROOT_FILE`: e.g. `out_hist_cfg.root` contains the histograms after applying the pipeline selections
- `CONFIG_JSON_FILE`: e.g. `pipeline.json` defines the selection pipeline and parameters

### How to add new selections

1) Implement the new selection class in `selections.cpp` (e.g. `HToMuTauESelection`), following the structure of existing selections.
2) Add the new selection class declaration in `selection.h`.
3) In `pipeline_config.cpp`, add a case in the selection factory function to create an instance of your new selection class when the corresponding name is encountered in the config.
4) In `pipeline.json`, add a new entry in the "selections" list with the name that matches the case you added in `pipeline_config.cpp` and specify any parameters needed for that selection.

### How to add new parameters

1) Define the new parameter in `type.h` (e.g. `double my_new_param;`).
2) In `pipeline_config.cpp`, update the code that reads the config file to parse the new parameter and store it in the `Parameters` struct.
3) In `pipeline.json`, add the new parameter under the "parameters" section with an appropriate value.
4) (Optional) Add print statements in `analyze_pipeline.cpp` to verify that the new parameter is being read correctly during execution.