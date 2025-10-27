#pragma once

#include <memory>
#include <string>
#include <vector>
#include "types.h"
#include "iselection.h"

namespace hlfv {

struct SelectionSpec {
    std::string name;
    bool enabled = true;
};

struct PipelineConfig {
    Parameters params;
    std::vector<SelectionSpec> selections; // ordered
};

// Parse JSON-like config file with keys: parameters{...}, selections[{ name, enabled }]
bool loadPipelineConfig(const std::string& filepath, PipelineConfig& out);

// Instantiate selection object from a name (e.g., "Z_to_ll", "H_to_mue", "MET_dphi")
std::unique_ptr<ISelection> makeSelectionByName(const std::string& name);

// Variables registry
struct Variables {
    static std::vector<class HistogramManager::VarSpec> getDefault();
};

} // namespace hlfv
