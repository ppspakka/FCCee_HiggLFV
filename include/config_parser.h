#pragma once

#include <string>
#include <vector>
#include "types.h"

struct SelectionSpec {
    std::string name;
    bool enabled = true;
};

struct PipelineConfig {
    Parameters params; // thresholds
    std::vector<SelectionSpec> selections; // in order
};

// Minimal JSON reader that expects keys: parameters{}, selections[] with objects {name, enabled}
// Returns true on success; false if file missing or parse error (in which case caller can fallback to defaults)
bool loadPipelineConfig(const std::string& filepath, PipelineConfig& out);
