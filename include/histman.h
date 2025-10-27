#pragma once

#include <functional>
#include <map>
#include <string>
#include <vector>
#include <TH1F.h>
#include <TFile.h>
#include "types.h"

namespace hlfv {

class HistogramManager {
public:
    using VarFunc = std::function<double(const Event&, const Meta&)>;

    struct VarSpec {
        std::string name;
        std::string title;
        int nbins; double xmin; double xmax;
        VarFunc compute;
    };

    HistogramManager(const std::vector<std::string>& stepNames,
                     const std::vector<VarSpec>& variables);

    void fill(size_t stepIndex, const Event& evt, const Meta& meta, double weight=1.0);
    void writeAll(TFile* fout);

private:
    std::vector<std::string> stepNames_;
    std::vector<VarSpec> variables_;
    std::map<std::string, TH1F*> hists_;
};

} // namespace hlfv
