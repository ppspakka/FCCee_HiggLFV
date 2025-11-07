#pragma once

#include <functional>
#include <map>
#include <string>
#include <vector>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include "types.h"

namespace hlfv {

class HistogramManager {
public:
    using VarFunc = std::function<double(const Event&, const Meta&)>;

    // 1D variable specification
    struct VarSpec {
        std::string name;     // histogram base name
        std::string title;    // ROOT histogram title (can include axis labels "Title;X;Y")
        int nbins; double xmin; double xmax;
        VarFunc compute;      // returns value or NaN if unavailable
    };

    // 2D variable specification
    struct Var2DSpec {
        std::string name;      // histogram base name
        std::string title;     // ROOT histogram title ("Title;X;Y")
        int nbinsX; double xmin; double xmax;
        int nbinsY; double ymin; double ymax;
        VarFunc computeX;      // X value or NaN
        VarFunc computeY;      // Y value or NaN
    };

    HistogramManager(const std::vector<std::string>& stepNames,
                     const std::vector<VarSpec>& variables1D,
                     const std::vector<Var2DSpec>& variables2D = {});

    void fill(size_t stepIndex, const Event& evt, const Meta& meta, double weight=1.0);
    void writeAll(TFile* fout);

private:
    std::vector<std::string> stepNames_;
    std::vector<VarSpec> variables1D_;
    std::vector<Var2DSpec> variables2D_;
    std::map<std::string, TH1F*> hists1D_;
    std::map<std::string, TH2F*> hists2D_;
};

} // namespace hlfv

