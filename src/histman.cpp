#include "../include/histman.h"

#include <cmath>

namespace hlfv {

HistogramManager::HistogramManager(const std::vector<std::string>& stepNames,
                                   const std::vector<VarSpec>& variables)
: stepNames_(stepNames), variables_(variables)
{
    for (size_t i = 0; i < stepNames_.size(); ++i) {
        for (const auto& v : variables_) {
            char hname[256];
            snprintf(hname, sizeof(hname), "%s_%s", stepNames_[i].c_str(), v.name.c_str());
            TH1F* h = new TH1F(hname, (v.title.size()? v.title.c_str(): v.name.c_str()), v.nbins, v.xmin, v.xmax);
            h->Sumw2();
            hists_[hname] = h;
        }
    }
}

void HistogramManager::fill(size_t stepIndex, const Event& evt, const Meta& meta, double weight) {
    if (stepIndex >= stepNames_.size()) return;
    const auto& step = stepNames_[stepIndex];
    for (const auto& v : variables_) {
        double val = v.compute(evt, meta);
        if (std::isnan(val)) continue; // skip unavailable values
        std::string hname = step + std::string("_") + v.name;
        auto it = hists_.find(hname);
        if (it != hists_.end()) it->second->Fill(val, weight);
    }
}

void HistogramManager::writeAll(TFile* fout) {
    if (!fout) return;
    fout->cd();
    for (auto& kv : hists_) kv.second->Write();
}

} // namespace hlfv
