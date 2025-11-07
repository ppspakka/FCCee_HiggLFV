#include "../include/histman.h"

#include <cmath>

namespace hlfv {

HistogramManager::HistogramManager(const std::vector<std::string>& stepNames,
                                   const std::vector<VarSpec>& variables1D,
                                   const std::vector<Var2DSpec>& variables2D)
: stepNames_(stepNames), variables1D_(variables1D), variables2D_(variables2D)
{
    // Create 1D histograms
    for (size_t i = 0; i < stepNames_.size(); ++i) {
        for (const auto& v : variables1D_) {
            char hname[256];
            snprintf(hname, sizeof(hname), "%s_%s", stepNames_[i].c_str(), v.name.c_str());
            TH1F* h = new TH1F(hname, (v.title.size()? v.title.c_str(): v.name.c_str()), v.nbins, v.xmin, v.xmax);
            h->Sumw2();
            hists1D_[hname] = h;
        }
        for (const auto& v2 : variables2D_) {
            char hname2[256];
            snprintf(hname2, sizeof(hname2), "%s_%s", stepNames_[i].c_str(), v2.name.c_str());
            TH2F* h2 = new TH2F(hname2, (v2.title.size()? v2.title.c_str(): v2.name.c_str()),
                                v2.nbinsX, v2.xmin, v2.xmax,
                                v2.nbinsY, v2.ymin, v2.ymax);
            h2->Sumw2();
            hists2D_[hname2] = h2;
        }
    }
}

void HistogramManager::fill(size_t stepIndex, const Event& evt, const Meta& meta, double weight) {
    if (stepIndex >= stepNames_.size()) return;
    const auto& step = stepNames_[stepIndex];
    // 1D
    for (const auto& v : variables1D_) {
        double val = v.compute(evt, meta);
        if (std::isnan(val)) continue; // skip unavailable values
        std::string hname = step + std::string("_") + v.name;
        auto it = hists1D_.find(hname);
        if (it != hists1D_.end()) it->second->Fill(val, weight);
    }
    // 2D
    for (const auto& v2 : variables2D_) {
        double x = v2.computeX(evt, meta);
        double y = v2.computeY(evt, meta);
        if (std::isnan(x) || std::isnan(y)) continue; // skip if either axis missing
        std::string hname2 = step + std::string("_") + v2.name;
        auto it2 = hists2D_.find(hname2);
        if (it2 != hists2D_.end()) it2->second->Fill(x, y, weight);
    }
}

void HistogramManager::writeAll(TFile* fout) {
    if (!fout) return;
    fout->cd();
    for (auto& kv : hists1D_) kv.second->Write();
    for (auto& kv : hists2D_) kv.second->Write();
}

} // namespace hlfv

