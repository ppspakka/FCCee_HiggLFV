#pragma once

#include "iselection.h"

namespace hlfv {

// Concrete selections
struct AtLeastOneMuonSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct ZToLLSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct HToMuESelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct METDphiSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

} // namespace hlfv
