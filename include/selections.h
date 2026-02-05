#pragma once

#include "iselection.h"

namespace hlfv {

// Concrete selections

struct EmptySelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct FinalState_NoCut : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};
struct LeptonSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct ZToLLSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct ZCandidateSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct HToMuTauESelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct HToETauMuSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct METEDphiSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct METMuDphiSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct RecoilMassSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct HCandidateSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

} // namespace hlfv
