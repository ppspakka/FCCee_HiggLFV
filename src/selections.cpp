#include "../include/selections.h"

#include <cmath>
#include <TLorentzVector.h>

namespace hlfv {

constexpr double Me = 0.000511;
constexpr double Mmu = 0.10565837;

std::string AtLeastOneMuonSelection::name() const { return "AtLeastOneMuon"; }

bool AtLeastOneMuonSelection::apply(const Event& evt, Meta&, const Parameters&) {
    return evt.d && (evt.d->Muon_size > 0);
}

std::string ZToLLSelection::name() const { return "Z_to_ll"; }

bool ZToLLSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    double bestDiff = 1e9;
    int best_i = -1, best_j = -1, best_flav = -1; // 0=e,1=mu
    double bestMass = std::numeric_limits<double>::quiet_NaN();

    // electrons
    for (int i = 0; i < evt.d->Electron_size; ++i) {
        for (int j = i+1; j < evt.d->Electron_size; ++j) {
            if (evt.d->Electron_Charge[i]*evt.d->Electron_Charge[j] >= 0) continue; // OS
            TLorentzVector l1, l2;
            l1.SetPtEtaPhiM(evt.d->Electron_PT[i], evt.d->Electron_Eta[i], evt.d->Electron_Phi[i], Me);
            l2.SetPtEtaPhiM(evt.d->Electron_PT[j], evt.d->Electron_Eta[j], evt.d->Electron_Phi[j], Me);
            double mass = (l1+l2).M();
            double diff = std::abs(mass - cfg.z_mass);
            if (diff < cfg.z_mass_window && diff < bestDiff) {
                bestDiff = diff; best_i = i; best_j = j; best_flav = 0; bestMass = mass;
            }
        }
    }

    // muons
    for (int i = 0; i < evt.d->Muon_size; ++i) {
        for (int j = i+1; j < evt.d->Muon_size; ++j) {
            if (evt.d->Muon_Charge[i]*evt.d->Muon_Charge[j] >= 0) continue; // OS
            TLorentzVector l1, l2;
            l1.SetPtEtaPhiM(evt.d->Muon_PT[i], evt.d->Muon_Eta[i], evt.d->Muon_Phi[i], Mmu);
            l2.SetPtEtaPhiM(evt.d->Muon_PT[j], evt.d->Muon_Eta[j], evt.d->Muon_Phi[j], Mmu);
            double mass = (l1+l2).M();
            double diff = std::abs(mass - cfg.z_mass);
            if (diff < cfg.z_mass_window && diff < bestDiff) {
                bestDiff = diff; best_i = i; best_j = j; best_flav = 1; bestMass = mass;
            }
        }
    }

    if (best_flav == -1) return false;
    meta.z_l1 = best_i; meta.z_l2 = best_j; meta.z_flavor = best_flav;
    meta.z_mass = bestMass; meta.z_mass_diff = bestDiff;
    return true;
}

std::string HToMuESelection::name() const { return "H_to_mue"; }

bool HToMuESelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    int zl1 = meta.z_l1, zl2 = meta.z_l2;

    // Handle for invalid indices
    if (zl1 < 0 || zl2 < 0) return false;

    int idx_mu = -1;
    for (int m = 0; m < evt.d->Muon_size; ++m) {
        if (meta.z_flavor == 1 && (m == zl1 || m == zl2)) continue; // exclude Z muons
        if (evt.d->Muon_PT[m] > cfg.mu_pt_min) {
            if (idx_mu == -1) idx_mu = m; else return false; // more than one
        }
    }
    if (idx_mu < 0) return false;

    int idx_e = -1;
    for (int e = 0; e < evt.d->Electron_size; ++e) {
        if (meta.z_flavor == 0 && (e == zl1 || e == zl2)) continue; // exclude Z electrons
        if (evt.d->Electron_PT[e] > cfg.e_pt_min) {
            if (idx_e == -1) idx_e = e; else return false; // more than one
        }
    }
    if (idx_e < 0) return false;

    if (evt.d->Muon_Charge[idx_mu] * evt.d->Electron_Charge[idx_e] >= 0) return false; // OS

    meta.h_mu = idx_mu;
    meta.h_e = idx_e;

    if (evt.d->MissingET_size > 0) {
        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(1.0, 0.0, evt.d->Electron_Phi[idx_e], 0.0);
        v2.SetPtEtaPhiM(1.0, 0.0, evt.d->MissingET_Phi[0], 0.0);
        meta.dphi_e_met = std::abs(v1.DeltaPhi(v2));
    }
    return true;
}

std::string METDphiSelection::name() const { return "MET_dphi"; }

bool METDphiSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    if (meta.h_e < 0) return false; // need H electron selected
    if (evt.d->MissingET_size <= 0) return false;
    if (std::isnan(meta.dphi_e_met)) {
        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(1.0, 0.0, evt.d->Electron_Phi[meta.h_e], 0.0);
        v2.SetPtEtaPhiM(1.0, 0.0, evt.d->MissingET_Phi[0], 0.0);
        meta.dphi_e_met = std::abs(v1.DeltaPhi(v2));
    }
    return meta.dphi_e_met < cfg.max_dphi_e_met;
}

} // namespace hlfv
