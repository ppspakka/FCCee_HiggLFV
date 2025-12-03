#include "../include/selections.h"

#include <cmath>
#include <TLorentzVector.h>

namespace hlfv {

constexpr double Me = 0.000511;
constexpr double Mmu = 0.10565837;

// helper functions extracted from selection logic
namespace {
    // compute absolute delta phi between two phi angles using TLorentzVector::DeltaPhi
    static double deltaPhiFromPhis(double phi1, double phi2) {
        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(1.0, 0.0, phi1, 0.0);
        v2.SetPtEtaPhiM(1.0, 0.0, phi2, 0.0);
        return std::abs(v1.DeltaPhi(v2));
    }

    // compute collinear mass for mu + e using missing ET (returns NaN if no MissingET)
    static double computeCollinearMassMuTauE(const Event& evt, int idx_mu, int idx_e) {
        if (!evt.d || evt.d->MissingET_size <= 0) return std::numeric_limits<double>::quiet_NaN();
        TLorentzVector vmu, ve;
        vmu.SetPtEtaPhiM(evt.d->Muon_PT[idx_mu], evt.d->Muon_Eta[idx_mu], evt.d->Muon_Phi[idx_mu], Mmu);
        ve.SetPtEtaPhiM(evt.d->Electron_PT[idx_e], evt.d->Electron_Eta[idx_e], evt.d->Electron_Phi[idx_e], Me);
        double px_miss = evt.d->MissingET_MET[0] * std::cos(evt.d->MissingET_Phi[0]);
        double py_miss = evt.d->MissingET_MET[0] * std::sin(evt.d->MissingET_Phi[0]);
        double pz_miss = (ve.Pz()/ve.Pt()) * std::sqrt(px_miss*px_miss + py_miss*py_miss);
        double e_miss = std::sqrt(px_miss*px_miss + py_miss*py_miss + pz_miss*pz_miss);
        TLorentzVector vmiss;
        vmiss.SetPxPyPzE(px_miss, py_miss, pz_miss, e_miss);
        return (vmu + ve + vmiss).M();
    }
    static double computeCollinearMassETauMu(const Event& evt, int idx_e, int idx_mu) {
        if (!evt.d || evt.d->MissingET_size <= 0) return std::numeric_limits<double>::quiet_NaN();
        TLorentzVector ve, vmu;
        ve.SetPtEtaPhiM(evt.d->Electron_PT[idx_e], evt.d->Electron_Eta[idx_e], evt.d->Electron_Phi[idx_e], Me);
        vmu.SetPtEtaPhiM(evt.d->Muon_PT[idx_mu], evt.d->Muon_Eta[idx_mu], evt.d->Muon_Phi[idx_mu], Mmu);
        double px_miss = evt.d->MissingET_MET[0] * std::cos(evt.d->MissingET_Phi[0]);
        double py_miss = evt.d->MissingET_MET[0] * std::sin(evt.d->MissingET_Phi[0]);
        double pz_miss = (vmu.Pz()/vmu.Pt()) * std::sqrt(px_miss*px_miss + py_miss*py_miss);
        double e_miss = std::sqrt(px_miss*px_miss + py_miss*py_miss + pz_miss*pz_miss);
        TLorentzVector vmiss;
        vmiss.SetPxPyPzE(px_miss, py_miss, pz_miss, e_miss);
        return (vmu + ve + vmiss).M();
    }

    // compute the transverse mass, mode=0 for e and 1 for mu
    static double computeTransverseMass(const Event& evt, int mode, int idx_lep) {
        if (!evt.d || evt.d->MissingET_size <= 0) return std::numeric_limits<double>::quiet_NaN();
        double pt_lep = 0.0;
        double phi_lep = 0.0;
        if (mode == 0) {
            // electron
            pt_lep = evt.d->Electron_PT[idx_lep];
            phi_lep = evt.d->Electron_Phi[idx_lep];
        } else if (mode == 1) {
            // muon
            pt_lep = evt.d->Muon_PT[idx_lep];
            phi_lep = evt.d->Muon_Phi[idx_lep];
        } else {
            return std::numeric_limits<double>::quiet_NaN();
        }
        double met = evt.d->MissingET_MET[0];
        double phi_met = evt.d->MissingET_Phi[0];
        double delta_phi = deltaPhiFromPhis(phi_lep, phi_met);
        return std::sqrt(2 * pt_lep * met * (1 - std::cos(delta_phi)));
    }
} // anonymous namespace

// test add empty cut
std::string EmptySelection::name() const { return "EmptySelection"; }
bool EmptySelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    return true;
}

// Final state: no cut applied
std::string FinalState_NoCut::name() const { return "FinalState_NoCut"; }
bool FinalState_NoCut::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    return true;
}


/*
    I) Lepton selections:
        1) Filter all leptons (e, mu) with pT > threshold (eg. 10 GeV)
        2) Require either
        - 1 muon + 3 electrons or 1 electron + 3 muons
        3) Conservative of charge requirement: net charge = 0
    Meta info to store:
        - l1flavor (0=e,1=mu) -> To store prompt lepton from H decay
        - l2flavor (0=e,1=mu)
        - l3flavor (0=e,1=mu)
        - l4flavor (0=e,1=mu)
        - l1 index
        - l2 index
        - l3 index
        - l4 index
*/

std::string LeptonSelection::name() const { return "LeptonSelection"; }
bool LeptonSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    // store selected leptons
    struct Lepton {
        int index;
        int flavor; // 0=e,1=mu
        double pt;
        int charge;
    };
    std::vector<Lepton> selected_muons;
    std::vector<Lepton> selected_electrons;

    // select electrons
    for (int i = 0; i < evt.d->Electron_size; ++i) {
        if (evt.d->Electron_PT[i] > cfg.lepton_pt_min) {
            selected_electrons.push_back(Lepton{i, 0, evt.d->Electron_PT[i], evt.d->Electron_Charge[i]});
        }
    }

    // select muons
    for (int i = 0; i < evt.d->Muon_size; ++i) {
        if (evt.d->Muon_PT[i] > cfg.lepton_pt_min) {
            selected_muons.push_back(Lepton{i, 1, evt.d->Muon_PT[i], evt.d->Muon_Charge[i]});
        }
    }

    // check for 1 mu + 3 e or 1 e + 3 mu
    if (!((selected_muons.size() == 1 && selected_electrons.size() == 3) ||
          (selected_muons.size() == 3 && selected_electrons.size() == 1))) {
        return false;
    }
    // check net charge = 0
    int net_charge = 0;
    for (const auto& mu : selected_muons) {
        net_charge += mu.charge;
    }
    for (const auto& ele : selected_electrons) {
        net_charge += ele.charge;
    }
    if (net_charge != 0) return false;

    // store in meta (1st lepton = prompt lepton from H decay, or lepton with only one copy)
    // if comb: 1 mu + 3 e, 1st = muon, else 1st = electron
    
    if (selected_muons.size() == 1) {
        // 1 mu + 3 e
        meta.l1flavor = 1;
        meta.l1_index = selected_muons[0].index;
        // sort electrons by pT descending
        std::sort(selected_electrons.begin(), selected_electrons.end(),
                  [](const Lepton& a, const Lepton& b) { return a.pt > b.pt; });
        meta.l2flavor = 0; meta.l2_index = selected_electrons[0].index;
        meta.l3flavor = 0; meta.l3_index = selected_electrons[1].index;
        meta.l4flavor = 0; meta.l4_index = selected_electrons[2].index;
    } else {
        // 1 e + 3 mu
        meta.l1flavor = 0;
        meta.l1_index = selected_electrons[0].index;
        // sort muons by pT descending
        std::sort(selected_muons.begin(), selected_muons.end(),
                  [](const Lepton& a, const Lepton& b) { return a.pt > b.pt; });
        meta.l2flavor = 1; meta.l2_index = selected_muons[0].index;
        meta.l3flavor = 1; meta.l3_index = selected_muons[1].index;
        meta.l4flavor = 1; meta.l4_index = selected_muons[2].index;
    }
    return true;
}

// Old Z to ll selection
std::string ZToLLSelection::name() const { return "Z_to_ll"; }
bool ZToLLSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    double bestDiff = 1e9;
    int best_i = -1, best_j = -1, best_flav = -1; // 0=e,1=mu
    double bestMass = std::numeric_limits<double>::quiet_NaN();

    // electrons
    for (int i = 0; i < evt.d->Electron_size; ++i) {
        if (evt.d->Electron_PT[i] < cfg.zl_pt_min) continue; // pT cut
        for (int j = i+1; j < evt.d->Electron_size; ++j) {
            if (evt.d->Electron_PT[j] < cfg.zl_pt_min) continue; // pT cut
            if (evt.d->Electron_Charge[i]*evt.d->Electron_Charge[j] >= 0) continue; // OS
            TLorentzVector l1, l2;
            l1.SetPtEtaPhiM(evt.d->Electron_PT[i], evt.d->Electron_Eta[i], evt.d->Electron_Phi[i], Me);
            l2.SetPtEtaPhiM(evt.d->Electron_PT[j], evt.d->Electron_Eta[j], evt.d->Electron_Phi[j], Me);
            double mass = (l1+l2).M();
            double diff = mass - cfg.z_mass;
            if (diff < cfg.z_mass_window_upper && diff > -cfg.z_mass_window_lower && diff < bestDiff) {
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
            double diff = mass - cfg.z_mass;
            if (diff < cfg.z_mass_window_upper && diff > -cfg.z_mass_window_lower && diff < bestDiff) {
                bestDiff = diff; best_i = i; best_j = j; best_flav = 1; bestMass = mass;
            }
        }
    }

    if (best_flav == -1) return false;
    meta.z_l1 = best_i; meta.z_l2 = best_j; meta.z_flavor = best_flav;
    meta.z_mass = bestMass; meta.z_mass_diff = bestDiff;
    return true;
}




/*
    II) Z candidate selection:
        1) Identify flavor of Z candidate (ee or mumu) from the selected leptons
        - if 1 muon + 3 electrons: Z->ee
        - if 1 electron + 3 muons: Z->mumu
        2) from 3 leptons:
            - select one that have same sign as the prompt lepton from H decay -> first lepton candidate
        3) Loop based on the Z from 2)
            - check OS for the two leptons from Z
            - compute invariant mass
            - check mass window
    Meta info to store:
        - z_l1 index
        - z_l2 index
        - z_flavor (0=e,1=mu)
        - z_mass
        - z_mass_diff
*/

std::string ZCandidateSelection::name() const { return "ZCandidateSelection"; }
bool ZCandidateSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    // determine Z flavor
    int z_flav = -1; // 0=e,1=mu
    if (meta.l1flavor == 1 && meta.l2flavor == 0 && meta.l3flavor == 0 && meta.l4flavor == 0) {
        z_flav = 0; // Z->ee
    } else if (meta.l1flavor == 0 && meta.l2flavor == 1 && meta.l3flavor == 1 && meta.l4flavor == 1) {
        z_flav = 1; // Z->mumu
    } else {
        return false; // invalid flavor combination
    }
    // identify first lepton candidate from Z
    int first_z_lep_index = -1;
    // get charge of prompt lepton from H decay
    int prompt_lep_charge = 0;
    if (meta.l1flavor == 0) {
        prompt_lep_charge = evt.d->Electron_Charge[meta.l1_index];
    } else if (meta.l1flavor == 1) {
        prompt_lep_charge = evt.d->Muon_Charge[meta.l1_index];
    }
    // identify first lepton candidate from Z (same sign as prompt lepton)
    if (z_flav == 0) {
        // Z->ee
        for (int i = 2; i <=4; ++i) {
            int idx = (i == 2) ? meta.l2_index : (i == 3) ? meta.l3_index : meta.l4_index;
            if (evt.d->Electron_Charge[idx] == prompt_lep_charge) {
                first_z_lep_index = idx;
                break;
            }
        }
    } else if (z_flav == 1) {
        // Z->mumu
        for (int i = 2; i <=4; ++i) {
            int idx = (i == 2) ? meta.l2_index : (i == 3) ? meta.l3_index : meta.l4_index;
            if (evt.d->Muon_Charge[idx] == prompt_lep_charge) {
                first_z_lep_index = idx;
                break;
            }
        }
    }
    if (first_z_lep_index == -1) return false; // could not find
    
    int second_z_lep_index = -1;
    double bestDiff = 1e9;
    double bestMass = std::numeric_limits<double>::quiet_NaN();
    if (z_flav == 0) {
        // Z->ee
        for (int i = 2; i <=4; ++i) {
            int idx = (i == 2) ? meta.l2_index : (i == 3) ? meta.l3_index : meta.l4_index;
            if (idx == first_z_lep_index) continue;
            // check OS
            if (evt.d->Electron_Charge[first_z_lep_index] * evt.d->Electron_Charge[idx] >= 0) continue;
            // compute invariant mass
            TLorentzVector l1, l2;
            l1.SetPtEtaPhiM(evt.d->Electron_PT[first_z_lep_index], evt.d->Electron_Eta[first_z_lep_index],
                            evt.d->Electron_Phi[first_z_lep_index], Me);
            l2.SetPtEtaPhiM(evt.d->Electron_PT[idx], evt.d->Electron_Eta[idx],
                            evt.d->Electron_Phi[idx], Me);
            double mass = (l1 + l2).M();
            double diff = mass - cfg.z_mass;
            if (diff < cfg.z_mass_window_upper && diff > -cfg.z_mass_window_lower) {
                if (diff < bestDiff) {
                    bestDiff = diff;
                    bestMass = mass;
                    second_z_lep_index = idx;
                }
            }
        }
    }
    if (z_flav == 1) {
        // Z->mumu
        for (int i = 2; i <=4; ++i) {
            int idx = (i == 2) ? meta.l2_index : (i == 3) ? meta.l3_index : meta.l4_index;
            if (idx == first_z_lep_index) continue;
            // check OS
            if (evt.d->Muon_Charge[first_z_lep_index] * evt.d->Muon_Charge[idx] >= 0) continue;
            // compute invariant mass
            TLorentzVector l1, l2;
            l1.SetPtEtaPhiM(evt.d->Muon_PT[first_z_lep_index], evt.d->Muon_Eta[first_z_lep_index],
                            evt.d->Muon_Phi[first_z_lep_index], Mmu);
            l2.SetPtEtaPhiM(evt.d->Muon_PT[idx], evt.d->Muon_Eta[idx],
                            evt.d->Muon_Phi[idx], Mmu);
            double mass = (l1 + l2).M();
            double diff = mass - cfg.z_mass;
            if (diff < cfg.z_mass_window_upper && diff > -cfg.z_mass_window_lower) {
                if (diff < bestDiff) {
                    bestDiff = diff;
                    bestMass = mass;
                    second_z_lep_index = idx;
                }
            }
        }   
    }
    if (second_z_lep_index != -1) {
        meta.z_l1 = first_z_lep_index;
        meta.z_l2 = second_z_lep_index;
        meta.z_flavor = z_flav;
        meta.z_mass = bestMass;
        meta.z_mass_diff = bestDiff;
        // return true;
    
    
        // Calculate the invariant mass of Z case 1 (use best pair)
        meta.m_z1 = meta.z_mass;
        // Calculate the invariant mass of Z case 2 (use 1st z candidate + second best lepton)
        // e.g., e1 (z 1st cand), e2 (z best pair), e3 (remaining) -> use e1 + e3 to calculate m_z2
        int remaining_lep_index = -1;
        if (z_flav == 0) {
            // Z->ee
            for (int i = 2; i <=4; ++i) {
                int idx = (i == 2) ? meta.l2_index : (i == 3) ? meta.l3_index : meta.l4_index;
                if (idx != meta.z_l1 && idx != meta.z_l2) {
                    remaining_lep_index = idx;
                    break;
                }
            }
            if (remaining_lep_index != -1) {
                TLorentzVector l1, l2;
                l1.SetPtEtaPhiM(evt.d->Electron_PT[meta.z_l1], evt.d->Electron_Eta[meta.z_l1],
                                evt.d->Electron_Phi[meta.z_l1], Me);
                l2.SetPtEtaPhiM(evt.d->Electron_PT[remaining_lep_index], evt.d->Electron_Eta[remaining_lep_index],
                                evt.d->Electron_Phi[remaining_lep_index], Me);
                meta.m_z2 = (l1 + l2).M();
            }
            meta.h_e_pt = evt.d->Electron_PT[remaining_lep_index];
            meta.h_mu_pt = evt.d->Muon_PT[meta.l1_index];

        } else if (z_flav == 1) {
            // Z->mumu
            for (int i = 2; i <=4; ++i) {
                int idx = (i == 2) ? meta.l2_index : (i == 3) ? meta.l3_index : meta.l4_index;
                if (idx != meta.z_l1 && idx != meta.z_l2) {
                    remaining_lep_index = idx;
                    break;
                }
            }
            if (remaining_lep_index != -1) {
                TLorentzVector l1, l2;
                l1.SetPtEtaPhiM(evt.d->Muon_PT[meta.z_l1], evt.d->Muon_Eta[meta.z_l1],
                                evt.d->Muon_Phi[meta.z_l1], Mmu);
                l2.SetPtEtaPhiM(evt.d->Muon_PT[remaining_lep_index], evt.d->Muon_Eta[remaining_lep_index],
                                evt.d->Muon_Phi[remaining_lep_index], Mmu);
                meta.m_z2 = (l1 + l2).M();
            }
            meta.h_mu_pt = evt.d->Muon_PT[remaining_lep_index];
            meta.h_e_pt = evt.d->Electron_PT[meta.l1_index];
        }
        return true;

    }
    return false;
}

std::string HToMuTauESelection::name() const { return "H_to_mutau_e"; }
bool HToMuTauESelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
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

    meta.h_mu_pt = evt.d->Muon_PT[idx_mu];
    meta.h_e_pt = evt.d->Electron_PT[idx_e];

    // precompute dphi_e_met
    if (evt.d->MissingET_size > 0) {
        meta.dphi_e_met = deltaPhiFromPhis(evt.d->Electron_Phi[idx_e], evt.d->MissingET_Phi[0]);
    }
    // precompute dphi_mu_met
    if (evt.d->MissingET_size > 0) {
        meta.dphi_mu_met = deltaPhiFromPhis(evt.d->Muon_Phi[idx_mu], evt.d->MissingET_Phi[0]);
    }
    // DeltaPhi(mu, e)
    if (evt.d->Muon_size > 0 && evt.d->Electron_size > 0) {
        meta.dphi_mu_e = deltaPhiFromPhis(evt.d->Muon_Phi[idx_mu], evt.d->Electron_Phi[idx_e]);
    }
    // collinear mass
    if (evt.d->MissingET_size > 0) {
        meta.m_collinear = computeCollinearMassMuTauE(evt, idx_mu, idx_e);
        meta.m_transverse_e = computeTransverseMass(evt, 0, idx_e);
        meta.m_transverse_mu = computeTransverseMass(evt, 1, idx_mu);
    }

    // calculate collinear mass in 2 scenarios
    // Scenario 1: consider Z candidate with closest mass to m_Z a Z's leptons-> same as above
    meta.m_h1 = meta.m_collinear;
    // Scenario 2: consider second best Z candidate's leptons as Z's leptons -> causing the H collinear mass to be
    // calculated using the remaining leptons (h_mu, zl2, MET) or m_h2 = collinear mass of system (h_mu, zl2, MET)
    if (meta.z_flavor == 0) {
        // Z->ee, remaining lepton is muon
        meta.m_h2 = computeCollinearMassMuTauE(evt, idx_mu, zl2);
    } else if (meta.z_flavor == 1) {
        // Z->mumu, remaining lepton is electron
        meta.m_h2 = computeCollinearMassMuTauE(evt, zl2, idx_e);
    }
    
    return true;
}

// H->e tau_mu selection
std::string HToETauMuSelection::name() const { return "H_to_etau_mu"; }
bool HToETauMuSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    // Implementation similar to HToMuESelection but for e and tau_mu
    int zl1 = meta.z_l1, zl2 = meta.z_l2;

    // Handle for invalid indices
    if (zl1 < 0 || zl2 < 0) return false;

    int idx_tau_mu = -1;
    for (int m = 0; m < evt.d->Muon_size; ++m) {
        if (meta.z_flavor == 1 && (m == zl1 || m == zl2)) continue; // exclude Z muons
        if (evt.d->Muon_PT[m] > cfg.mu_pt_min) {
            if (idx_tau_mu == -1) idx_tau_mu = m; else return false; // more than one
        }
    }
    if (idx_tau_mu < 0) return false;

    int idx_e = -1;
    for (int e = 0; e < evt.d->Electron_size; ++e) {
        if (meta.z_flavor == 0 && (e == zl1 || e == zl2)) continue; // exclude Z electrons
        if (evt.d->Electron_PT[e] > cfg.e_pt_min) {
            if (idx_e == -1) idx_e = e; else return false; // more than one
        }
    }
    if (idx_e < 0) return false;

    if (evt.d->Muon_Charge[idx_tau_mu] * evt.d->Electron_Charge[idx_e] >= 0) return false; // OS

    meta.h_mu = idx_tau_mu;
    meta.h_e = idx_e;

    meta.h_mu_pt = evt.d->Muon_PT[idx_tau_mu];
    meta.h_e_pt = evt.d->Electron_PT[idx_e];

    // precompute dphi_e_met
    if (evt.d->MissingET_size > 0) {
        meta.dphi_e_met = deltaPhiFromPhis(evt.d->Electron_Phi[idx_e], evt.d->MissingET_Phi[0]);
    }
    // precompute dphi_mu_met
    if (evt.d->MissingET_size > 0) {
        meta.dphi_mu_met = deltaPhiFromPhis(evt.d->Muon_Phi[idx_tau_mu], evt.d->MissingET_Phi[0]);
    }
    // DeltaPhi(mu, e)
    if (evt.d->Muon_size > 0 && evt.d->Electron_size > 0) {
        meta.dphi_mu_e = deltaPhiFromPhis(evt.d->Muon_Phi[idx_tau_mu], evt.d->Electron_Phi[idx_e]);
    }
    // collinear mass
    if (evt.d->MissingET_size > 0) {
        meta.m_collinear = computeCollinearMassETauMu(evt, idx_e, idx_tau_mu);
        meta.m_transverse_e = computeTransverseMass(evt, 0, idx_e);  // wrong calculate transverse mass for electron
        meta.m_transverse_mu = computeTransverseMass(evt, 1, idx_tau_mu);
    }

    // calculate collinear mass in 2 scenarios
    // Scenario 1: consider Z candidate with closest mass to m_Z a Z's leptons-> same as above
    meta.m_h1 = meta.m_collinear;
    // Scenario 2: consider second best Z candidate's leptons as Z's leptons -> causing the H collinear mass to be
    // calculated using the remaining leptons (h_e, zl2, MET) or m_h2 = collinear mass of system (h_e, zl2, MET)
    if (meta.z_flavor == 0) {
        // Z->ee, remaining lepton is muon
        meta.m_h2 = computeCollinearMassETauMu(evt, zl2, idx_tau_mu);
    } else if (meta.z_flavor == 1) {
        // Z->mumu, remaining lepton is electron
        meta.m_h2 = computeCollinearMassETauMu(evt, idx_e, zl2);
    }

    return true;
}

std::string METEDphiSelection::name() const { return "MET_E_dphi"; }
bool METEDphiSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    if (meta.h_e < 0) return false; // need H electron selected
    if (evt.d->MissingET_size <= 0) return false;
    // if (std::isnan(meta.dphi_e_met)) {
    //     TLorentzVector v1, v2;
    //     v1.SetPtEtaPhiM(1.0, 0.0, evt.d->Electron_Phi[meta.h_e], 0.0);
    //     v2.SetPtEtaPhiM(1.0, 0.0, evt.d->MissingET_Phi[0], 0.0);
    //     meta.dphi_e_met = std::abs(v1.DeltaPhi(v2));
    // }
    // return meta.dphi_e_met < cfg.max_dphi_e_met;
    // Use precomputed value
    if (std::isnan(meta.dphi_e_met)) return false;
    return meta.dphi_e_met < cfg.max_dphi_e_met;
}

std::string METMuDphiSelection::name() const { return "MET_Mu_dphi"; }
bool METMuDphiSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    if (meta.h_mu < 0) return false; // need H muon selected
    if (evt.d->MissingET_size <= 0) return false;
    if (std::isnan(meta.dphi_mu_met)) return false;
    return meta.dphi_mu_met < cfg.max_dphi_mu_met;
}

} // namespace hlfv

