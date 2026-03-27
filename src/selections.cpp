#include "../include/selections.h"

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <TLorentzVector.h>
#include <TVector3.h>

namespace hlfv {

// --- Constants ---
constexpr double MASS_E   = 0.000511;
constexpr double MASS_MU  = 0.10565837;
constexpr double ECM      = 240.0; // Center of mass energy
constexpr double Z_MASS   = 91.1876;

enum Flavor { ELECTRON = 0, MUON = 1 };

struct LeptonObj {
    int index;
    int flavor; 
    double pt;
    int charge;
};

// --- Helper Functions ---
namespace {

    // Generic P4 retriever
    TLorentzVector get_p4(const Event& evt, int idx, int flavor) {
        TLorentzVector v;
        if (flavor == ELECTRON) {
            v.SetPtEtaPhiM(evt.d->Electron_PT[idx], evt.d->Electron_Eta[idx], 
                           evt.d->Electron_Phi[idx], MASS_E);
        } else {
            v.SetPtEtaPhiM(evt.d->Muon_PT[idx], evt.d->Muon_Eta[idx], 
                           evt.d->Muon_Phi[idx], MASS_MU);
        }
        return v;
    }

    double delta_phi(double phi1, double phi2) {
        double dphi = std::abs(phi1 - phi2);
        if (dphi > M_PI) dphi = 2 * M_PI - dphi;
        return dphi;
    }

    // Generic Collinear Mass Calculation
    // Calculates mass of (visible_lep + invisible_tau_products)
    // Assumes neutrinos go in direction of visible tau product (approx)
    double compute_collinear_mass(const Event& evt, const TLorentzVector& v_tau_vis, const TLorentzVector& v_other) {
        if (!evt.d || evt.d->MissingET_size <= 0) return std::numeric_limits<double>::quiet_NaN();
        if (v_tau_vis.Pt() <= 0) return std::numeric_limits<double>::quiet_NaN();

        double met = evt.d->MissingET_MET[0];
        double met_phi = evt.d->MissingET_Phi[0];
        double met_px = met * std::cos(met_phi);
        double met_py = met * std::sin(met_phi);

        // Project MET onto the visible tau direction
        double dir_x = v_tau_vis.Px() / v_tau_vis.Pt();
        double dir_y = v_tau_vis.Py() / v_tau_vis.Pt();
        double proj = met_px * dir_x + met_py * dir_y;

        if (proj <= 0) return std::numeric_limits<double>::quiet_NaN();

        // Visible fraction of tau momentum
        double x_tau_vis = v_tau_vis.Pt() / (v_tau_vis.Pt() + proj);

        if (x_tau_vis <= 0 || x_tau_vis > 1) return std::numeric_limits<double>::quiet_NaN();

        double m_vis = (v_tau_vis + v_other).M();
        return m_vis / std::sqrt(x_tau_vis);
    }

    double compute_transverse_mass(double pt, double phi, double met, double met_phi) {
        double dphi = delta_phi(phi, met_phi);
        return std::sqrt(2 * pt * met * (1 - std::cos(dphi)));
    }
// Post-calculation routine to compute recoil, alternative masses, and impact parameters
    void compute_z_post_calculations(const Event& evt, Meta& meta, const Parameters& cfg) {
        
        // 1. Invariant Mass of the Selected Z (Case 1)
        // Already computed as meta.z_mass, but we assign to m_z1 for consistency
        meta.m_z1 = meta.z_mass;

        // 2. Invariant Mass of the Alternative Pair (Case 2)
        // Pair: Z_L1 (fixed L2) + The remaining OS lepton (other_idx)
        if (meta.otherZ_idx != -1) {
            TLorentzVector p4_z1 = get_p4(evt, meta.z_l1, meta.z_flavor);
            TLorentzVector p4_other = get_p4(evt, meta.otherZ_idx, meta.z_flavor);
            meta.m_z2 = (p4_z1 + p4_other).M();
        }

        // 3. Recoil Mass Calculation (against Beam)
        // Based on the SELECTED Z candidate
        TLorentzVector p4_z_selected = get_p4(evt, meta.z_l1, meta.z_flavor) + 
                                       get_p4(evt, meta.z_l2, meta.z_flavor);

        double beam_E = ECM / 2.0; 
        TLorentzVector p_beam_total(0.0, 0.0, 0.0, ECM); // px, py, pz, E (approx for 0 crossing angle)
        
        TLorentzVector p_recoil = p_beam_total - p4_z_selected;
        meta.m_recoil = p_recoil.M();

        // 4. Boost Vectors (Recoil Frame)
        TVector3 boost_vec = -p_recoil.BoostVector();
        meta.beta_x = boost_vec.X();
        meta.beta_y = boost_vec.Y();
        meta.beta_z = boost_vec.Z();

        // 5. Impact Parameters (D0/DZ) for Higgs Candidates
        // The Higgs candidates are: 
        //   A. L1 (The singleton, prompt lepton)
        //   B. The 'other' lepton (The one NOT selected for Z)
        
        // Helper to set D0/DZ safely
        auto set_d0dz = [&](int idx, int flavor) {
            double d0 = (flavor == ELECTRON) ? evt.d->Electron_D0[idx] : evt.d->Muon_D0[idx];
            double dz = (flavor == ELECTRON) ? evt.d->Electron_DZ[idx] : evt.d->Muon_DZ[idx];
            
            if (flavor == ELECTRON) {
                meta.h_e_d0 = std::abs(d0);
                meta.h_e_dz = std::abs(dz);
            } else {
                meta.h_mu_d0 = std::abs(d0);
                meta.h_mu_dz = std::abs(dz);
            }
        };

        // Set for L1 (Prompt)
        set_d0dz(meta.l1_index, meta.l1flavor);

        // Set for 'Other' (From Z triplet, but assigned to Higgs)
        // This has the same flavor as the Z pair
        if (meta.otherZ_idx != -1) {
            set_d0dz(meta.otherZ_idx, meta.z_flavor);
        }

        // Collinear mass calculation
        
        // Higgs
        if (meta.otherH_idx != -1) {
            // 1. Calculate Correct Pair (h_e, h_mu)
            TLorentzVector p4_e  = get_p4(evt, meta.h_e, ELECTRON);
            TLorentzVector p4_mu = get_p4(evt, meta.h_mu, MUON);

            // Lambda or simple ternary macro prevents repeating the mode logic twice
            auto get_mass = [&](TLorentzVector& e, TLorentzVector& mu) {
                return compute_collinear_mass(evt, 
                    (cfg.mode == 0) ? mu : e,  // v_tau (mode 0: mu is tau)
                    (cfg.mode == 0) ? e  : mu  // v_other
                );
            };

            meta.m_h1 = get_mass(p4_e, p4_mu);

            // 2. Calculate Incorrect Pair (Swap one index)
            // Update the specific vector that corresponds to the "other" flavor
            if (meta.l1flavor == MUON) p4_e  = get_p4(evt, meta.otherH_idx, ELECTRON);
            else                       p4_mu = get_p4(evt, meta.otherH_idx, MUON);

            meta.m_h2 = get_mass(p4_e, p4_mu);
        }

        // Z boson
        if (meta.otherZ_idx != -1) {
            TLorentzVector p4_z1 = get_p4(evt, meta.z_l1, meta.z_flavor);

            // 1. Correct Pair
            TLorentzVector p4_z2 = get_p4(evt, meta.z_l2, meta.z_flavor);
            meta.m_z1 = (p4_z1 + p4_z2).M();

            // 2. Incorrect Pair
            TLorentzVector p4_other = get_p4(evt, meta.otherZ_idx, meta.z_flavor);
            meta.m_z2 = (p4_z1 + p4_other).M();
        }

    }

} // anonymous namespace

// --- Selections ---

std::string EmptySelection::name() const { return "EmptySelection"; }
bool EmptySelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    // Example test
    if (meta.m_transverse_e >= 2.0) return false;
    return true;
}

std::string FinalState_NoCut::name() const { return "FinalState_NoCut"; }
bool FinalState_NoCut::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    return true;
}

// -----------------------------------------------------------------------------
// Lepton Selection
// -----------------------------------------------------------------------------
std::string LeptonSelection::name() const { return "LeptonSelection"; }

bool LeptonSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;

    // save met
    if (evt.d->MissingET_size > 0) {
        meta.MET = evt.d->MissingET_MET[0];
    }

    // ----- Stricter all leptons must equal 4 ------
    if ((evt.d->Electron_size + evt.d->Muon_size) != 4) return false;

    // Helper to extract passing leptons
    auto extract = [&](int n, const float* pts, const int* charges, int flav) {
        std::vector<LeptonObj> leps;
        for (int i = 0; i < n; ++i) {
            if (pts[i] > cfg.lepton_pt_min) leps.push_back({i, flav, (double)pts[i], charges[i]});
        }
        return leps;
    };

    auto electrons = extract(evt.d->Electron_size, evt.d->Electron_PT, evt.d->Electron_Charge, ELECTRON);
    auto muons     = extract(evt.d->Muon_size, evt.d->Muon_PT, evt.d->Muon_Charge, MUON);

    // Topology: (1e + 3mu) OR (1mu + 3e)
    bool is_1e_3mu = (electrons.size() == 1 && muons.size() == 3);
    bool is_1mu_3e = (muons.size() == 1 && electrons.size() == 3);

    if (!is_1e_3mu && !is_1mu_3e) return false;

    // Net Charge Check
    int net_charge = 0;
    for(auto& l : electrons) net_charge += l.charge;
    for(auto& l : muons)     net_charge += l.charge;
    if (net_charge != 0) return false;

    // Identify Singleton (L1) and Triplet (L2, L3, L4)
    // If 1mu+3e: Singleton is muons, Triplet is electrons.
    const auto& singleton = is_1mu_3e ? muons : electrons;
    const auto& triplet   = is_1mu_3e ? electrons : muons;

    // Assign L1 (The unique flavor lepton)
    const auto& l1 = singleton[0];
    meta.l1flavor = l1.flavor;
    meta.l1_index = l1.index;

    // Sort Triplet based on charge relative to L1
    // Logic: If Net Charge is 0, and we have 1 vs 3.
    // Let L1 charge be Q. 
    // The triplet must contain one lepton with charge Q (SS) and two with -Q (OS).
    // Or else net charge wouldn't be 0.
    std::vector<const LeptonObj*> ss_leptons; // Same Sign as L1
    std::vector<const LeptonObj*> os_leptons; // Opposite Sign to L1

    for (const auto& l : triplet) {
        if (l.charge == l1.charge) ss_leptons.push_back(&l);
        else                       os_leptons.push_back(&l);
    }

    if (ss_leptons.size() != 1 || os_leptons.size() != 2) return false;

    // Assign L2 (The SS lepton) -> this must be z boson candidate
    meta.l2flavor = ss_leptons[0]->flavor;
    meta.l2_index = ss_leptons[0]->index;

    // Sort OS leptons by pT (Descending)
    std::sort(os_leptons.begin(), os_leptons.end(), [](const LeptonObj* a, const LeptonObj* b) {
        return a->pt > b->pt;
    });

    // Assign L3, L4 (The OS leptons)
    meta.l3flavor = os_leptons[0]->flavor; meta.l3_index = os_leptons[0]->index;
    meta.l4flavor = os_leptons[1]->flavor; meta.l4_index = os_leptons[1]->index;

    return true;
}

// -----------------------------------------------------------------------------
// Z Candidate Selection
// -----------------------------------------------------------------------------

// std::string ZCandidateSelection::name() const { return "ZCandidateSelection"; }
// bool ZCandidateSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
//     if (!evt.d) return false;

//     // 1. Determine Z flavor (inverse of L1 flavor)
//     int z_flav = (meta.l1flavor == ELECTRON) ? MUON : ELECTRON;

//     // 2. L2 is GUARANTEED to be the first Z candidate 
//     // (It is the only lepton in the triplet with the same sign as L1)
//     int z1_idx = meta.l2_index;
    
//     // 3. Test L2 against L3 and L4 to find the best Z partner
//     // We only need to check the OS candidates (L3, L4)
//     int candidates[2] = {meta.l3_index, meta.l4_index};
    
//     int best_z2_idx  = -1;
//     int h2_idx       = -1;
//     int otherZ_idx   = -1; // The one NOT selected
//     int otherH_idx   = -1; // The one NOT selected (for Higgs)
//     double best_diff = 1e9;
//     double best_mass = std::numeric_limits<double>::quiet_NaN();

//     // Get L2 P4 once
//     TLorentzVector p4_z1 = get_p4(evt, z1_idx, z_flav);

//     for (int idx : candidates) {
//         TLorentzVector p4_z2 = get_p4(evt, idx, z_flav);
//         double mass = (p4_z1 + p4_z2).M();
//         double diff = mass - cfg.z_mass;

//         if (diff < cfg.z_mass_window_upper && diff > -cfg.z_mass_window_lower) {
//             if (diff < best_diff) {
//                 best_diff = diff;
//                 best_mass = mass;
//                 best_z2_idx = idx;
//                 h2_idx = (idx == candidates[0]) ? candidates[1] : candidates[0];
//                 // Alternative pair; this serves as check if our selection is correct
//                 otherZ_idx = h2_idx; // swap
//                 otherH_idx = idx;
//             }
//         }
//     }

//     if (best_z2_idx != -1) {
//         meta.z_l1 = z1_idx;
//         meta.z_l2 = best_z2_idx;
//         meta.z_flavor = z_flav;
//         meta.z_mass = best_mass;
//         meta.z_mass_diff = best_diff;
//         meta.otherZ_idx = otherZ_idx;
//         meta.otherH_idx = otherH_idx;
//         meta.h_e = (meta.l1flavor == ELECTRON) ? meta.l1_index : h2_idx;
//         meta.h_mu = (meta.l1flavor == MUON) ? meta.l1_index : h2_idx;
//         meta.h_e_pt = evt.d->Electron_PT[meta.h_e];
//         meta.h_mu_pt = evt.d->Muon_PT[meta.h_mu];

//         // Post-calculations
//         compute_z_post_calculations(evt, meta, cfg);

//         return true;
//     }

//     return false;
// }

// Hard coded for two signal regions [21, 81) and [81, 101) GeV
// If two pairs fall within the two windows, assign the one in [81, 101] GeV as the Z candidate

std::string ZCandidateSelection::name() const { return "ZCandidateSelection"; }
bool ZCandidateSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;

    // 1. Determine Z flavor (inverse of L1 flavor)
    int z_flav = (meta.l1flavor == ELECTRON) ? MUON : ELECTRON;

    // 2. L2 is GUARANTEED to be the first Z candidate 
    // (It is the only lepton in the triplet with the same sign as L1)
    int z1_idx = meta.l2_index;
    
    // 3. Test L2 against L3 and L4 to find the best Z partner
    // We only need to check the OS candidates (L3, L4)
    int candidates[2] = {meta.l3_index, meta.l4_index};
    
    int best_z2_idx  = -1;
    int h2_idx       = -1;
    int otherZ_idx   = -1; // The one NOT selected
    int otherH_idx   = -1; // The one NOT selected (for Higgs)
    double best_diff = 1e9;
    double best_mass = std::numeric_limits<double>::quiet_NaN();

    bool use_window_21_81 = false;
    if (cfg.z_mass == 81 && cfg.z_mass_window_lower == 60 && cfg.z_mass_window_upper == 0) {
        use_window_21_81 = true;
    }
    // Get L2 P4 once
    TLorentzVector p4_z1 = get_p4(evt, z1_idx, z_flav);

    for (int idx : candidates) {
        TLorentzVector p4_z2 = get_p4(evt, idx, z_flav);
        double mass = (p4_z1 + p4_z2).M();
        double diff = mass - cfg.z_mass;

        if (diff < cfg.z_mass_window_upper && diff > -cfg.z_mass_window_lower) {
            if (std::abs(diff) < best_diff) {
                best_diff = std::abs(diff);
                best_mass = mass;
                best_z2_idx = idx;
                h2_idx = (idx == candidates[0]) ? candidates[1] : candidates[0];
                // Alternative pair; this serves as check if our selection is correct
                otherZ_idx = h2_idx; // swap
                otherH_idx = idx;

                if (use_window_21_81) {
                    double alt_pair_mass = (p4_z1 + get_p4(evt, otherZ_idx, z_flav)).M();
                    if (alt_pair_mass > 81.0 && alt_pair_mass < 101.0) {
                        // If the alternative pair falls within the Z mass window, reject the event
                        // debug print
                        return false;
                    }
                }
            }
        }
    }

    if (best_z2_idx != -1) {
        meta.z_l1 = z1_idx;
        meta.z_l2 = best_z2_idx;
        meta.z_flavor = z_flav;
        meta.z_mass = best_mass;
        meta.z_mass_diff = best_diff;
        meta.otherZ_idx = otherZ_idx;
        meta.otherH_idx = otherH_idx;
        meta.h_e = (meta.l1flavor == ELECTRON) ? meta.l1_index : h2_idx;
        meta.h_mu = (meta.l1flavor == MUON) ? meta.l1_index : h2_idx;
        meta.h_e_pt = evt.d->Electron_PT[meta.h_e];
        meta.h_mu_pt = evt.d->Muon_PT[meta.h_mu];

        // Post-calculations
        compute_z_post_calculations(evt, meta, cfg);

        return true;
    }

    return false;
}


// -----------------------------------------------------------------------------
// H -> Mu E selection (for mH < 150 GeV) (no MET in signal)
// - dR(z leptons) < 2.0
// - Higgs m_inv > 75 GeV (no MET)
// -----------------------------------------------------------------------------
std::string HToMuESelection::name() const { return "H_to_mu_e"; }
bool HToMuESelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    
    int zl1 = meta.z_l1;
    int zl2 = meta.z_l2;
    int idx_mu = meta.h_mu;
    int idx_e  = meta.h_e;
    // OS Check
    if (evt.d->Muon_Charge[idx_mu] * evt.d->Electron_Charge[idx_e] >= 0) {
        // std::cout << "HToMuESelection: Failed OS check between mu and e." << std::endl;
        return false;
    }

    // PT requirements
    if (evt.d->Muon_PT[idx_mu] < cfg.mu_pt_min || evt.d->Electron_PT[idx_e] < cfg.e_pt_min) {
        // std::cout << "HToMuESelection: Failed pT requirement. mu_pt=" << evt.d->Muon_PT[idx_mu] 
        //           << ", e_pt=" << evt.d->Electron_PT[idx_e] << std::endl;
        return false;
    }

    // dR of Higgs pair
    // TLorentzVector p4_mu = get_p4(evt, idx_mu, MUON);
    // TLorentzVector p4_e  = get_p4(evt, idx_e, ELECTRON);
    // double deltaR_mu_e = p4_mu.DeltaR(p4_e);
    // meta.deltaR_mu_e = deltaR_mu_e;

    // dR Check between Z leptons
    // TLorentzVector p4_zl1 = get_p4(evt, zl1, meta.z_flavor);
    // TLorentzVector p4_zl2 = get_p4(evt, zl2, meta.z_flavor);
    // double deltaR = p4_zl1.DeltaR(p4_zl2);
    // if (deltaR >= 2.0) {
    //     std::cout << "HToMuESelection: Failed dR check between Z leptons. deltaR=" << deltaR << std::endl;
    //     std::cout << "  Z1 (idx=" << zl1 << ", pt=" << p4_zl1.Pt() << ", eta=" << p4_zl1.Eta() << ", phi=" << p4_zl1.Phi() << ")" << std::endl;
    //     std::cout << "  Z2 (idx=" << zl2 << ", pt=" << p4_zl2.Pt() << ", eta=" << p4_zl2.Eta() << ", phi=" << p4_zl2.Phi() << ")" << std::endl;
    //     return false;
    // }

    // MET cut < 10 GeV
    if (evt.d->MissingET_MET[0] > 10.0) {
        // std::cout << "HToMuESelection: Failed MET cut. MET=" << evt.d->MissingET_MET[0] << std::endl;
        return false;
    }

    // Invariant Mass of Higgs Candidate (mu + e)
    TLorentzVector p4_mu = get_p4(evt, idx_mu, MUON);
    TLorentzVector p4_e  = get_p4(evt, idx_e, ELECTRON);

    double m_h_inv = (p4_mu + p4_e).M();
    meta.m_h_invariant = m_h_inv;
    // if (m_h_inv < 75.0) return false;

    return true;
}

// -----------------------------------------------------------------------------
// H -> Mu E with PT selection first
std::string PT_ZCandidateSelection::name() const { return "PT_ZCandidateSelection"; }
bool PT_ZCandidateSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    // evt.d delphes class pointer
    // meta contains the calculated information
    // cfg contains selection information
    // l1 is the singleton, l2 is Z-prompt, l3,l4 are the triplet
    // z_l1 is SCDF of l1flavor, z_l2 is the best Z prompt, z_flavor is their flavor

    // ---------- Prepare p4 ----------
    meta.z_flavor = (meta.l1flavor == ELECTRON) ? MUON : ELECTRON;
    TLorentzVector h1_p4 = get_p4(evt, meta.l1_index, meta.l1flavor);
    // ---------- PT cut on the singleton ----------
    if (h1_p4.Pt() <= (meta.l1flavor == ELECTRON ? cfg.e_pt_min : cfg.mu_pt_min)) {return false;}
    TLorentzVector z1_p4 = get_p4(evt, meta.l2_index, meta.l2flavor);
    vector<TLorentzVector> z_candidates;
    vector<float> z_candidate_masses;
    vector<int> PhSp_21t81_candidateIdx, PhSp_81t101_candidateIdx;
    // ---------- Set vectors for Z candidates and their masses ----------
    z_candidates.push_back(get_p4(evt, meta.l3_index, meta.l3flavor));
    z_candidates.push_back(get_p4(evt, meta.l4_index, meta.l4flavor));
    z_candidate_masses.push_back((z1_p4 + z_candidates[0]).M());
    z_candidate_masses.push_back((z1_p4 + z_candidates[1]).M());
    // ---------- Looking at phase space ----------
    bool use_window_21_81 = (cfg.z_mass == 81 && cfg.z_mass_window_lower == 60 && cfg.z_mass_window_upper == 0);
    for (int vector_idx = 0; vector_idx < z_candidates.size(); vector_idx++) {
        if (z_candidate_masses[vector_idx] > 81.0 && z_candidate_masses[vector_idx] < 101.0) {
            if (use_window_21_81) return false; // If the alternative pair falls within the Z mass window, reject the event
            PhSp_81t101_candidateIdx.push_back(vector_idx);
        }
        else if (z_candidate_masses[vector_idx] > 21.0 && z_candidate_masses[vector_idx] < 81.0) {
            PhSp_21t81_candidateIdx.push_back(vector_idx);
        }
    }
    // Logic here:
    // 1. if 21t81 : the event will not pass to this section if there is a candidate in 81t101. so the loop will
    // only run for candidates in 21t81, and assign the best one as Z candidate
    // 2. if 81t101: the loop will only run for candidates in 81t101, and assign the best one as Z candidate, no matter if there is a candidate in 21t81 or not
    // 3. This two phase space is completely non-overlapping, so no need to worry about the case where there are candidates in both windows.
    
    float min_diff = 1e9;
    for (int pair_idx : (use_window_21_81 ? PhSp_21t81_candidateIdx : PhSp_81t101_candidateIdx)) {
        float diff = std::abs(z_candidate_masses[pair_idx] - cfg.z_mass);
        if ((diff < min_diff) && (z_candidates[1-pair_idx].Pt() > (meta.l1flavor == ELECTRON ? cfg.e_pt_min : cfg.mu_pt_min))) {
            // The unselected z-candidate which will be assigned as Higgs prompt lepton must also pass PT cut.
            min_diff = diff;
            // ----------- assign variables in meta -----------
            // pair_idx is the index of the selected pair in the z_candidates vector (0 or 1) and always has 2 members.
            meta.otherZ_idx = (pair_idx == 0) ? meta.l4_index : meta.l3_index;
            meta.otherH_idx = (pair_idx == 0) ? meta.l3_index : meta.l4_index;
            meta.z_l1 = meta.l2_index;
            meta.z_l2 = (pair_idx == 0) ? meta.l3_index : meta.l4_index;
            meta.z_mass = z_candidate_masses[pair_idx];
            meta.z_mass_diff = diff;
            meta.h_mu = (meta.l1flavor == MUON) ? meta.l1_index : ((pair_idx == 0) ? meta.l4_index : meta.l3_index);
            meta.h_e = (meta.l1flavor == ELECTRON) ? meta.l1_index : ((pair_idx == 0) ? meta.l4_index : meta.l3_index);
            meta.h_mu_pt = evt.d->Muon_PT[meta.h_mu];
            meta.h_e_pt = evt.d->Electron_PT[meta.h_e];
            // meta.m_z1 = z_candidate_masses[pair_idx];
            // meta.m_z2 = z_candidate_masses[1-pair_idx];
            // meta.m_h1 = (h1_p4+z_candidates[1-pair_idx]).M();
            // meta.m_h2 = (h1_p4+z_candidates[pair_idx]).M();
            // meta.mh_invaraint = (h1_p4 + z_candidates[1-pair_idx]).M();
             // Post-calculations
            compute_z_post_calculations(evt, meta, cfg);
        }
    }   
    return (min_diff != 1e9);
}

// -----------------------------------------------------------------------------
// H candidate Selections (for mH > 150 GeV)
// -----------------------------------------------------------------------------
std::string HCandidateSelection::name() const { return "HCandidateSelection"; }
bool HCandidateSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    // Choose the pair in which collinear mass is highest -> assign to higgs
    // The pair must > 145 GeV
    int idx_h1 = meta.l1_index; // singleton
    bool is_h1_tau_vis = false;

    int candidates[2] = {meta.l3_index, meta.l4_index};
    double best_mcol = -1.0;
    int best_h2_idx = -1;

    int otherZ_idx = -1;
    int otherH_idx = -1;
    int z2_idx = -1;

    // use 'cfg.mode' to indicate collinear mass selection either mutaue or etaumu
    if (cfg.mode == 1) { // mutaue
        is_h1_tau_vis = (meta.l1flavor == MUON) ? false : true;
    } else if (cfg.mode == 0) { // etamu
        is_h1_tau_vis = (meta.l1flavor == ELECTRON) ? false : true;
    } else {
        std::cout << "HCandidateSelection: Unknown mode '" << cfg.mode << "'. Expected 'mutaue' or 'etamu'." << std::endl;
        return false;
    }
    

    TLorentzVector p4_h1 = get_p4(evt, idx_h1, meta.l1flavor);

    for (int idx_h2 : candidates) {
        TLorentzVector p4_h2 = get_p4(evt, idx_h2, meta.l3flavor);
        // double m_collinear = compute_collinear_mass(evt, tau_vis, other_lep);
        TLorentzVector v_tau_vis = is_h1_tau_vis ? p4_h1 : p4_h2;
        TLorentzVector v_other   = is_h1_tau_vis ? p4_h2 : p4_h1;
        double m_collinear = compute_collinear_mass(evt, v_tau_vis, v_other);
        if (m_collinear > cfg.mcol_min && m_collinear > best_mcol) {
            best_mcol = m_collinear;
            best_h2_idx = idx_h2;
            z2_idx = (idx_h2 == candidates[0]) ? candidates[1] : candidates[0];
            // Alternative pair; this serves as check if our selection is correct
            otherZ_idx = idx_h2;
            otherH_idx = z2_idx;

        }
    }

    // Z mis-identification rejection
    // if the meta.l2_index + best_h2_idx forms a Z in mass window, reject
    // 86-96 GeV
    if (best_h2_idx != -1) {
        TLorentzVector p4_z1 = get_p4(evt, meta.l2_index, meta.l2flavor);
        TLorentzVector p4_z2 = get_p4(evt, best_h2_idx, meta.l2flavor);
        double z_mass = (p4_z1 + p4_z2).M();
        if (z_mass > 86.0 && z_mass < 96.0) return false;
    }   
    
    // Final assignment
    if (best_h2_idx != -1) {
        meta.h_mu = (meta.l1flavor == MUON) ? meta.l1_index : best_h2_idx;
        meta.h_e  = (meta.l1flavor == ELECTRON) ? meta.l1_index : best_h2_idx;
        meta.h_mu_pt = evt.d->Muon_PT[meta.h_mu];
        meta.h_e_pt  = evt.d->Electron_PT[meta.h_e];
        meta.otherZ_idx = otherZ_idx;
        meta.otherH_idx = otherH_idx;
        meta.m_collinear = best_mcol;

        // Z candidate info
        meta.z_l1 = meta.l2_index;
        meta.z_l2 = z2_idx;
        meta.z_flavor = (meta.l1flavor == ELECTRON) ? MUON : ELECTRON;
        meta.z_mass = (get_p4(evt, meta.z_l1, meta.z_flavor) + get_p4(evt, meta.z_l2, meta.z_flavor)).M();
        meta.z_mass_diff = meta.z_mass - cfg.z_mass;

        // Post-calculations
        compute_z_post_calculations(evt, meta, cfg);

        return true;
    }

    return false;
}


// -----------------------------------------------------------------------------
// H -> Mu Tau -> E selection (Mu is prompt, Tau->E)
// -----------------------------------------------------------------------------
std::string HToMuTauESelection::name() const { return "H_to_mutau_e"; }

bool HToMuTauESelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;

    // Require MET > 10 GeV
    if (evt.d->MissingET_MET[0] <= 10.0) return false;

    int zl1 = meta.z_l1;
    int zl2 = meta.z_l2;
    int idx_mu = meta.h_mu;
    int idx_e  = meta.h_e;
    // OS Check
    if (evt.d->Muon_Charge[idx_mu] * evt.d->Electron_Charge[idx_e] >= 0) return false;

    // PT requirements
    if (evt.d->Muon_PT[idx_mu] < cfg.mu_pt_min) return false;
    if (evt.d->Electron_PT[idx_e] < cfg.e_pt_min) return false;
    
    // Compute Kinematics
    TLorentzVector v_mu = get_p4(evt, idx_mu, MUON);
    TLorentzVector v_e  = get_p4(evt, idx_e, ELECTRON);
    
    if (evt.d->MissingET_size > 0) {
        double met_phi = evt.d->MissingET_Phi[0];
        double met     = evt.d->MissingET_MET[0];
        
        meta.dphi_e_met  = delta_phi(v_e.Phi(), met_phi);
        meta.dphi_mu_met = delta_phi(v_mu.Phi(), met_phi);
        meta.m_transverse_e  = compute_transverse_mass(v_e.Pt(), v_e.Phi(), met, met_phi);
        meta.m_transverse_mu = compute_transverse_mass(v_mu.Pt(), v_mu.Phi(), met, met_phi);

        // Collinear Mass: Tau decays to Electron. v_tau_vis = v_e. Other = v_mu.
        meta.m_collinear = compute_collinear_mass(evt, v_e, v_mu);
        meta.m_h1 = meta.m_collinear;
        
        // Invariant Mass (Visible + MET)
        TLorentzVector v_met;
        v_met.SetPtEtaPhiM(met, 0.0, met_phi, 0.0);
        meta.m_h_invariant = (v_mu + v_e + v_met).M();
    }

    meta.dphi_mu_e = std::abs(v_mu.DeltaPhi(v_e));
    meta.deltaR_mu_e = v_mu.DeltaR(v_e);

    // Boosted kinematics
    TVector3 boost_vec(meta.beta_x, meta.beta_y, meta.beta_z);
    v_mu.Boost(boost_vec);
    v_e.Boost(boost_vec);
    meta.h_mu_boosted_pt = v_mu.Pt();
    meta.h_e_boosted_pt = v_e.Pt();

    return true;
}

// -----------------------------------------------------------------------------
// H -> E Tau -> Mu selection (E is prompt, Tau->Mu)
// -----------------------------------------------------------------------------
std::string HToETauMuSelection::name() const { return "H_to_etau_mu"; }

bool HToETauMuSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;

    // Require MET > 10 GeV
    if (evt.d->MissingET_MET[0] <= 10.0) return false;

    int zl1 = meta.z_l1;
    int zl2 = meta.z_l2;
    int idx_tau_mu = meta.h_mu;
    int idx_e      = meta.h_e;
    // OS Check
    if (evt.d->Muon_Charge[idx_tau_mu] * evt.d->Electron_Charge[idx_e] >= 0) return false;
    // PT requirements
    if (evt.d->Muon_PT[idx_tau_mu] < cfg.mu_pt_min) return false;
    if (evt.d->Electron_PT[idx_e] < cfg.e_pt_min) return false;

    TLorentzVector v_mu = get_p4(evt, idx_tau_mu, MUON);
    TLorentzVector v_e  = get_p4(evt, idx_e, ELECTRON);

    if (evt.d->MissingET_size > 0) {
        double met_phi = evt.d->MissingET_Phi[0];
        double met     = evt.d->MissingET_MET[0];

        meta.dphi_e_met  = delta_phi(v_e.Phi(), met_phi);
        meta.dphi_mu_met = delta_phi(v_mu.Phi(), met_phi);
        meta.m_transverse_e  = compute_transverse_mass(v_e.Pt(), v_e.Phi(), met, met_phi);
        meta.m_transverse_mu = compute_transverse_mass(v_mu.Pt(), v_mu.Phi(), met, met_phi);

        // Collinear Mass: Tau decays to Muon. v_tau_vis = v_mu. Other = v_e.
        meta.m_collinear = compute_collinear_mass(evt, v_mu, v_e); 
        meta.m_h1 = meta.m_collinear;
        
        TLorentzVector v_met;
        v_met.SetPtEtaPhiM(met, 0.0, met_phi, 0.0);
        meta.m_h_invariant = (v_mu + v_e + v_met).M();
    }

    meta.dphi_mu_e = std::abs(v_mu.DeltaPhi(v_e));
    meta.deltaR_mu_e = v_mu.DeltaR(v_e);

    TVector3 boost_vec(meta.beta_x, meta.beta_y, meta.beta_z);
    v_mu.Boost(boost_vec);
    v_e.Boost(boost_vec);
    meta.h_mu_boosted_pt = v_mu.Pt();
    meta.h_e_boosted_pt = v_e.Pt();

    return true;
}

// -----------------------------------------------------------------------------
// MET Cuts
// -----------------------------------------------------------------------------
std::string METEDphiSelection::name() const { return "MET_E_dphi"; }
bool METEDphiSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d || meta.h_e < 0 || evt.d->MissingET_size <= 0) return false;
    if (std::isnan(meta.dphi_e_met)) return false;
    return meta.dphi_e_met < cfg.max_dphi_e_met;
}

std::string METMuDphiSelection::name() const { return "MET_Mu_dphi"; }
bool METMuDphiSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d || meta.h_mu < 0 || evt.d->MissingET_size <= 0) return false;
    if (std::isnan(meta.dphi_mu_met)) return false;
    return meta.dphi_mu_met < cfg.max_dphi_mu_met;
}

// -----------------------------------------------------------------------------
// Recoil Mass Selection
// -----------------------------------------------------------------------------
std::string RecoilMassSelection::name() const { return "RecoilMassSelection"; }
bool RecoilMassSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;

    // Helper to calc recoil for a pair
    auto get_recoil = [&](int i1, int flav1, int i2, int flav2) {
        TLorentzVector v1 = get_p4(evt, i1, flav1);
        TLorentzVector v2 = get_p4(evt, i2, flav2);
        TLorentzVector beam(0, 0, 0, ECM);
        return (beam - (v1 + v2)).M();
    };

    // Calculate recoil for both possible Z pairs within the triplet
    // (Recall: l2 is SS, l3/l4 are OS relative to l1)
    double recoil_23 = get_recoil(meta.l2_index, meta.l2flavor, meta.l3_index, meta.l3flavor);
    double recoil_24 = get_recoil(meta.l2_index, meta.l2flavor, meta.l4_index, meta.l4flavor);

    meta.m_recoil1 = std::max(recoil_23, recoil_24);
    meta.m_recoil2 = std::min(recoil_23, recoil_24);

    if (meta.m_recoil1 > cfg.recoil_mass_min) {
        // Optional: Update Z candidate to be the one with best recoil?
        // Current logic in ZCandidateSelection uses Z mass window.
        // If recoil selection is preferred, logic goes here.
        return true;
    }
    return false;
}

// Unused function -----------------------------------
// Old Z to ll selection
std::string ZToLLSelection::name() const { return "Z_to_ll"; }
bool ZToLLSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    return true;
}


} // namespace hlfv
