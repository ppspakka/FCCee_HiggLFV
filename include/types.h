#pragma once

#include <limits>

class Delphes; // forward declaration

namespace hlfv {

struct Parameters {
    // Lepton selection
    double lepton_pt_min;
    // Z->ll parameters
    double z_mass;
    double zl_pt_min;
    double z_mass_window_upper;
    double z_mass_window_lower;
    // H->mue parameters
    double mu_pt_min;
    double e_pt_min;
    // MET selection
    double max_dphi_e_met;
    double max_dphi_mu_met;

    // Recoil mass selection
    double recoil_mass_min;

    // H candidate selection mode: "mutaue" or "etamu"
    double mode;
    double mcol_min;
};

struct Event {
    Delphes* d = nullptr; // bound Delphes tree object
};

struct Meta {
    // Lepton selection info (unset -> NaN / -1)
    int l1flavor = -1; // 0=e,1=mu
    int l2flavor = -1;
    int l3flavor = -1;
    int l4flavor = -1;
    int l1_index = -1;
    int l2_index = -1;
    int l3_index = -1;
    int l4_index = -1;

    int otherZ_idx = -1;
    int otherH_idx = -1;

    // Z candidate info (unset -> NaN / -1)
    int z_l1 = -1;
    int z_l2 = -1;
    int z_flavor = -1; // 0=e, 1=mu
    double z_mass = std::numeric_limits<double>::quiet_NaN();
    double z_mass_diff = std::numeric_limits<double>::quiet_NaN();
    // H->mu tau_e or H->e tau_mu candidate info
    int h_mu = -1;
    int h_e = -1;
    double h_mu_pt = std::numeric_limits<double>::quiet_NaN();
    double h_mu_boosted_pt = std::numeric_limits<double>::quiet_NaN();
    double h_e_pt = std::numeric_limits<double>::quiet_NaN();
    double h_e_boosted_pt = std::numeric_limits<double>::quiet_NaN();
    double dphi_e_met = std::numeric_limits<double>::quiet_NaN();
    double dphi_mu_met = std::numeric_limits<double>::quiet_NaN();
    double dphi_mu_e = std::numeric_limits<double>::quiet_NaN();
    // collinear mass
    double m_collinear = std::numeric_limits<double>::quiet_NaN();
    // transverse mass
    double m_transverse_e = std::numeric_limits<double>::quiet_NaN();
    double m_transverse_mu = std::numeric_limits<double>::quiet_NaN();
    // recoil mass
    double m_recoil = std::numeric_limits<double>::quiet_NaN();
    double m_recoil1 = std::numeric_limits<double>::quiet_NaN();
    double m_recoil2 = std::numeric_limits<double>::quiet_NaN();

    double deltaR_mu_e = std::numeric_limits<double>::quiet_NaN();

    // m_Z candidate case 1 (closest pair)
    double m_z1 = std::numeric_limits<double>::quiet_NaN();
    // m_Z candidate case 2 (second pair)
    double m_z2 = std::numeric_limits<double>::quiet_NaN();
    // m_h case 1 (choose Z closest pair)
    double m_h1 = std::numeric_limits<double>::quiet_NaN();
    // m_h case 2 choose Z second pair, remaining used to calculate col mass
    double m_h2 = std::numeric_limits<double>::quiet_NaN();

    // displaced track: abs(D0), abs(DZ) (only for higgs daughters)
    double h_e_d0 = std::numeric_limits<double>::quiet_NaN();
    double h_e_dz = std::numeric_limits<double>::quiet_NaN();
    double h_mu_d0 = std::numeric_limits<double>::quiet_NaN();
    double h_mu_dz = std::numeric_limits<double>::quiet_NaN();

    // boost vector
    double beta_x = std::numeric_limits<double>::quiet_NaN();
    double beta_y = std::numeric_limits<double>::quiet_NaN();
    double beta_z = std::numeric_limits<double>::quiet_NaN();

    // m_h reconstructed using invariant mass of 3 objects (2l from H + MET)
    double m_h_invariant = std::numeric_limits<double>::quiet_NaN();

    double MET = std::numeric_limits<double>::quiet_NaN();
};

} // namespace hlfv
