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
    double z_mass_window;
    // H->mue parameters
    double mu_pt_min;
    double e_pt_min;
    // MET selection
    double max_dphi_e_met;
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

    // Z candidate info (unset -> NaN / -1)
    int z_l1 = -1;
    int z_l2 = -1;
    int z_flavor = -1; // 0=e, 1=mu
    double z_mass = std::numeric_limits<double>::quiet_NaN();
    double z_mass_diff = std::numeric_limits<double>::quiet_NaN();
    // H->mu e candidate (not from Z)
    int h_mu = -1;
    int h_e = -1;
    double h_mu_pt = std::numeric_limits<double>::quiet_NaN();
    double h_e_pt = std::numeric_limits<double>::quiet_NaN();
    double dphi_e_met = std::numeric_limits<double>::quiet_NaN();
    double dphi_mu_e = std::numeric_limits<double>::quiet_NaN();
    // collinear mass
    double m_collinear = std::numeric_limits<double>::quiet_NaN();
};

} // namespace hlfv
