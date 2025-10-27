#pragma once

#include <limits>

class Delphes; // forward declaration

namespace hlfv {

struct Parameters {
    double z_mass = 91.1876;     // GeV
    double z_mass_window = 10.0; // |Mll - MZ| < window (GeV)
    // H->mue thresholds
    double mu_pt_min = 30.0;     // GeV
    double e_pt_min = 20.0;      // GeV
    // MET selection
    double max_dphi_e_met = 0.7; // selection threshold
};

struct Event {
    Delphes* d = nullptr; // bound Delphes tree object
};

struct Meta {
    // Z candidate info (unset -> NaN / -1)
    int z_l1 = -1;
    int z_l2 = -1;
    int z_flavor = -1; // 0=e, 1=mu
    double z_mass = std::numeric_limits<double>::quiet_NaN();
    double z_mass_diff = std::numeric_limits<double>::quiet_NaN();
    // H->mu e candidate (not from Z)
    int h_mu = -1;
    int h_e = -1;
    double dphi_e_met = std::numeric_limits<double>::quiet_NaN();
};

} // namespace hlfv
