#include "../include/pipeline_config.h"
#include "../include/selections.h"
#include "../include/histman.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

namespace hlfv {

// ---- selection factory ----
std::unique_ptr<ISelection> makeSelectionByName(const std::string& name) {
    std::string key = name;
    std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c){ return std::tolower(c); });

    if (key == "empty_selection") return std::make_unique<EmptySelection>();

    if (key == "finalstate_nocut" || key == "final_state_nocut" || key == "finalstate-nocut") return std::make_unique<FinalState_NoCut>();

    if (key == "lepton_selection" || key == "leptonselection") return std::make_unique<LeptonSelection>();

    if (key == "z_candidate_selection" || key == "zcandidate_selection") return std::make_unique<ZCandidateSelection>();

    if (key == "z_to_ll" || key == "z->ll" || key == "z_ll") return std::make_unique<ZToLLSelection>();

    // H to mu tau_e
    if (key == "h_to_mutau_e" || key == "h_to_mutaue" || key == "h_to_mue" || key == "h->mutau_e") return std::make_unique<HToMuTauESelection>();

    if (key == "h_to_etau_mu") return std::make_unique<HToETauMuSelection>();
    
    // if (key == "met_dphi" || key == "met-dphi" || key == "metdphi") return std::make_unique<METDphiSelection>();
    if (key == "met_e_dphi" || key == "metedphi" || key == "met_e-dphi") return std::make_unique<METEDphiSelection>();
    if (key == "met_mu_dphi" || key == "metmudphi" || key == "met_mu-dphi") return std::make_unique<METMuDphiSelection>();

    if (key == "recoil_mass_selection" || key == "recoilmass_selection" || key == "recoil-mass-selection") return std::make_unique<RecoilMassSelection>();

    if (key == "h_candidate_selection" || key == "hcandidate_selection" || key == "h-candidate-selection") return std::make_unique<HCandidateSelection>();
    return nullptr;
}

// ---- minimal config parser (schema-specific) ----
namespace {
static inline void ltrim(std::string &s) { s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch){return !std::isspace(ch);})); }
static inline void rtrim(std::string &s) { s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch){return !std::isspace(ch);}).base(), s.end()); }
static inline void trim(std::string &s) { ltrim(s); rtrim(s); }

static bool parseNumber(const std::string& v, double& out) {
    try { size_t idx=0; out = std::stod(v, &idx); return idx==v.size(); } catch(...) { return false; }
}

static std::string extractBlock(const std::string& text, const std::string& key, char open, char close) {
    size_t k = text.find(key);
    if (k == std::string::npos) return {};
    size_t lb = text.find(open, k);
    if (lb == std::string::npos) return {};
    int depth = 0; for (size_t i = lb; i < text.size(); ++i) {
        if (text[i]==open) depth++;
        else if (text[i]==close) depth--;
        if (depth==0) { return text.substr(lb+1, i - lb - 1); }
    }
    return {};
}
}

bool loadPipelineConfig(const std::string& filepath, PipelineConfig& out) {
    std::ifstream in(filepath);
    if (!in.is_open()) return false;
    std::stringstream buf; buf << in.rdbuf();
    std::string text = buf.str();

    PipelineConfig cfg; // defaults

    // parameters
    std::string p = extractBlock(text, "\"parameters\"", '{', '}');
    if (!p.empty()) {
        auto parseParam = [&](const char* key, double& dst){
            std::string k = std::string("\"") + key + "\"";
            size_t pos = p.find(k);
            if (pos != std::string::npos) {
                size_t colon = p.find(':', pos);
                if (colon != std::string::npos) {
                    size_t end = p.find_first_of(",}\n\r", colon+1);
                    std::string val = p.substr(colon+1, end - (colon+1));
                    trim(val);
                    parseNumber(val, dst);
                }
            }
        };
        parseParam("lepton_pt_min", cfg.params.lepton_pt_min);
        parseParam("z_mass", cfg.params.z_mass);
        parseParam("zl_pt_min", cfg.params.zl_pt_min);
        parseParam("z_mass_window_upper", cfg.params.z_mass_window_upper);
        parseParam("z_mass_window_lower", cfg.params.z_mass_window_lower);
        parseParam("mu_pt_min", cfg.params.mu_pt_min);
        parseParam("e_pt_min", cfg.params.e_pt_min);
        parseParam("max_dphi_e_met", cfg.params.max_dphi_e_met);
        parseParam("max_dphi_mu_met", cfg.params.max_dphi_mu_met);
        parseParam("mode", cfg.params.mode);
        parseParam("mcol_min", cfg.params.mcol_min);

    }

    // selections array
    std::string arr = extractBlock(text, "\"selections\"", '[', ']');
    if (!arr.empty()) {
        size_t i = 0; while (i < arr.size()) {
            if (arr[i] == '{') {
                int depth = 0; size_t start = i;
                while (i < arr.size()) { if (arr[i]=='{') depth++; else if (arr[i]=='}') depth--; i++; if (depth==0) break; }
                std::string obj = arr.substr(start, i - start);
                SelectionSpec spec; spec.enabled = true;
                // name
                size_t nk = obj.find("\"name\"");
                if (nk != std::string::npos) {
                    size_t colon = obj.find(':', nk);
                    size_t q1 = obj.find('"', colon+1);
                    size_t q2 = obj.find('"', q1+1);
                    if (colon!=std::string::npos && q1!=std::string::npos && q2!=std::string::npos) {
                        spec.name = obj.substr(q1+1, q2 - q1 - 1);
                    }
                }
                // enabled
                size_t ek = obj.find("\"enabled\"");
                if (ek != std::string::npos) {
                    size_t colon = obj.find(':', ek);
                    if (colon != std::string::npos) {
                        size_t end = obj.find_first_of(",}\n\r", colon+1);
                        std::string val = obj.substr(colon+1, end - (colon+1));
                        trim(val);
                        if (val.find("true") != std::string::npos) spec.enabled = true;
                        else if (val.find("false") != std::string::npos) spec.enabled = false;
                    }
                }
                if (!spec.name.empty()) cfg.selections.push_back(spec);
            } else { i++; }
        }
    }

    if (cfg.selections.empty()) {
        // error: at least one selection
        return false;
    }

    out = cfg;
    return true;
}

// ---- variables registry ----
std::vector<HistogramManager::VarSpec> Variables::getDefault() {
    std::vector<HistogramManager::VarSpec> vars;
    vars.push_back({
        "n_muons", "Number of muons;N_{#mu};Events", 6, -0.5, 5.5,
        [](const Event& evt, const Meta&) -> double { return evt.d ? (double)evt.d->Muon_size : std::numeric_limits<double>::quiet_NaN(); }
    });
    vars.push_back({
        "z_mass", "M_{ll} (GeV);M_{ll} [GeV];Events", 60, 60.0, 120.0,
        [](const Event&, const Meta& m) -> double { return m.z_mass; }
    });
    vars.push_back({
        "z_mass_diff", "|M_{ll}-M_{Z}| (GeV);|M_{ll}-M_{Z}| [GeV];Events", 50, -25.0, 25.0,
        [](const Event&, const Meta& m) -> double { return m.z_mass_diff; }
    });
    vars.push_back({
        "h_mu_pt", "p_{T}(#mu) (GeV);p_{T}(#mu) [GeV];Events", 50, 0.0, 100.0,
        [](const Event&, const Meta& m) -> double { return m.h_mu_pt; }
    });
    vars.push_back({
        "h_e_pt", "p_{T}(e) (GeV);p_{T}(e) [GeV];Events", 50, 0.0, 100.0,
        [](const Event&, const Meta& m) -> double { return m.h_e_pt; }
    });
    vars.push_back({
        "dphi_e_met", "|#Delta#phi(e,MET)|;|#Delta#phi|;Events", 64, 0.0, 3.2,
        [](const Event&, const Meta& m) -> double { return m.dphi_e_met; }
    });
    vars.push_back({
        "dphi_mu_met", "|#Delta#phi(#mu,MET)|;|#Delta#phi|;Events", 64, 0.0, 3.2,
        [](const Event&, const Meta& m) -> double { return m.dphi_mu_met; }
    });
    vars.push_back({
        "dphi_mu_e", "|#Delta#phi(#mu,e)|;|#Delta#phi|;Events", 64, 0.0, 3.2,
        [](const Event&, const Meta& m) -> double { return m.dphi_mu_e; }
    });
    vars.push_back({
        "m_collinear", "Collinear mass (GeV);M_{collinear} [GeV];Events", 200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars.push_back({
        "m_transverse_e_fine", "Transverse mass (e) (GeV);M_{T}(e) [GeV];Events", 300, 0.01, 3.01,
        [](const Event&, const Meta& m) -> double { return m.m_transverse_e; }
    });
    vars.push_back({
        "m_transverse_mu_fine", "Transverse mass (#mu) (GeV);M_{T}(#mu) [GeV];Events", 300, 0.01, 3.01,
        [](const Event&, const Meta& m) -> double { return m.m_transverse_mu; }
    });

    // rough transverse mass histograms
    vars.push_back({
        "m_transverse_e", "Transverse mass (e) (GeV);M_{T}(e) [GeV];Events", 120, 0.0, 120.0,
        [](const Event&, const Meta& m) -> double { return m.m_transverse_e; }
    });
    vars.push_back({
        "m_transverse_mu", "Transverse mass (#mu) (GeV);M_{T}(#mu) [GeV];Events", 120, 0.0, 120.0,
        [](const Event&, const Meta& m) -> double { return m.m_transverse_mu; }
    });
    vars.push_back({
        "deltaR_mu_e", "#DeltaR(#mu,e);#DeltaR;Events", 60, 0.0, 6.0,
        [](const Event&, const Meta& m) -> double { return m.deltaR_mu_e; }
    });
    vars.push_back({
        "m_recoil", "Recoil mass (GeV);M_{recoil} [GeV];Events", 200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_recoil; }
    });

    // additional variables for Z and H candidates
    vars.push_back({
        "m_z1", "M_{Z1} (GeV);M_{Z1} [GeV];Events", 200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_z1; }
    });
    vars.push_back({
        "m_z2", "M_{Z2} (GeV);M_{Z2} [GeV];Events", 200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_z2; }
    });
    vars.push_back({
        "m_h1", "M_{H1} (GeV);M_{H1} [GeV];Events", 200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_h1; }
    });
    vars.push_back({
        "m_h2", "M_{H2} (GeV);M_{H2} [GeV];Events", 200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_h2; }
    });
    vars.push_back({
        "h_mu_boosted_pt", "Boosted p_{T}(#mu) (GeV);p_{T}(#mu) [GeV];Events", 50, 0.0, 100.0,
        [](const Event&, const Meta& m) -> double { return m.h_mu_boosted_pt; }
    });
    vars.push_back({
        "h_e_boosted_pt", "Boosted p_{T}(e) (GeV);p_{T}(e) [GeV];Events", 50, 0.0, 100.0,
        [](const Event&, const Meta& m) -> double { return m.h_e_boosted_pt; }
    });

    vars.push_back({
        "h_e_d0", "|D0(e)| (Unknow unit);|D0(e)| [Unknow Unit];Events", 100, 0.0, 0.1,
        [](const Event&, const Meta& m) -> double { return m.h_e_d0; }
    });

    vars.push_back({
        "h_e_dz", "|DZ(e)| (Unknow unit);|DZ(e)| [Unknow Unit];Events", 100, 0.0, 0.5,
        [](const Event&, const Meta& m) -> double { return m.h_e_dz; }
    });

    vars.push_back({
        "h_mu_d0", "|D0(#mu)| (Unknow unit);|D0(#mu)| [Unknow Unit];Events", 100, 0.0, 0.1,
        [](const Event&, const Meta& m) -> double { return m.h_mu_d0; }
    });

    vars.push_back({
        "h_mu_dz", "|DZ(#mu)| (Unknow unit);|DZ(#mu)| [Unknow Unit];Events", 100, 0.0, 0.5,
        [](const Event&, const Meta& m) -> double { return m.h_mu_dz; }
    });

    // Recoil mass related variables
    // highest recoil mass pair
    vars.push_back({
        "recoil_mass_1", "Recoil mass 1 (GeV);M_{recoil 1} [GeV];Events", 200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_recoil1; }
    });
    // second highest recoil mass pair
    vars.push_back({
        "recoil_mass_2", "Recoil mass 2 (GeV);M_{recoil 2} [GeV];Events", 200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_recoil2; }
    });

    // invariant mass of 3 objects (2l from H + MET)
    vars.push_back({
        "m_h_invariant", "Invariant mass of 3 objects (GeV);M_{H inv} [GeV];Events", 200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_h_invariant; }
    });

    return vars;
}

// ---- 2D variables registry ----
std::vector<HistogramManager::Var2DSpec> Variables::getDefault2D() {
    std::vector<HistogramManager::Var2DSpec> vars2d;
    // Example correlations using existing Meta fields
    vars2d.push_back({
        "mZ_vs_mH1", "M_{Z1} vs M_{H1};M_{Z1} [GeV];M_{H1} [GeV]",
        100, 0.5, 200.5,
        100, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_z1; },
        [](const Event&, const Meta& m) -> double { return m.m_h1; }
    });
    vars2d.push_back({
        "mZ2_vs_mH2", "M_{Z2} vs M_{H2};M_{Z2} [GeV];M_{H2} [GeV]",
        100, 0.5, 200.5,
        100, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_z2; },
        [](const Event&, const Meta& m) -> double { return m.m_h2; }
    });
    vars2d.push_back({
        "pt_mu_vs_pt_e", "p_{T}(#mu) vs p_{T}(e);p_{T}(#mu) [GeV];p_{T}(e) [GeV]",
        50, 0.0, 100.0,
        50, 0.0, 100.0,
        [](const Event&, const Meta& m) -> double { return m.h_mu_pt; },
        [](const Event&, const Meta& m) -> double { return m.h_e_pt; }
    });

    // Add all 1D hist + col mass
    vars2d.push_back({
        "h_mu_pt_vs_m_collinear", "p_{T}(#mu) vs M_{collinear};p_{T}(#mu) [GeV];M_{collinear} [GeV]",
        50, 0.0, 100.0,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.h_mu_pt; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "h_e_pt_vs_m_collinear", "p_{T}(e) vs M_{collinear};p_{T}(e) [GeV];M_{collinear} [GeV]",
        50, 0.0, 100.0,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.h_e_pt; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "dphi_e_met_vs_m_collinear", "|#Delta#phi(e,MET)| vs M_{collinear};|#Delta#phi(e,MET)|;M_{collinear} [GeV]",
        64, 0.0, 3.2,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.dphi_e_met; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "dphi_mu_met_vs_m_collinear", "|#Delta#phi(#mu,MET)| vs M_{collinear};|#Delta#phi(#mu,MET)|;M_{collinear} [GeV]",
        64, 0.0, 3.2,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.dphi_mu_met; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "dphi_mu_e_vs_m_collinear", "|#Delta#phi(#mu,e)| vs M_{collinear};|#Delta#phi(#mu,e)|;M_{collinear} [GeV]",
        64, 0.0, 3.2,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.dphi_mu_e; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "m_transverse_e_vs_m_collinear", "M_{T}(e) vs M_{collinear};M_{T}(e) [GeV];M_{collinear} [GeV]",
        120, 0.0, 120.0,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_transverse_e; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "m_transverse_mu_vs_m_collinear", "M_{T}(#mu) vs M_{collinear};M_{T}(#mu) [GeV];M_{collinear} [GeV]",
        120, 0.0, 120.0,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_transverse_mu; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "m_transverse_e_fine_vs_m_collinear", "M_{T}(e) fine vs M_{collinear};M_{T}(e) [GeV];M_{collinear} [GeV]",
        300, 0.01, 3.01,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_transverse_e; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "m_transverse_mu_fine_vs_m_collinear", "M_{T}(#mu) fine vs M_{collinear};M_{T}(#mu) [GeV];M_{collinear} [GeV]",
        300, 0.01, 3.01,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.m_transverse_mu; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "deltaR_mu_e_vs_m_collinear", "#DeltaR(#mu,e) vs M_{collinear};#DeltaR(#mu,e);M_{collinear} [GeV]",
        60, 0.0, 6.0,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.deltaR_mu_e; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "h_mu_boosted_pt_vs_m_collinear", "Boosted p_{T}(#mu) vs M_{collinear};p_{T}(#mu) [GeV];M_{collinear} [GeV]",
        50, 0.0, 100.0,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.h_mu_boosted_pt; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "h_e_boosted_pt_vs_m_collinear", "Boosted p_{T}(e) vs M_{collinear};p_{T}(e) [GeV];M_{collinear} [GeV]",
        50, 0.0, 100.0,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.h_e_boosted_pt; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    // displaced track variables
    vars2d.push_back({
        "h_e_d0_vs_m_collinear", "|D0(e)| vs M_{collinear};|D0(e)| [Unknow Unit];M_{collinear} [GeV]",
        100, 0.0, 0.1,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.h_e_d0; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "h_e_dz_vs_m_collinear", "|DZ(e)| vs M_{collinear};|DZ(e)| [Unknow Unit];M_{collinear} [GeV]",
        100, 0.0, 0.5,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.h_e_dz; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "h_mu_d0_vs_m_collinear", "|D0(#mu)| vs M_{collinear};|D0(#mu)| [Unknow Unit];M_{collinear} [GeV]",
        100, 0.0, 0.1,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.h_mu_d0; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    vars2d.push_back({
        "h_mu_dz_vs_m_collinear", "|DZ(#mu)| vs M_{collinear};|DZ(#mu)| [Unknow Unit];M_{collinear} [GeV]",
        100, 0.0, 0.5,
        200, 0.5, 200.5,
        [](const Event&, const Meta& m) -> double { return m.h_mu_dz; },
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });


    return vars2d;
}

} // namespace hlfv
