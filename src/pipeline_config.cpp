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
    if (key == "z_to_ll" || key == "z->ll" || key == "z_ll") return std::make_unique<ZToLLSelection>();
    if (key == "h_to_mue" || key == "h->mue" || key == "hmue") return std::make_unique<HToMuESelection>();
    if (key == "met_dphi" || key == "met-dphi" || key == "metdphi") return std::make_unique<METDphiSelection>();
    if (key == "empty_selection") return std::make_unique<EmptySelection>();
    if (key == "finalstate_nocut" || key == "final_state_nocut" || key == "finalstate-nocut") return std::make_unique<FinalState_NoCut>();
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
        parseParam("z_mass", cfg.params.z_mass);
        parseParam("zl_pt_min", cfg.params.zl_pt_min);
        parseParam("z_mass_window", cfg.params.z_mass_window);
        parseParam("mu_pt_min", cfg.params.mu_pt_min);
        parseParam("e_pt_min", cfg.params.e_pt_min);
        parseParam("max_dphi_e_met", cfg.params.max_dphi_e_met);
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
        "z_mass_diff", "|M_{ll}-M_{Z}| (GeV);|M_{ll}-M_{Z}| [GeV];Events", 50, 0.0, 25.0,
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
        "dphi_mu_e", "|#Delta#phi(#mu,e)|;|#Delta#phi|;Events", 64, 0.0, 3.2,
        [](const Event&, const Meta& m) -> double { return m.dphi_mu_e; }
    });
    vars.push_back({
        "m_collinear", "Collinear mass (GeV);M_{collinear} [GeV];Events", 100, 0.0, 200.0,
        [](const Event&, const Meta& m) -> double { return m.m_collinear; }
    });
    return vars;
}

} // namespace hlfv
