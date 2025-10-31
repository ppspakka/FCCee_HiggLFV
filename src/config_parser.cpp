#include "../include/config_parser.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

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

    PipelineConfig cfg; // defaults already set by struct defaults

    // Parse parameters object
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
    }

    // Parse selections array
    std::string arr = extractBlock(text, "\"selections\"", '[', ']');
    if (!arr.empty()) {
        // Iterate objects { ... } in arr
        size_t i = 0; while (i < arr.size()) {
            if (arr[i] == '{') {
                int depth = 0; size_t start = i;
                while (i < arr.size()) {
                    if (arr[i]=='{') depth++; else if (arr[i]=='}') depth--; i++;
                    if (depth==0) break;
                }
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
            } else {
                i++;
            }
        }
    }

    if (cfg.selections.empty()) {
        cfg.selections.push_back({"Z_to_ll", true});
        cfg.selections.push_back({"H_to_mue", true});
        cfg.selections.push_back({"MET_dphi", true});
    }

    out = cfg;
    return true;
}
