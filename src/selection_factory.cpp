#include "../include/selection_factory.h"
#include "../include/selections.h"

#include <algorithm>

std::unique_ptr<ISelection> makeSelectionByName(const std::string& name) {
    std::string key = name;
    // Normalize a bit: allow case-insensitive and hyphen/space variations if needed
    std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c){ return std::tolower(c); });

    if (key == "z_to_ll" || key == "z->ll" || key == "z_ll") {
        return std::make_unique<ZToLLSelection>();
    } else if (key == "h_to_mue" || key == "h->mue" || key == "hmue") {
        return std::make_unique<HToMuESelection>();
    } else if (key == "met_dphi" || key == "met-dphi" || key == "metdphi") {
        return std::make_unique<METDphiSelection>();
    } else if (key == "empty_selection") {
        return std::make_unique<EmptySelection>();
    }
    return nullptr;
}
