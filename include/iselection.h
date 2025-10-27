#pragma once

#include <string>
#include "types.h"

namespace hlfv {

struct ISelection {
    virtual ~ISelection() {}
    virtual std::string name() const = 0;
    virtual bool apply(const Event& evt, Meta& meta, const Parameters& cfg) = 0;
};

} // namespace hlfv
