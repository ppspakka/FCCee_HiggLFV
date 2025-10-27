#pragma once

#include <memory>
#include <string>
#include "iselection.h"

// Factory to create selections by name
std::unique_ptr<ISelection> makeSelectionByName(const std::string& name);
