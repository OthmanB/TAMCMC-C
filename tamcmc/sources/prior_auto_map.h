#pragma once
#include <string>
#include <vector>

// Curated mapping from model_fullname to prior_fct_name.
// Source of truth for MS_Global variants: TAMCMC-C/Config/default/models_ctrl.list
// Only non-obsolete models (cases 6,7,8,12,13,18,19,21,23) are included.

// Sentinel value recognized by is_auto_sentinel()
extern const std::string PRIOR_AUTO_SENTINEL;

// Returns the canonical prior_fct_name for a given model_fullname.
// Returns empty string "" if model_fullname is not in the curated table.
// Caller is responsible for handling the unknown case (fatal exit).
std::string resolve_prior_from_model(const std::string& model_fullname);

// Returns true iff prior_fct_name is the auto sentinel after:
//   1. Trimming leading/trailing whitespace
//   2. Stripping a trailing ';' if present
//   3. Exact case-sensitive comparison with "auto"
// Returns false for "AUTO", "Auto", "", "auto;extra", etc.
bool is_auto_sentinel(const std::string& prior_fct_name);

// Returns the complete list of model_fullname strings in the curated table.
// Useful for error messages listing supported models.
std::vector<std::string> supported_model_fullnames();
