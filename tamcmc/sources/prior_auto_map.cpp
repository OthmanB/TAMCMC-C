#include "prior_auto_map.h"
#include <map>
#include <algorithm>

// Sentinel value for auto-detection
const std::string PRIOR_AUTO_SENTINEL = "auto";

// Curated mapping table: model_fullname -> prior_fct_name
// 15 entries covering non-obsolete models (cases 6,7,8,12,13,18,19,21,23)
static const std::map<std::string, std::string> prior_auto_map = {
    {"model_Harvey_Gaussian", "priors_Harvey_Gaussian"},
    {"model_Kallinger2014_Gaussian", "priors_Kallinger2014_Gaussian"},
    {"model_MS_Global_a1l_etaa3_HarveyLike", "io_MS_Global"},
    {"model_MS_Global_a1n_etaa3_HarveyLike", "io_MS_Global"},
    {"model_MS_Global_a1nl_etaa3_HarveyLike", "io_MS_Global"},
    {"model_MS_Global_a1etaa3_HarveyLike_Classic_v2", "io_MS_Global"},
    {"model_MS_Global_a1etaa3_HarveyLike_Classic_v3", "io_MS_Global"},
    {"model_MS_Global_a1n_a2a3_HarveyLike", "io_MS_Global"},
    {"model_MS_Global_a1nl_a2a3_HarveyLike", "io_MS_Global"},
    {"model_MS_Global_ajAlm_HarveyLike", "io_MS_Global"},
    {"model_MS_Global_aj_HarveyLike", "io_MS_Global"},
    {"model_RGB_asympt_aj_AppWidth_HarveyLike_v4", "io_asymptotic"},
    {"model_RGB_asympt_aj_CteWidth_HarveyLike_v4", "io_asymptotic"},
    {"model_MS_local_basic", "io_local"},
    {"model_MS_local_Hnlm", "io_local"}
};

std::string resolve_prior_from_model(const std::string& model_fullname) {
    auto it = prior_auto_map.find(model_fullname);
    if (it != prior_auto_map.end()) {
        return it->second;
    }
    return "";
}

bool is_auto_sentinel(const std::string& prior_fct_name) {
    // Trim leading and trailing whitespace (spaces and tabs)
    size_t start = prior_fct_name.find_first_not_of(" \t");
    if (start == std::string::npos) {
        return false; // Empty or all whitespace
    }
    size_t end = prior_fct_name.find_last_not_of(" \t");
    std::string trimmed = prior_fct_name.substr(start, end - start + 1);
    
    // Strip one trailing semicolon if present
    if (!trimmed.empty() && trimmed.back() == ';') {
        trimmed.pop_back();
    }
    
    // Exact case-sensitive comparison with "auto"
    return trimmed == PRIOR_AUTO_SENTINEL;
}

std::vector<std::string> supported_model_fullnames() {
    std::vector<std::string> result;
    for (const auto& entry : prior_auto_map) {
        result.push_back(entry.first);
    }
    return result;
}
