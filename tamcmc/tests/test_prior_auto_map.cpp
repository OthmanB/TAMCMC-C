// Standalone unit tests for prior_auto_map.{h,cpp}
// No external test framework — compile and run directly.
// Exit 0 if all 26 tests pass, exit 1 otherwise.

#include <iostream>
#include <string>
#include <vector>
#include "../sources/prior_auto_map.h"

static int total_tests = 0;
static int failed_tests = 0;

static void check_bool(int test_num, const std::string& description,
                       bool got, bool expected) {
    ++total_tests;
    if (got == expected) {
        std::cout << "PASS: " << description << "\n";
    } else {
        ++failed_tests;
        std::cout << "FAIL: test #" << test_num << " " << description
                  << " expected " << (expected ? "true" : "false")
                  << ", got " << (got ? "true" : "false") << "\n";
    }
}

static void check_str(int test_num, const std::string& description,
                      const std::string& got, const std::string& expected) {
    ++total_tests;
    if (got == expected) {
        std::cout << "PASS: " << description << "\n";
    } else {
        ++failed_tests;
        std::cout << "FAIL: test #" << test_num << " " << description
                  << " expected \"" << expected << "\", got \"" << got << "\"\n";
    }
}

static void check_size(int test_num, const std::string& description,
                       std::size_t got, std::size_t expected) {
    ++total_tests;
    if (got == expected) {
        std::cout << "PASS: " << description << "\n";
    } else {
        ++failed_tests;
        std::cout << "FAIL: test #" << test_num << " " << description
                  << " expected " << expected << ", got " << got << "\n";
    }
}

int main() {
    // ---- is_auto_sentinel ----
    check_bool(1,  "is_auto_sentinel(\"auto\") == true",
               is_auto_sentinel("auto"), true);

    check_bool(2,  "is_auto_sentinel(\"auto;\") == true",
               is_auto_sentinel("auto;"), true);

    check_bool(3,  "is_auto_sentinel(\" auto \") == true (whitespace trimmed)",
               is_auto_sentinel(" auto "), true);

    check_bool(4,  "is_auto_sentinel(\"AUTO\") == false (case-sensitive)",
               is_auto_sentinel("AUTO"), false);

    check_bool(5,  "is_auto_sentinel(\"Auto\") == false",
               is_auto_sentinel("Auto"), false);

    check_bool(6,  "is_auto_sentinel(\"\") == false",
               is_auto_sentinel(""), false);

    check_bool(7,  "is_auto_sentinel(\"io_MS_Global\") == false",
               is_auto_sentinel("io_MS_Global"), false);

    // ---- resolve_prior_from_model ----
    check_str(8,  "resolve_prior_from_model(\"model_Harvey_Gaussian\") == \"priors_Harvey_Gaussian\"",
              resolve_prior_from_model("model_Harvey_Gaussian"),
              "priors_Harvey_Gaussian");

    check_str(9,  "resolve_prior_from_model(\"model_Kallinger2014_Gaussian\") == \"priors_Kallinger2014_Gaussian\"",
              resolve_prior_from_model("model_Kallinger2014_Gaussian"),
              "priors_Kallinger2014_Gaussian");

    check_str(10, "resolve_prior_from_model(\"model_MS_Global_a1l_etaa3_HarveyLike\") == \"io_MS_Global\"",
              resolve_prior_from_model("model_MS_Global_a1l_etaa3_HarveyLike"),
              "io_MS_Global");

    check_str(11, "resolve_prior_from_model(\"model_MS_Global_a1n_etaa3_HarveyLike\") == \"io_MS_Global\"",
              resolve_prior_from_model("model_MS_Global_a1n_etaa3_HarveyLike"),
              "io_MS_Global");

    check_str(12, "resolve_prior_from_model(\"model_MS_Global_a1nl_etaa3_HarveyLike\") == \"io_MS_Global\"",
              resolve_prior_from_model("model_MS_Global_a1nl_etaa3_HarveyLike"),
              "io_MS_Global");

    check_str(13, "resolve_prior_from_model(\"model_MS_Global_a1etaa3_HarveyLike_Classic_v2\") == \"io_MS_Global\"",
              resolve_prior_from_model("model_MS_Global_a1etaa3_HarveyLike_Classic_v2"),
              "io_MS_Global");

    check_str(14, "resolve_prior_from_model(\"model_MS_Global_a1etaa3_HarveyLike_Classic_v3\") == \"io_MS_Global\"",
              resolve_prior_from_model("model_MS_Global_a1etaa3_HarveyLike_Classic_v3"),
              "io_MS_Global");

    check_str(15, "resolve_prior_from_model(\"model_MS_Global_a1n_a2a3_HarveyLike\") == \"io_MS_Global\"",
              resolve_prior_from_model("model_MS_Global_a1n_a2a3_HarveyLike"),
              "io_MS_Global");

    check_str(16, "resolve_prior_from_model(\"model_MS_Global_a1nl_a2a3_HarveyLike\") == \"io_MS_Global\"",
              resolve_prior_from_model("model_MS_Global_a1nl_a2a3_HarveyLike"),
              "io_MS_Global");

    check_str(17, "resolve_prior_from_model(\"model_MS_Global_ajAlm_HarveyLike\") == \"io_MS_Global\"",
              resolve_prior_from_model("model_MS_Global_ajAlm_HarveyLike"),
              "io_MS_Global");

    check_str(18, "resolve_prior_from_model(\"model_MS_Global_aj_HarveyLike\") == \"io_MS_Global\"",
              resolve_prior_from_model("model_MS_Global_aj_HarveyLike"),
              "io_MS_Global");

    check_str(19, "resolve_prior_from_model(\"model_RGB_asympt_aj_AppWidth_HarveyLike_v4\") == \"io_asymptotic\"",
              resolve_prior_from_model("model_RGB_asympt_aj_AppWidth_HarveyLike_v4"),
              "io_asymptotic");

    check_str(20, "resolve_prior_from_model(\"model_RGB_asympt_aj_CteWidth_HarveyLike_v4\") == \"io_asymptotic\"",
              resolve_prior_from_model("model_RGB_asympt_aj_CteWidth_HarveyLike_v4"),
              "io_asymptotic");

    check_str(21, "resolve_prior_from_model(\"model_MS_local_basic\") == \"io_local\"",
              resolve_prior_from_model("model_MS_local_basic"),
              "io_local");

    check_str(22, "resolve_prior_from_model(\"model_MS_local_Hnlm\") == \"io_local\"",
              resolve_prior_from_model("model_MS_local_Hnlm"),
              "io_local");

    check_str(23, "resolve_prior_from_model(\"model_ajfit\") == \"io_ajfit\"",
              resolve_prior_from_model("model_ajfit"),
              "io_ajfit");

    check_str(24, "resolve_prior_from_model(\"model_UNKNOWN\") == \"\"",
              resolve_prior_from_model("model_UNKNOWN"),
              "");

    check_str(25, "resolve_prior_from_model(\"\") == \"\"",
              resolve_prior_from_model(""),
              "");

    // ---- supported_model_fullnames ----
    check_size(26, "supported_model_fullnames().size() == 16",
               supported_model_fullnames().size(),
               16u);

    // ---- Summary ----
    if (failed_tests == 0) {
        std::cout << "All " << total_tests << " tests passed\n";
        return 0;
    } else {
        std::cout << failed_tests << " tests FAILED\n";
        return 1;
    }
}
