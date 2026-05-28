# Auto-Prior Detection Guide

This document describes the auto-detection feature for `prior_fct_name`, implemented in Wave 2. This feature allows the system to automatically resolve the correct prior function based on the model being used.

## Section 1: How to opt in

To enable auto-prior detection, set the `prior_fct_name` to `auto` in your configuration file:

```text
# In config file:
prior_fct_name=auto;
```

**How it works:**
1. The system reads the `.model` file associated with the run.
2. It extracts the `model_fullname`.
3. It looks up this name in a curated mapping table.
4. If a match is found, it uses the corresponding `prior_fct_name`.

**Important Notes:**
- The keyword `auto` is **CASE-SENSITIVE**. Using `AUTO` or `Auto` will be treated as an unknown prior name and will trigger an error.
- This feature is strictly opt-in. If `prior_fct_name` is set to anything else (except empty), the explicit value is used.

## Section 2: Complete mapping table

The following table lists all 15 supported models and their corresponding prior functions:

| model_fullname | prior_fct_name | Notes |
|---|---|---|
| model_Harvey_Gaussian | priors_Harvey_Gaussian | Harvey Gaussian fit |
| model_Kallinger2014_Gaussian | priors_Kallinger2014_Gaussian | Kallinger 2014 Gaussian fit |
| model_MS_Global_a1l_etaa3_HarveyLike | io_MS_Global | MS Global (case N from models_ctrl.list) |
| model_MS_Global_a1n_etaa3_HarveyLike | io_MS_Global | MS Global (case N from models_ctrl.list) |
| model_MS_Global_a1nl_etaa3_HarveyLike | io_MS_Global | MS Global (case N from models_ctrl.list) |
| model_MS_Global_a1etaa3_HarveyLike_Classic_v2 | io_MS_Global | MS Global (case N from models_ctrl.list) |
| model_MS_Global_a1etaa3_HarveyLike_Classic_v3 | io_MS_Global | MS Global (case N from models_ctrl.list) |
| model_MS_Global_a1n_a2a3_HarveyLike | io_MS_Global | MS Global (case N from models_ctrl.list) |
| model_MS_Global_a1nl_a2a3_HarveyLike | io_MS_Global | MS Global (case N from models_ctrl.list) |
| model_MS_Global_ajAlm_HarveyLike | io_MS_Global | MS Global (case N from models_ctrl.list) |
| model_MS_Global_aj_HarveyLike | io_MS_Global | MS Global (case N from models_ctrl.list) |
| model_RGB_asympt_aj_AppWidth_HarveyLike_v4 | io_asymptotic | RGB asymptotic (AppWidth) |
| model_RGB_asympt_aj_CteWidth_HarveyLike_v4 | io_asymptotic | RGB asymptotic (CteWidth) |
| model_MS_local_basic | io_local | MS Local basic |
| model_MS_local_Hnlm | io_local | MS Local Hnlm |

## Section 3: Error scenarios

The system is designed to fail fast with actionable messages in the following cases:

- **`.model` file missing or `model_fullname` undefined**: The system will emit a FATAL error indicating that auto-resolution is impossible without a model definition.
- **`model_fullname` not in the table**: If the model is not among the 15 listed above, the system will emit a FATAL error listing all currently supported models.
- **`prior_fct_name=auto;` in ajfit workflow**: This is not supported. The system will fail with: `not supported for ajfit; use io_ajfit explicitly`.
- **Wrong case for `auto`**: Using `AUTO` or `Auto` will trigger the standard unknown prior name error (typically at `config.cpp:2106`).

## Section 4: Mismatch behavior

When a user provides an explicit `prior_fct_name` (other than `auto`), the system still checks it against the `model_fullname` for consistency:

- **Conflict Detected**: If `prior_fct_name=X` but `model_fullname=M` typically maps to `Y` (where `X != Y`), the system emits a WARNING on stderr:
  `WARNING: prior_fct_name=X is explicitly set, but model_fullname=M typically uses Y`
- **User Preference**: Execution proceeds using the user-specified `X`. The user's explicit setting always wins. No data is modified, and the run does not fail.

## Section 5: Worked examples

### Example 1 — Harvey Gaussian
```text
# In config file:
prior_fct_name=auto;

# In .model file:
model_fullname              model_Harvey_Gaussian

# Result in log:
INFO: prior_fct_name auto-resolved from model_fullname=model_Harvey_Gaussian to priors_Harvey_Gaussian
```

### Example 2 — Kallinger
```text
# In config file:
prior_fct_name=auto;

# In .model file:
model_fullname              model_Kallinger2014_Gaussian

# Result in log:
INFO: prior_fct_name auto-resolved from model_fullname=model_Kallinger2014_Gaussian to priors_Kallinger2014_Gaussian
```

### Example 3 — MS Global
```text
# In config file:
prior_fct_name=auto;

# In .model file:
model_fullname              model_MS_Global_ajAlm_HarveyLike

# Result in log:
INFO: prior_fct_name auto-resolved from model_fullname=model_MS_Global_ajAlm_HarveyLike to io_MS_Global
```

### Example 4 — RGB Asymptotic
```text
# In config file:
prior_fct_name=auto;

# In .model file:
model_fullname              model_RGB_asympt_aj_AppWidth_HarveyLike_v4

# Result in log:
INFO: prior_fct_name auto-resolved from model_fullname=model_RGB_asympt_aj_AppWidth_HarveyLike_v4 to io_asymptotic
```

## Section 6: Template files

Three template files demonstrate the use of the `auto` mode. They are located in the templates directory:

- `.cfg.MS_Global_fit_auto`
- `.cfg.MS_local_fit_auto`
- `.cfg.RGB_asympt_fit_auto`
