# IMPLEMENTED (Wave 2): `model_fullname` to `prior_fct_name` Heuristic Mapping

## Implementation Notes (Wave 2)
Implemented in Wave 2. Key files:
- `TAMCMC-C/tamcmc/sources/prior_auto_map.h` / `.cpp` — mapping table and helper functions
- `TAMCMC-C/tamcmc/sources/config.cpp` — dispatch logic at Config::read_inputs_files()
- Trigger: `prior_fct_name=auto;` in config file

## 1. Executive Summary
This specification outlines the auto-mapping heuristic for deriving `prior_fct_name` from `model_fullname`. 
This functionality is implemented in Wave 2.
In Wave 1, `prior_fct_name` behavior remained unchanged and required explicit settings.

## 2. Status
- **Implemented in Wave 2.**
- Wave 1 focused on `.model`-wins conflict detection and strict plotting preflight verification.
- The semantics of `prior_fct_name` are now extended with the `auto` keyword.

## 3. Heuristic Mapping
When `prior_fct_name=auto;` is set in the config file, the system auto-derives the appropriate `prior_fct_name` from the `.model` file's `model_fullname`.

### Mapping Table
| model_fullname | prior_fct_name |
|---|---|
| model_Harvey_Gaussian | priors_Harvey_Gaussian |
| model_Kallinger2014_Gaussian | priors_Kallinger2014_Gaussian |
| model_MS_Global_a1l_etaa3_HarveyLike | io_MS_Global |
| model_MS_Global_a1n_etaa3_HarveyLike | io_MS_Global |
| model_MS_Global_a1nl_etaa3_HarveyLike | io_MS_Global |
| model_MS_Global_a1etaa3_HarveyLike_Classic_v2 | io_MS_Global |
| model_MS_Global_a1etaa3_HarveyLike_Classic_v3 | io_MS_Global |
| model_MS_Global_a1n_a2a3_HarveyLike | io_MS_Global |
| model_MS_Global_a1nl_a2a3_HarveyLike | io_MS_Global |
| model_MS_Global_ajAlm_HarveyLike | io_MS_Global |
| model_MS_Global_aj_HarveyLike | io_MS_Global |
| model_RGB_asympt_aj_AppWidth_HarveyLike_v4 | io_asymptotic |
| model_RGB_asympt_aj_CteWidth_HarveyLike_v4 | io_asymptotic |
| model_MS_local_basic | io_local |
| model_MS_local_Hnlm | io_local |

#### `io_local` Policy
`io_local` is operationally tied to local/slice-by-slice workflows (high-HNR and/or crowded-mode use cases) rather than a broad family wildcard. To avoid accidental misrouting:
- `model_MS_local_basic` → `io_local`
- `model_MS_local_Hnlm` → `io_local`
- Any other model should **NOT** auto-map to `io_local`.
- If users intend local analysis with a non-local model family, they must set `prior_fct_name` explicitly.

- **Unknown models**: The system FAILS with a clear error message. Silent fallbacks are prohibited to avoid incorrect prior associations.

## 4. Risks
- **Inconsistent Naming**: Model and prior names don't always follow a 1:1 derivable pattern.
- **Silent Errors**: Auto-mapping might incorrectly pair a model with a prior function if rules are too loose.
- **Ambiguity**: Must be strictly opt-in via `prior_fct_name=auto;`.
- **ajfit Exemption**: The `ajfit` path uses distinct naming conventions and is exempt from this heuristic to prevent regression.

## 5. Precedence
The following priority order is used for determining the prior function:
1. **Explicit `prior_fct_name` (not "auto")**: User setting in `.priors` always wins.
2. **Heuristic Mapping**: Derived from `model_fullname` if `prior_fct_name=auto;`.
3. **Error**: If neither is available and no mapping exists.

When a heuristic is applied, the system emits an INFO message to the log. If the heuristic fails to map, it produces an actionable error listing the required `prior_fct_name`.

## 6. Implementation Notes
- **Insertion Point**: Dispatch logic is in `Config::read_inputs_files()` in `config.cpp`.
- **Mechanism**: Uses a lookup table in `prior_auto_map.cpp` for `model_*` → `prior_*` mappings.
- **ajfit Path**: `prior_fct_name=auto;` is rejected for `ajfit` workflows.
- **User Intent**: The heuristic does NOT overwrite an explicitly set `prior_fct_name` (unless it is "auto").
