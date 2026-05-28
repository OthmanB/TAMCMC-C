# Follow-up Spec: `model_fullname` to `prior_fct_name` Heuristic Mapping

## 1. Executive Summary
This specification outlines the proposed auto-mapping heuristic for deriving `prior_fct_name` from `model_fullname`. 
This functionality is DEFERRED and is NOT part of the Wave 1 rollout. 
In Wave 1, `prior_fct_name` behavior remains unchanged and must be explicitly set or it will follow existing defaults.

## 2. Out-of-Scope Statement (Wave 1)
- **NOT implemented in this rollout.**
- Wave 1 focuses on `.model`-wins conflict detection and strict plotting preflight verification.
- The semantics of `prior_fct_name` are UNCHANGED in Wave 1.

## 3. Proposed Heuristic
When `.model` `model_fullname` is set and `prior_fct_name` is absent or blank in the `.priors` file, the system should auto-derive the appropriate `prior_fct_name`.

### Mapping Proposals
- `model_Harvey_Gaussian` → `priors_Harvey_Gaussian`
- `model_Kallinger2014_Gaussian` → `priors_Kallinger2014_Gaussian`
- `model_MS_Global_*` (any variant) → `io_MS_Global`
- `model_RGB_asympt_aj_AppWidth_HarveyLike_v4` → `io_asymptotic`
- `model_RGB_asympt_aj_CteWidth_HarveyLike_v4` → `io_asymptotic`

#### `io_local` Policy Proposal
`io_local` is operationally tied to local/slice-by-slice workflows (high-HNR and/or crowded-mode use cases) rather than a broad family wildcard. To avoid accidental misrouting:
- `model_MS_local_basic` → `io_local`
- `model_MS_local_Hnlm` → `io_local`
- Any other model should **NOT** auto-map to `io_local`.
- If users intend local analysis with a non-local model family, they must set `prior_fct_name` explicitly.

- **Unknown models**: The system must FAIL with a clear error message. Silent fallbacks are prohibited to avoid incorrect prior associations.

## 4. Risks
- **Inconsistent Naming**: Model and prior names don't always follow a 1:1 derivable pattern.
- **Silent Errors**: Auto-mapping might incorrectly pair a model with a prior function if rules are too loose.
- **Ambiguity**: Must be strictly opt-in or fail-fast.
- **ajfit Exemption**: The `ajfit` path uses distinct naming conventions and must remain exempt from this heuristic to prevent regression.

## 5. Precedence Proposal (Future Wave)
The following priority order is proposed for determining the prior function:
1. **Explicit `prior_fct_name`**: User setting in `.priors` always wins.
2. **Heuristic Mapping**: Derived from `model_fullname` if (1) is missing.
3. **Error**: If neither is available and no mapping exists.

When a heuristic is applied, the system should emit an INFO message to the log. If the heuristic fails to map, it must produce an actionable error listing the required `prior_fct_name`.

## 6. Implementation Notes
- **Insertion Point**: Logic should be placed after the conflict check in `read_inputs_priors_MS_Global()` (approx. `config.cpp:700`).
- **Mechanism**: Requires a lookup table or pattern-matching engine for `model_*` → `prior_*` mappings.
- **ajfit Path**: Ensure `prior_fct_name` logic remains untouched for `ajfit` workflows.
- **User Intent**: The heuristic must NOT overwrite an explicitly set `prior_fct_name`.
