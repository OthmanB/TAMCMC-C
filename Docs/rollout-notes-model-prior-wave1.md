# Rollout Notes: Model-Prior Integration (Wave 1)

## 1. Summary
This rollout (Wave 1) delivers the core infrastructure for model-to-prior integration, focusing on resolving conflicts between configuration files and `.model` files.

### Wave 1 Scope
- **Conflict Resolution Logic**: Implements `.model`-wins logic when both `config_default.cfg` and the `.model` file specify a model name.
- **Strict Plotting Preflight**: Enforces gnuplot availability when diagnostic flags are enabled.

### Explicitly Out of Scope
- **Heuristic Mapping**: The automatic mapping of `model_fullname` to `prior_fct_name` is deferred to a follow-up rollout.

## 2. Legacy Field: `model_fct_name`
The `model_fct_name` field in `config_default.cfg` is now considered a **legacy field**. 

### Guidance
Users should **comment out** `model_fct_name` in their configuration files. The `model_fullname` field within the `.model` file now takes precedence. 

### Conflict Behavior
When both are present, `cpptamcmc` emits the following warning and uses the value from the `.model` file:

```
WARNING: model_fct_name in config file is legacy and should be commented out. Using model_fullname from .model file instead.
         Config model_fct_name: <config-value>
         .model model_fullname: <model-file-value>
```

The config value is discarded. If the `.model` file lacks a `model_fullname` (or it is blank), the config `model_fct_name` serves as a fallback, and no warning is issued.

## 3. Migration Guide
To migrate your configuration:
1. Open your `config_default.cfg` or project-specific config file.
2. Locate `model_fct_name` and **comment it out** using `#`.
3. Verify that your `.model` file contains the correct model name in the `model_fullname` field:
   ```
   model_fullname = <correct-model-name>
   ```

## 4. Deferred: `model_fullname` → `prior_fct_name` Heuristic
The automatic heuristic for mapping models to priors is **deferred** and is not part of Wave 1.
- See the follow-up specification (T14) for details on this planned feature.
- The behavior of `prior_fct_name` remains **unchanged** in this wave. Users must still specify `prior_fct_name` manually.

## 5. Strict Plotting Preflight
Wave 1 introduces a strict preflight check for plotting. If any of the following 6 diagnostic flags are set to `1`, `gnuplot` must be available in the system path:
- `chains_diags`
- `evidence_diags`
- `pdfs_diags`
- `model_initial_diags`
- `model_buffer_diags`
- `model_final_diags`

On failure, `cpptamcmc` will print an error listing all 6 flags and instructions to either install `gnuplot` or disable the diagnostics.
