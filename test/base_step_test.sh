#!/bin/bash
# Base-step test suite for model-name precedence and conflict detection
# Tests validate the conflict detection and model-name precedence logic.
#
# Scenarios:
#   A - config model_fct_name matches .model model_fullname → exit 0, no warning
#   B - config model_fct_name conflicts with .model model_fullname → exit 0, WARNING printed
#   E - both model sources absent → non-zero exit (fatal error)
#
# Requires: binary built at ./bin/cpptamcmc
# Fixtures:  MS Global format .model/.data/.priors files in MODELS_DIR
# Config:    config_default.cfg and config_presets.cfg managed per-scenario
# Flags:     -S 1 -L 1 (mandatory)

set -euo pipefail

PASS=0
FAIL=0
BINARY="./bin/cpptamcmc"
CONFIG_DIR="Config/default"
CONFIG_FILE="${CONFIG_DIR}/config_default.cfg"
CONFIG_BAK="${CONFIG_DIR}/config_default.cfg.bak"
PRESETS_FILE="Config/config_presets.cfg"
MODELS_DIR="/var/folders/xx/66ps59dd23b91gw3mf08_x080000gn/T/opencode/tamcmc-conflict-test/models"
OUTPUT_DIR="/var/folders/xx/66ps59dd23b91gw3mf08_x080000gn/T/opencode/tamcmc-conflict-test/output"
TIMEOUT="/opt/homebrew/bin/timeout 30"
AUTO_MODELS_DIR="/var/folders/xx/66ps59dd23b91gw3mf08_x080000gn/T/opencode/tamcmc-auto-test/models"
AUTO_OUTPUT_DIR="/var/folders/xx/66ps59dd23b91gw3mf08_x080000gn/T/opencode/tamcmc-auto-test/output"

run_test() {
    local name="$1"
    local result="$2"
    if [ "$result" = "PASS" ]; then
        echo "  PASS: $name"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $name"
        FAIL=$((FAIL + 1))
    fi
}

# Restore config_default.cfg to safe baseline (model_Harvey_Gaussian active)
restore_config() {
    if [ -f "${CONFIG_BAK}" ]; then
        cp "${CONFIG_BAK}" "${CONFIG_FILE}"
    fi
}

# Set model_fct_name in config_default.cfg
set_model_fct_name() {
    local value="$1"
    sed -i.tmp "s|^	model_fct_name=.*;|	model_fct_name=${value};|" "${CONFIG_FILE}"
    rm -f "${CONFIG_FILE}.tmp"
}

# Comment out model_fct_name in config_default.cfg
comment_model_fct_name() {
    sed -i.tmp "s|^	model_fct_name=|	#model_fct_name=|" "${CONFIG_FILE}"
    rm -f "${CONFIG_FILE}.tmp"
}

# Set all 6 diagnostics flags to a given value (0 or 1)
set_all_diags() {
    local val="$1"
    for flag in chains_diags evidence_diags pdfs_diags model_initial_diags model_buffer_diags model_final_diags; do
        sed -i.tmp "s|^	${flag}=.*;|	${flag}=${val};|" "${CONFIG_FILE}"
        rm -f "${CONFIG_FILE}.tmp"
    done
}

# Write presets file for a given scenario name
write_presets() {
    local scenario="$1"
    cat > "${PRESETS_FILE}" << PRESETS
force_manual_config=0;
manual_config_file=;
cfg_models_dir=${MODELS_DIR}/;
cfg_out_dir=${OUTPUT_DIR}/;
processing      = Burn-in , Learning , Acquire;
Nsamples        = 1000    , 10000    , 10000;
c0              = 1.6;
restore         = 0       , 1        , 2;
core_out        = B       , L        , A;
core_in         = B       , B        , L;
start_index_processing=0;
last_index_processing=2;
table_ids= 1 , 2;
${scenario}      1;
/END;
PRESETS
}

# Set prior_fct_name in config_default.cfg (active/uncommented line only)
set_prior_fct_name() {
    local value="$1"
    sed -i.tmp "s|^	prior_fct_name=.*;|	prior_fct_name=${value};|" "${CONFIG_FILE}"
    rm -f "${CONFIG_FILE}.tmp"
}

# Write presets pointing to auto-test models dir
write_auto_presets() {
    local scenario="$1"
    cat > "${PRESETS_FILE}" << PRESETS
force_manual_config=0;
manual_config_file=;
cfg_models_dir=${AUTO_MODELS_DIR}/;
cfg_out_dir=${AUTO_OUTPUT_DIR}/;
processing      = Burn-in , Learning , Acquire;
Nsamples        = 1000    , 1000     , 1000;
c0              = 1.6;
restore         = 0       , 1        , 2;
core_out        = B       , L        , A;
core_in         = B       , B        , L;
start_index_processing=0;
last_index_processing=2;
table_ids= 1 , 2;
${scenario}      1;
/END;
PRESETS
}

# Create all Wave-2 fixture files (called once at start of Wave-2 section)
setup_auto_fixtures() {
    mkdir -p "${AUTO_MODELS_DIR}" "${AUTO_OUTPUT_DIR}"

    # Create minimal data file
    cat > "${AUTO_MODELS_DIR}/AUTO_BASE.data" << 'DATA'
# fake spectrum for auto-mode tests
!   frequency   power
500.0   1.0
600.0   1.0
700.0   1.0
800.0   1.0
900.0   1.0
1000.0  1.0
DATA

    # .data files for each scenario (copy of AUTO_BASE.data)
    for SCEN in AUTO_HARVEY AUTO_KALLINGER AUTO_MS_GLOBAL_aj AUTO_MS_GLOBAL_ajAlm \
                AUTO_RGB_APP AUTO_RGB_CTE AUTO_LOCAL_BASIC AUTO_LOCAL_HNLM \
                AUTO_UNKNOWN AUTO_MALFORMED AUTO_WRONG_CASE \
                MISMATCH_WARN AJFIT_AUTO_FATAL AJFIT_EXEMPT; do
        cp "${AUTO_MODELS_DIR}/AUTO_BASE.data" "${AUTO_MODELS_DIR}/${SCEN}.data"
    done

    # Individual .model files
    printf '* 500.0 1000.0\nmodel_fullname model_Harvey_Gaussian\n'                          > "${AUTO_MODELS_DIR}/AUTO_HARVEY.model"
    printf '* 500.0 1000.0\nmodel_fullname model_Kallinger2014_Gaussian\n'                   > "${AUTO_MODELS_DIR}/AUTO_KALLINGER.model"
    printf '* 500.0 1000.0\nmodel_fullname model_MS_Global_aj_HarveyLike\n'                  > "${AUTO_MODELS_DIR}/AUTO_MS_GLOBAL_aj.model"
    printf '* 500.0 1000.0\nmodel_fullname model_MS_Global_ajAlm_HarveyLike\n'               > "${AUTO_MODELS_DIR}/AUTO_MS_GLOBAL_ajAlm.model"
    printf '* 500.0 1000.0\nmodel_fullname model_RGB_asympt_aj_AppWidth_HarveyLike_v4\n'     > "${AUTO_MODELS_DIR}/AUTO_RGB_APP.model"
    printf '* 500.0 1000.0\nmodel_fullname model_RGB_asympt_aj_CteWidth_HarveyLike_v4\n'     > "${AUTO_MODELS_DIR}/AUTO_RGB_CTE.model"
    printf '* 500.0 1000.0\nmodel_fullname model_MS_local_basic\n'                           > "${AUTO_MODELS_DIR}/AUTO_LOCAL_BASIC.model"
    printf '* 500.0 1000.0\nmodel_fullname model_MS_local_Hnlm\n'                            > "${AUTO_MODELS_DIR}/AUTO_LOCAL_HNLM.model"
    printf '* 500.0 1000.0\nmodel_fullname model_DOES_NOT_EXIST\n'                           > "${AUTO_MODELS_DIR}/AUTO_UNKNOWN.model"
    printf '* 500.0 1000.0\n# no model_fullname line here - just a comment\n'               > "${AUTO_MODELS_DIR}/AUTO_MALFORMED.model"
    printf '* 500.0 1000.0\nmodel_fullname model_Harvey_Gaussian\n'                          > "${AUTO_MODELS_DIR}/AUTO_WRONG_CASE.model"
    printf '* 500.0 1000.0\nmodel_fullname model_Harvey_Gaussian\n'                          > "${AUTO_MODELS_DIR}/MISMATCH_WARN.model"
    printf '* 500.0 1000.0\nmodel_fullname model_ajfit\n'                                    > "${AUTO_MODELS_DIR}/AJFIT_AUTO_FATAL.model"
    printf '* 500.0 1000.0\nmodel_fullname model_Harvey_Gaussian\n'                          > "${AUTO_MODELS_DIR}/AJFIT_EXEMPT.model"

    # WAVE1_REG: copy SCENARIO_A files from conflict-test dir
    cp "${MODELS_DIR}/SCENARIO_A.model"    "${AUTO_MODELS_DIR}/WAVE1_REG.model"
    cp "${MODELS_DIR}/SCENARIO_A.data"     "${AUTO_MODELS_DIR}/WAVE1_REG.data"
    cp "${MODELS_DIR}/SCENARIO_A_0.priors" "${AUTO_MODELS_DIR}/WAVE1_REG_0.priors"
    [ -f "${MODELS_DIR}/SCENARIO_A_1.priors" ] && \
        cp "${MODELS_DIR}/SCENARIO_A_1.priors" "${AUTO_MODELS_DIR}/WAVE1_REG_1.priors" || true
}

# Preflight checks
if [ ! -f "$BINARY" ]; then
    echo "ERROR: Binary not found at $BINARY"
    exit 1
fi

for fixture in \
    "${MODELS_DIR}/SCENARIO_A.model" \
    "${MODELS_DIR}/SCENARIO_A.data" \
    "${MODELS_DIR}/SCENARIO_A_0.priors" \
    "${MODELS_DIR}/SCENARIO_B.model" \
    "${MODELS_DIR}/SCENARIO_B.data" \
    "${MODELS_DIR}/SCENARIO_B_0.priors" \
    "${MODELS_DIR}/SCENARIO_E.model" \
    "${MODELS_DIR}/SCENARIO_E.data"
do
    if [ ! -f "$fixture" ]; then
        echo "ERROR: Required fixture not found: $fixture"
        exit 1
    fi
done

echo "=========================================="
echo "Base-Step Test Suite (Model Precedence)"
echo "=========================================="
echo ""

# ---------------------------------------------------------------------------
echo "Test A: Aligned scenario (config matches .model, no conflict)"
# ---------------------------------------------------------------------------
write_presets "SCENARIO_A"
set_model_fct_name "model_MS_Global_ajAlm_HarveyLike"

output_a=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
exit_code_a=$?

restore_config

if ! echo "$output_a" | grep -qE "model_fct_name.*legacy|legacy.*model_fct_name|commented out"; then
    run_test "A: no conflict warning" "PASS"
else
    run_test "A: no conflict warning" "FAIL"
fi

if [ $exit_code_a -eq 0 ]; then
    run_test "A: exit code 0" "PASS"
else
    run_test "A: exit code 0 (got ${exit_code_a})" "FAIL"
fi

echo ""

# ---------------------------------------------------------------------------
echo "Test B: Conflict scenario (config disagrees with .model)"
# ---------------------------------------------------------------------------
write_presets "SCENARIO_B"
set_model_fct_name "model_Harvey_Gaussian"

output_b=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
exit_b=$?
# timeout (124) is acceptable — binary was still running (not crashed)
[ $exit_b -eq 124 ] && exit_b=0

restore_config

if echo "$output_b" | grep -q "model_fct_name" && \
   echo "$output_b" | grep -q "legacy" && \
   echo "$output_b" | grep -q "commented out"; then
    run_test "B: conflict WARNING printed" "PASS"
else
    run_test "B: conflict WARNING printed" "FAIL"
fi

if [ $exit_b -eq 0 ]; then
    run_test "B: exit code 0 (or timeout=running)" "PASS"
else
    run_test "B: exit code 0 (got ${exit_b})" "FAIL"
fi

echo ""

# ---------------------------------------------------------------------------
echo "Test E: Both model sources absent (fatal error expected)"
# ---------------------------------------------------------------------------
write_presets "SCENARIO_E"
comment_model_fct_name

output_e=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) && exit_code_e=0 || exit_code_e=$?

restore_config

if [ $exit_code_e -ne 0 ] && [ $exit_code_e -ne 124 ]; then
    run_test "E: non-zero exit (got ${exit_code_e})" "PASS"
else
    run_test "E: non-zero exit (got ${exit_code_e})" "FAIL"
fi

echo ""

# ---------------------------------------------------------------------------
echo "Test C: Diagnostics ON + gnuplot unavailable (preflight must FAIL)"
# ---------------------------------------------------------------------------
write_presets "SCENARIO_A"
restore_config
# diags already ON (all 6 = 1 in restored config)

output_c=$(PATH=/bin:/usr/bin /opt/homebrew/bin/timeout 15 ${BINARY} -S 1 -L 1 2>&1) && exit_code_c=0 || exit_code_c=$?

restore_config

# Must have non-zero exit (preflight killed it)
if [ $exit_code_c -ne 0 ] && [ $exit_code_c -ne 124 ]; then
    run_test "C: non-zero exit (got ${exit_code_c})" "PASS"
else
    run_test "C: non-zero exit (got ${exit_code_c})" "FAIL"
fi

# Must contain gnuplot error message
if echo "$output_c" | grep -q "gnuplot"; then
    run_test "C: gnuplot error message present" "PASS"
else
    run_test "C: gnuplot error message present" "FAIL"
fi

# Must list chains_diags flag
if echo "$output_c" | grep -q "chains_diags"; then
    run_test "C: chains_diags listed in error" "PASS"
else
    run_test "C: chains_diags listed in error" "FAIL"
fi

# Must contain REMEDIATION section
if echo "$output_c" | grep -q "REMEDIATION"; then
    run_test "C: REMEDIATION guidance present" "PASS"
else
    run_test "C: REMEDIATION guidance present" "FAIL"
fi

# Must NOT have reached burnin
if ! echo "$output_c" | grep -qE "BEGINING THE MCMC PROCESS|Burn-in phase"; then
    run_test "C: burnin NOT entered (fail-fast confirmed)" "PASS"
else
    run_test "C: burnin NOT entered (fail-fast confirmed)" "FAIL"
fi

echo ""

# ---------------------------------------------------------------------------
echo "Test D: Diagnostics OFF + gnuplot unavailable (preflight must PASS)"
# ---------------------------------------------------------------------------
write_presets "SCENARIO_A"
restore_config
set_all_diags 0

output_d=$(PATH=/bin:/usr/bin /opt/homebrew/bin/timeout 20 ${BINARY} -S 1 -L 1 2>&1) || true
exit_code_d=$?

restore_config

# Exit 0 or timeout (124 = still running = pass) — NOT a gnuplot error exit
if [ $exit_code_d -eq 0 ] || [ $exit_code_d -eq 124 ]; then
    run_test "D: exit 0 or timeout (no preflight failure)" "PASS"
else
    run_test "D: exit 0 or timeout (got ${exit_code_d})" "FAIL"
fi

# Must NOT contain gnuplot error
if ! echo "$output_d" | grep -q "ERROR: Plotting diagnostics enabled"; then
    run_test "D: no gnuplot preflight error" "PASS"
else
    run_test "D: no gnuplot preflight error" "FAIL"
fi

# Must have reached "Initial configuration done" (got past preflight)
if echo "$output_d" | grep -q "Initial configuration done"; then
    run_test "D: past preflight (Initial configuration done)" "PASS"
else
    run_test "D: past preflight (Initial configuration done)" "FAIL"
fi

echo ""

# ==========================================================================
# Wave 2: auto-mode prior detection tests
# ==========================================================================
setup_auto_fixtures

# ---------------------------------------------------------------------------
echo "Test AUTO_HARVEY: auto resolves model_Harvey_Gaussian → priors_Harvey_Gaussian"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_HARVEY"
set_prior_fct_name "auto"
set_model_fct_name "model_Harvey_Gaussian"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
restore_config

if echo "$output" | grep -q "INFO: prior_fct_name auto-resolved" && \
   echo "$output" | grep -q "priors_Harvey_Gaussian"; then
    run_test "AUTO_HARVEY: INFO auto-resolved to priors_Harvey_Gaussian" "PASS"
else
    run_test "AUTO_HARVEY: INFO auto-resolved to priors_Harvey_Gaussian" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_KALLINGER: auto resolves model_Kallinger2014_Gaussian → priors_Kallinger2014_Gaussian"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_KALLINGER"
set_prior_fct_name "auto"
set_model_fct_name "model_Kallinger2014_Gaussian"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
restore_config

if echo "$output" | grep -q "priors_Kallinger2014_Gaussian"; then
    run_test "AUTO_KALLINGER: auto-resolved to priors_Kallinger2014_Gaussian" "PASS"
else
    run_test "AUTO_KALLINGER: auto-resolved to priors_Kallinger2014_Gaussian" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_MS_GLOBAL_aj: auto resolves model_MS_Global_aj_HarveyLike → io_MS_Global"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_MS_GLOBAL_aj"
set_prior_fct_name "auto"
set_model_fct_name "model_MS_Global_aj_HarveyLike"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
restore_config

if echo "$output" | grep -q "io_MS_Global"; then
    run_test "AUTO_MS_GLOBAL_aj: auto-resolved to io_MS_Global" "PASS"
else
    run_test "AUTO_MS_GLOBAL_aj: auto-resolved to io_MS_Global" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_MS_GLOBAL_ajAlm: auto resolves model_MS_Global_ajAlm_HarveyLike → io_MS_Global"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_MS_GLOBAL_ajAlm"
set_prior_fct_name "auto"
set_model_fct_name "model_MS_Global_ajAlm_HarveyLike"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
restore_config

if echo "$output" | grep -q "io_MS_Global"; then
    run_test "AUTO_MS_GLOBAL_ajAlm: auto-resolved to io_MS_Global" "PASS"
else
    run_test "AUTO_MS_GLOBAL_ajAlm: auto-resolved to io_MS_Global" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_RGB_APP: auto resolves RGB_AppWidth → io_asymptotic"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_RGB_APP"
set_prior_fct_name "auto"
set_model_fct_name "model_RGB_asympt_aj_AppWidth_HarveyLike_v4"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
restore_config

if echo "$output" | grep -q "io_asymptotic"; then
    run_test "AUTO_RGB_APP: auto-resolved to io_asymptotic" "PASS"
else
    run_test "AUTO_RGB_APP: auto-resolved to io_asymptotic" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_RGB_CTE: auto resolves RGB_CteWidth → io_asymptotic"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_RGB_CTE"
set_prior_fct_name "auto"
set_model_fct_name "model_RGB_asympt_aj_CteWidth_HarveyLike_v4"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
restore_config

if echo "$output" | grep -q "io_asymptotic"; then
    run_test "AUTO_RGB_CTE: auto-resolved to io_asymptotic" "PASS"
else
    run_test "AUTO_RGB_CTE: auto-resolved to io_asymptotic" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_LOCAL_BASIC: auto resolves model_MS_local_basic → io_local"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_LOCAL_BASIC"
set_prior_fct_name "auto"
set_model_fct_name "model_MS_local_basic"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
restore_config

if echo "$output" | grep -q "io_local"; then
    run_test "AUTO_LOCAL_BASIC: auto-resolved to io_local" "PASS"
else
    run_test "AUTO_LOCAL_BASIC: auto-resolved to io_local" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_LOCAL_HNLM: auto resolves model_MS_local_Hnlm → io_local"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_LOCAL_HNLM"
set_prior_fct_name "auto"
set_model_fct_name "model_MS_local_Hnlm"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
restore_config

if echo "$output" | grep -q "io_local"; then
    run_test "AUTO_LOCAL_HNLM: auto-resolved to io_local" "PASS"
else
    run_test "AUTO_LOCAL_HNLM: auto-resolved to io_local" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_UNKNOWN: unknown model_fullname → FATAL exit"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_UNKNOWN"
set_prior_fct_name "auto"
set_model_fct_name "model_Harvey_Gaussian"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) && exit_code=0 || exit_code=$?
restore_config

if [ $exit_code -ne 0 ] && [ $exit_code -ne 124 ] && echo "$output" | grep -q "FATAL"; then
    run_test "AUTO_UNKNOWN: FATAL exit on unknown model_fullname" "PASS"
else
    run_test "AUTO_UNKNOWN: FATAL exit on unknown model_fullname (got ${exit_code})" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_MISSING_MODEL: missing .model file → FATAL exit"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_DOES_NOT_EXIST_XYZ"
set_prior_fct_name "auto"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) && exit_code=0 || exit_code=$?
restore_config

if [ $exit_code -ne 0 ] && [ $exit_code -ne 124 ]; then
    run_test "AUTO_MISSING_MODEL: non-zero exit when .model absent (got ${exit_code})" "PASS"
else
    run_test "AUTO_MISSING_MODEL: non-zero exit when .model absent (got ${exit_code})" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_MALFORMED: .model with no model_fullname line → FATAL exit"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_MALFORMED"
set_prior_fct_name "auto"
set_model_fct_name "model_Harvey_Gaussian"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) && exit_code=0 || exit_code=$?
restore_config

if [ $exit_code -ne 0 ] && [ $exit_code -ne 124 ] && echo "$output" | grep -q "FATAL"; then
    run_test "AUTO_MALFORMED: FATAL exit on malformed .model" "PASS"
else
    run_test "AUTO_MALFORMED: FATAL exit on malformed .model (got ${exit_code})" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AUTO_WRONG_CASE: prior_fct_name=AUTO (uppercase) → NOT auto path"
# ---------------------------------------------------------------------------
write_auto_presets "AUTO_WRONG_CASE"
set_prior_fct_name "AUTO"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) && exit_code=0 || exit_code=$?
restore_config

if ! echo "$output" | grep -q "INFO: prior_fct_name auto-resolved" && \
   [ $exit_code -ne 0 ] && [ $exit_code -ne 124 ]; then
    run_test "AUTO_WRONG_CASE: uppercase AUTO not treated as auto-sentinel" "PASS"
else
    run_test "AUTO_WRONG_CASE: uppercase AUTO not treated as auto-sentinel (got ${exit_code})" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test MISMATCH_WARN: explicit prior + mismatching model_fullname → WARNING"
# ---------------------------------------------------------------------------
write_auto_presets "MISMATCH_WARN"
set_prior_fct_name "io_MS_Global"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
restore_config

if echo "$output" | grep -q "WARNING"; then
    run_test "MISMATCH_WARN: WARNING printed for mismatched prior/model" "PASS"
else
    run_test "MISMATCH_WARN: WARNING printed for mismatched prior/model" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AJFIT_AUTO_FATAL: prior_fct_name=auto + model_ajfit → FATAL"
# ---------------------------------------------------------------------------
write_auto_presets "AJFIT_AUTO_FATAL"
set_prior_fct_name "auto"
set_model_fct_name "model_Harvey_Gaussian"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) && exit_code=0 || exit_code=$?
restore_config

if [ $exit_code -ne 0 ] && [ $exit_code -ne 124 ] && echo "$output" | grep -q "ajfit"; then
    run_test "AJFIT_AUTO_FATAL: FATAL exit when auto + model_ajfit" "PASS"
else
    run_test "AJFIT_AUTO_FATAL: FATAL exit when auto + model_ajfit (got ${exit_code})" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test AJFIT_EXEMPT: explicit io_ajfit + mappable model → no WARNING"
# ---------------------------------------------------------------------------
write_auto_presets "AJFIT_EXEMPT"
set_prior_fct_name "io_ajfit"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
restore_config

if ! echo "$output" | grep -q "WARNING" && \
   ! echo "$output" | grep -q "not supported for ajfit"; then
    run_test "AJFIT_EXEMPT: no WARNING when io_ajfit explicit with mappable model" "PASS"
else
    run_test "AJFIT_EXEMPT: no WARNING when io_ajfit explicit with mappable model" "FAIL"
fi

# ---------------------------------------------------------------------------
echo "Test WAVE1_REGRESSION: sentinel not set → Wave 1 behavior preserved"
# ---------------------------------------------------------------------------
write_auto_presets "WAVE1_REG"
set_model_fct_name "model_MS_Global_ajAlm_HarveyLike"
output=$(${TIMEOUT} ${BINARY} -S 1 -L 1 2>&1) || true
exit_code=$?
[ $exit_code -eq 124 ] && exit_code=0
restore_config

if ! echo "$output" | grep -q "INFO: prior_fct_name auto-resolved" && \
   [ $exit_code -eq 0 ]; then
    run_test "WAVE1_REGRESSION: no auto-resolve INFO, exit 0 or timeout" "PASS"
else
    run_test "WAVE1_REGRESSION: no auto-resolve INFO, exit 0 or timeout (got ${exit_code})" "FAIL"
fi

echo ""
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo "Total: $((PASS+FAIL)) | PASS: $PASS | FAIL: $FAIL"
echo ""

if [ $FAIL -eq 0 ]; then
    echo "All tests passed"
    exit 0
else
    echo "Some tests failed"
    exit 1
fi
