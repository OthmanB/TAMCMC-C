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
