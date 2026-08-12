#! /bin/bash
# Perform fits to define the process shapes

set -euo pipefail

Help() {
    echo "Create the signal and background model workspaces"
    echo "Usage: ./perform_fits.sh <process, default \"mumem\"> <selection, default 20> <optional tag>"
}

PROCESS=${1:-}
SELECTION=${2:-}
TAG=${3:-}
EVTANA=${4:-}

if [[ "${PROCESS}" == "-h" ]] || [[ "${PROCESS}" == "--help" ]]; then
    Help
    exit 0
fi

if [[ "${PROCESS}" == "" ]]; then
    PROCESS="mumem"
fi
if [[ "${SELECTION}" == "" ]]; then
    SELECTION="20"
fi

run_fit() {
    local macro=$1
    root.exe -q -b "${macro}.C(\"${PROCESS}\", ${SELECTION}, \"${TAG}\")"
}

run_fit signal_fit
run_fit cosmic_fit
run_fit rpc_ext_fit
run_fit rpc_int_fit
if [[ "${EVTANA}" == "" ]]; then
    run_fit pbar_fit
fi
if [[ "${PROCESS}" == "mumem" ]]; then
    run_fit dio_fit
fi
if [[ "${PROCESS}" == "mumep" ]]; then
    run_fit rmc_int_fit
    run_fit rmc_ext_fit
elif [[ "${EVTANA}" != "" ]]; then
    run_fit rmc_ext_fit
fi
run_fit systematics
