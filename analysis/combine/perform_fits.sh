#! /bin/bash
# Perform fits to define the process shapes

set -euo pipefail

Help() {
    echo "Create the signal and background model workspaces"
    echo "Usage: ./perform_fits.sh <process, default \"mumem\"> <selection, default 20> <optional tag> <evt ana> <no sys>"
    echo "Environment overrides:"
    echo "  FIT_PDF_TYPE        Global PDF type (default, auto, hist, uniform, poly, polyN, cb, analytic)"
    echo "  FIT_TAIL_MODEL      Global tail model (default, exp, power, convolution, none)"
    echo "  FIT_PDF_TYPE_<COMP> Component PDF override, e.g. FIT_PDF_TYPE_SIGNAL=cb"
    echo "  FIT_TAIL_MODEL_<COMP> Component tail override for signal/dio/rmc_ext/rmc_int"
    echo "    Components: SIGNAL, DIO, COSMIC, RPC_EXT, RPC_INT, PBAR, RMC_EXT, RMC_INT"
    echo "  FIT_SHAPE_SETS      Global comma-separated shape set override"
    echo "  FIT_SHAPE_SETS_<COMP> Component shape set override"
    echo "  FIT_CONTROL_SETS    Global comma-separated control-region set override"
    echo "  FIT_CONTROL_SETS_<COMP> Component control-region set override"
}

PROCESS=${1:-}
SELECTION=${2:-}
TAG=${3:-}
EVTANA=${4:-}
NOSYS=${5:-}

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
    local component=""
    local pdf_default="${FIT_PDF_TYPE:-default}"
    local tail_default="${FIT_TAIL_MODEL:-default}"
    local shape_sets_default="${FIT_SHAPE_SETS:-}"
    local control_sets_default="${FIT_CONTROL_SETS:-}"

    case "${macro}" in
        signal_fit) component="SIGNAL" ;;
        dio_fit)    component="DIO"    ;;
        cosmic_fit) component="COSMIC" ;;
        rpc_ext_fit) component="RPC_EXT" ;;
        rpc_int_fit) component="RPC_INT" ;;
        pbar_fit) component="PBAR" ;;
        rmc_ext_fit) component="RMC_EXT" ;;
        rmc_int_fit) component="RMC_INT" ;;
    esac

    local pdf_var="FIT_PDF_TYPE_${component}"
    local tail_var="FIT_TAIL_MODEL_${component}"
    local shape_sets_var="FIT_SHAPE_SETS_${component}"
    local control_sets_var="FIT_CONTROL_SETS_${component}"
    local pdf="${!pdf_var:-${pdf_default}}"
    local tail="${!tail_var:-${tail_default}}"
    local shape_sets="${!shape_sets_var:-${shape_sets_default}}"
    local control_sets="${!control_sets_var:-${control_sets_default}}"

    local shape_vec="std::vector<int>{}"
    local ctrl_vec="std::vector<int>{}"
    if [[ "${shape_sets}" != "" ]]; then
        shape_vec="std::vector<int>{${shape_sets}}"
    fi
    if [[ "${control_sets}" != "" ]]; then
        ctrl_vec="std::vector<int>{${control_sets}}"
    fi

    if [[ "${macro}" == "signal_fit" ]]; then
        root.exe -q -b "${macro}.C(\"${PROCESS}\", ${SELECTION}, \"${TAG}\", -1, \"${pdf}\", \"${tail}\", ${shape_vec})"
    elif [[ "${macro}" == "dio_fit" ]]; then
        root.exe -q -b "${macro}.C(\"${PROCESS}\", ${SELECTION}, \"${TAG}\", -1, \"${pdf}\", \"${tail}\", ${shape_vec})"
    elif [[ "${macro}" == "cosmic_fit" ]]; then
        root.exe -q -b "${macro}.C(\"${PROCESS}\", ${SELECTION}, \"${TAG}\", -1, \"${pdf}\", ${ctrl_vec})"
    elif [[ "${macro}" == "rmc_ext_fit" ]] || [[ "${macro}" == "rmc_int_fit" ]]; then
        root.exe -q -b "${macro}.C(\"${PROCESS}\", ${SELECTION}, \"${TAG}\", \"${pdf}\", \"${tail}\", ${shape_vec}, ${ctrl_vec})"
    elif [[ "${macro}" == "rpc_ext_fit" ]] || [[ "${macro}" == "rpc_int_fit" ]] || [[ "${macro}" == "pbar_fit" ]]; then
        root.exe -q -b "${macro}.C(\"${PROCESS}\", ${SELECTION}, \"${TAG}\", \"${pdf}\", ${shape_vec}, ${ctrl_vec})"
    else
        root.exe -q -b "${macro}.C(\"${PROCESS}\", ${SELECTION}, \"${TAG}\")"
    fi
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
if [[ "${NOSYS}" == "" ]]; then
    run_fit systematics
else
    echo "Skipping systematics fit (--no-sys enabled)"
fi
