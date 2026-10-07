#! /bin/bash
# Compare selection set sensitivities

Help() {
    echo "Compare selection set sensitivities"
    echo "Usage: ./selection_comparisons.sh [options] [card ...]"
    echo "  card                   Merged card name(s), datacards/combine_total_<card>.txt"
    echo "                         (default: mumem_75_evt_r0101 mumem_77_78_evt_r0101)"
    echo "  --cl <value>           Confidence level of the CLs limits (default: 0.9)"
    echo "  --single-test          Only the nominal shape limit with all uncertainties, for each card"
    echo "  --skip-toys            Skip the toy-based CLs limits"
    echo "  --fit-args \"<args>\"    Additional Combine arguments for every fit, quoted"
    echo "  -h, --help             Show this help"
}

CARDS=""
CL="0.9"
SINGLETEST=""
SKIPTOYS=""
FITARGS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            Help
            exit 0
            ;;
        --cl)
            CL="${2:-}"
            shift 2
            ;;
        --single-test)
            SINGLETEST="1"
            shift
            ;;
        --skip-toys)
            SKIPTOYS="1"
            shift
            ;;
        --fit-args)
            FITARGS="${2:-}"
            shift 2
            ;;
        -*)
            echo "Unknown option: $1"
            Help
            exit 1
            ;;
        *)
            CARDS="${CARDS} $1"
            shift
            ;;
    esac
done

if [[ "${CARDS}" == "" ]]; then
    CARDS="mumem_75_evt_r0101 mumem_77_78_evt_r0101"
fi
echo ">>> Using ${CL} CL CLs"

RANGE="100"
if [[ "${CARDS}" == *"mumep"* ]]; then
    RANGE="1000"
fi

FREEZE="--freezeParameters allConstrainedNuisances"
ARGS="-t -1 --rMin 0. --rMax ${RANGE} --cl ${CL} --cminDefaultMinimizerStrategy=1 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rAbsAcc 0.001 --rRelAcc 0.001"
ARGS="${ARGS} --cminRunAllDiscreteCombinations --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams"
ARGS="${ARGS} --X-rtd MINIMIZER_multiMin_hideConstants --noFitAsimov"
if [[ "${FITARGS}" != "" ]]; then
    ARGS="${ARGS} ${FITARGS}"
fi

# Print and run a command
run() {
    echo "$@"
    "$@"
}

for CARD in ${CARDS}; do
    echo "============================================"
    echo "Processing card ${CARD}"
    echo "============================================"
    # nominal: the shape card with all uncertainties
    run combine -d datacards/combine_total_${CARD}.txt ${ARGS}
    if [[ "${SINGLETEST}" != "" ]]; then
        continue
    fi
    run combine -d datacards/combine_total_${CARD}.txt ${ARGS} ${FREEZE}
    run combine -d datacards/combine_total_${CARD}_cc.txt ${ARGS}
    run combine -d datacards/combine_total_${CARD}_cc.txt ${ARGS} ${FREEZE}

    if [[ "${SKIPTOYS}" == "" ]]; then
        echo "======> Performing toy-based CLs calculations"
        run "$(dirname "$0")/cls_toys.sh" datacards/combine_total_${CARD}.txt "${FITARGS}"
    fi
done

echo "Done."
