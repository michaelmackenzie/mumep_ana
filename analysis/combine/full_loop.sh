#! /bin/bash
# Create the search setup: fit backgrounds, create cards, merge categories

set -euo pipefail

Help() {
    echo "Create the search setup: fit backgrounds, create cards, merge categories"
    echo "Usage: ./full_loop.sh [options]"
    echo "  -p, --process <name>         Process name (default: mumem)"
    echo "  -s, --selections <list>      Selection list, quoted (default: \"20\")"
    echo "  -t, --tag <tag>              Optional output tag"
    echo "  -c, --control-regions <list> Control region set offsets added as extra channels to each"
    echo "                               selection, quoted (default: none). Offsets are sums of the flags"
    echo "                               1000 (inverted CRV veto), 2000 (inverted time cut, early time"
    echo "                               for set 75), and 4000 (charge-inverted), e.g. \"1000 2000 4000\""
    echo "      --mc-stat                Add MC statistical yield uncertainties to the cards, for"
    echo "                               processes whose PDF is the MC histogram (not fit or smoothed)"
    echo "      --skip-fits              Skip per-selection fit stage"
    echo "      --skip-model             Skip model-building stage"
    echo "      --skip-cards             Skip datacard merge stage"
    echo "      --test-combine           Test combine execution stage"
    echo "      --no-sys                 Skip fitting systematics"
    echo "      --evt-ana                Assume Mu2eEvtAna processing"
    echo "  -h, --help                   Show this help"
    echo ""
    echo "Environment fit controls forwarded to perform_fits.sh:"
    echo "  FIT_PDF_TYPE, FIT_TAIL_MODEL, FIT_SHAPE_SETS, FIT_CONTROL_SETS"
    echo "  FIT_PDF_TYPE_<COMP>, FIT_TAIL_MODEL_<COMP>, FIT_SHAPE_SETS_<COMP>, FIT_CONTROL_SETS_<COMP>"
    echo "    Components: SIGNAL DIO COSMIC RPC_EXT RPC_INT PBAR RMC_EXT RMC_INT"
    echo "    mumep splits RMC by neutron knockout: RMC_EXT_0N RMC_EXT_1N RMC_INT_0N RMC_INT_1N"
    echo ""
    echo "Backward-compatible positional mode is still supported:"
    echo "  ./full_loop.sh <process> <selection list> <tag> <skip fits> <skip model> <skip cards> <skip combine> <no sys>"
}

PROCESS=""
SELECTIONS=""
TAG=""
CONTROLREGIONS=""
MCSTAT="false"
SKIPFITS=""
SKIPMODEL=""
SKIPCARDS=""
TESTCOMBINE=""
NOSYS=""
EVTANA=""

# Parse flags first. If no leading flag is found, fall back to legacy positional args.
if [[ $# -gt 0 ]] && [[ "${1}" == -* ]]; then
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                Help
                exit 0
                ;;
            -p|--process)
                PROCESS="${2:-}"
                shift 2
                ;;
            -s|--selections)
                SELECTIONS="${2:-}"
                shift 2
                ;;
            -t|--tag)
                TAG="${2:-}"
                shift 2
                ;;
            -c|--control-regions)
                CONTROLREGIONS="${2:-}"
                shift 2
                ;;
            --mc-stat)
                MCSTAT="true"
                shift
                ;;
            --skip-fits)
                SKIPFITS="1"
                shift
                ;;
            --skip-model)
                SKIPMODEL="1"
                shift
                ;;
            --skip-cards)
                SKIPCARDS="1"
                shift
                ;;
            --test-combine)
                TESTCOMBINE="1"
                shift
                ;;
            --no-sys)
                NOSYS="1"
                shift
                ;;
            --evt-ana)
                EVTANA="1"
                shift
                ;;
            --)
                shift
                break
                ;;
            *)
                echo "Unknown option: $1"
                Help
                exit 1
                ;;
        esac
    done
else
    PROCESS=${1:-}
    SELECTIONS=${2:-}
    TAG=${3:-}
    SKIPFITS=${4:-}
    SKIPMODEL=${5:-}
    SKIPCARDS=${6:-}
    SKIPCOMBINE=${7:-}
    NOSYS=${8:-}
fi

if [[ "${PROCESS}" == "" ]]; then
    PROCESS="mumem"
fi
if [[ "${SELECTIONS}" == "" ]]; then
    SELECTIONS="20"
fi

for OFFSET in ${CONTROLREGIONS}; do
    if ! [[ "${OFFSET}" =~ ^[1-7]000$ ]]; then
        echo "Control region offset ${OFFSET} is not a sum of the flags 1000, 2000, and 4000"
        exit 1
    fi
done

# Fit the components and build the model and cards of one selection set
build_set() {
    local set=$1
    if [[ "${SKIPFITS}" == "" ]]; then
        [[ "${FIT_PDF_TYPE:-}" != "" ]] && export FIT_PDF_TYPE
        [[ "${FIT_TAIL_MODEL:-}" != "" ]] && export FIT_TAIL_MODEL
        [[ "${FIT_SHAPE_SETS:-}" != "" ]] && export FIT_SHAPE_SETS
        [[ "${FIT_CONTROL_SETS:-}" != "" ]] && export FIT_CONTROL_SETS
        for comp in SIGNAL DIO COSMIC RPC_EXT RPC_INT PBAR RMC_EXT RMC_INT \
                    RMC_EXT_0N RMC_EXT_1N RMC_INT_0N RMC_INT_1N; do
            var="FIT_PDF_TYPE_${comp}";      [[ "${!var:-}" != "" ]] && export "$var"
            var="FIT_TAIL_MODEL_${comp}";    [[ "${!var:-}" != "" ]] && export "$var"
            var="FIT_SHAPE_SETS_${comp}";    [[ "${!var:-}" != "" ]] && export "$var"
            var="FIT_CONTROL_SETS_${comp}";  [[ "${!var:-}" != "" ]] && export "$var"
        done
        ./perform_fits.sh "${PROCESS}" "${set}" "${TAG}" "${EVTANA}" "${NOSYS}"
    fi
    if [[ "${SKIPMODEL}" == "" ]]; then
        root.exe -q -b "build_model.C(\"${PROCESS}\", ${set}, \"${TAG}\", ${MCSTAT})"
    fi
}

# Card of one selection set, as written by build_model
card_name() {
    local card="combine_${PROCESS}_$1"
    if [[ "${TAG}" != "" ]]; then
        card="${card}_${TAG}"
    fi
    echo "datacards/${card}"
}

MERGECARD="datacards/combine_total_${PROCESS}"
declare -a MERGELIST=()
declare -a MERGELISTCC=()
for SELECTION in ${SELECTIONS}; do
    build_set "${SELECTION}"
    CARD=$(card_name "${SELECTION}")
    MERGELIST+=("set_${SELECTION}=${CARD}.txt")
    MERGELISTCC+=("set_${SELECTION}=${CARD}_cc.txt")
    MERGECARD="${MERGECARD}_${SELECTION}"
    # Control regions enter as extra channels, sharing the process names and nuisances
    for OFFSET in ${CONTROLREGIONS}; do
        CRSET=$((SELECTION + OFFSET))
        build_set "${CRSET}"
        CARD=$(card_name "${CRSET}")
        MERGELIST+=("cr${OFFSET}_${SELECTION}=${CARD}.txt")
        MERGELISTCC+=("cr${OFFSET}_${SELECTION}=${CARD}_cc.txt")
    done
done
for OFFSET in ${CONTROLREGIONS}; do
    MERGECARD="${MERGECARD}_cr${OFFSET}"
done
if [[ "${TAG}" != "" ]]; then
    MERGECARD="${MERGECARD}_${TAG}"
fi
MERGECARDCC="${MERGECARD}_cc.txt"
MERGECARD="${MERGECARD}.txt"

if [[ "${SKIPCARDS}" == "" ]]; then
    echo "Merging combine cards:"
    echo "combineCards.py ${MERGELIST[*]} >| ${MERGECARD}"
    head -n 7 `echo ${MERGELIST[0]} | cut -d '=' -f 2` >| "${MERGECARD}"
    echo "" >> "${MERGECARD}"
    combineCards.py "${MERGELIST[@]}" >> "${MERGECARD}"
    sed -i 's|datacards/||g' "${MERGECARD}"
    echo "combineCards.py ${MERGELISTCC[*]} >| ${MERGECARDCC}"
    head -n 8 `echo ${MERGELISTCC[0]} | cut -d '=' -f 2` >| "${MERGECARDCC}"
    echo "" >> "${MERGECARDCC}"
    combineCards.py "${MERGELISTCC[@]}" >> "${MERGECARDCC}"
fi

if [[ "${TESTCOMBINE}" != "" ]]; then
    echo "Running combine:"
    ADDON=""
    if [[ "${TAG}" != "" ]]; then
        ADDON="-n .${TAG}"
        TAG2="_${TAG}"
    fi

    echo "Running the CLs limits"
    COMMAND="combine -d ${MERGECARD} --rMin -50. --rMax 50. ${ADDON} -t -1 --cl 0.9 --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1"
    echo ${COMMAND}
    ${COMMAND}

    echo "Running a single fit:"
    FITFILE="fitDiagnosticsTest.root"
    if [[ "${TAG}" != "" ]]; then
        FITFILE="fitDiagnostics.${TAG}.root"
    fi
    COMMAND="combine -d ${MERGECARD} -M FitDiagnostics --rMin -50. --rMax 50. --saveShapes --saveWithUncertainties  --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 ${ADDON}"
    echo ${COMMAND}
    ${COMMAND}
    root.exe -q -b '../tools/plot_combine_fit.C("'${FITFILE}'", "figures/'${PROCESS}${TAG2}'", true)'
fi
