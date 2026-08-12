#! /bin/bash
# Create the search setup: fit backgrounds, create cards, merge categories

set -euo pipefail

Help() {
    echo "Create the search setup: fit backgrounds, create cards, merge categories"
    echo "Usage: ./full_loop.sh [options]"
    echo "  -p, --process <name>         Process name (default: mumem)"
    echo "  -s, --selections <list>      Selection list, quoted (default: \"20\")"
    echo "  -t, --tag <tag>              Optional output tag"
    echo "      --skip-fits              Skip per-selection fit stage"
    echo "      --skip-model             Skip model-building stage"
    echo "      --skip-cards             Skip datacard merge stage"
    echo "      --skip-combine           Skip combine execution stage"
    echo "      --evt-ana                Assume Mu2eEvtAna processing"
    echo "  -h, --help                   Show this help"
    echo ""
    echo "Backward-compatible positional mode is still supported:"
    echo "  ./full_loop.sh <process> <selection list> <tag> <skip fits> <skip model> <skip cards> <skip combine>"
}

PROCESS=""
SELECTIONS=""
TAG=""
SKIPFITS=""
SKIPMODEL=""
SKIPCARDS=""
SKIPCOMBINE=""
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
            --skip-combine)
                SKIPCOMBINE="1"
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
fi

if [[ "${PROCESS}" == "" ]]; then
    PROCESS="mumem"
fi
if [[ "${SELECTIONS}" == "" ]]; then
    SELECTIONS="20"
fi

MERGECARD="datacards/combine_total_${PROCESS}"
declare -a MERGELIST=()
declare -a MERGELISTCC=()
for SELECTION in ${SELECTIONS}; do
    if [[ "${SKIPFITS}" == "" ]]; then
        ./perform_fits.sh "${PROCESS}" "${SELECTION}" "${TAG}" "${EVTANA}"
    fi
    if [[ "${SKIPMODEL}" == "" ]]; then
        root.exe -q -b "build_model.C(\"${PROCESS}\", ${SELECTION}, \"${TAG}\")"
    fi
    CARD="combine_${PROCESS}_${SELECTION}"
    if [[ "${TAG}" != "" ]]; then
        CARD="${CARD}_${TAG}"
    fi
    MERGELIST+=("set_${SELECTION}=datacards/${CARD}.txt")
    MERGELISTCC+=("set_${SELECTION}=datacards/${CARD}_cc.txt")
    MERGECARD="${MERGECARD}_${SELECTION}"
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

if [[ "${SKIPCOMBINE}" == "" ]]; then
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
