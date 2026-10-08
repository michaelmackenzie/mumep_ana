#! /bin/bash
# Estimate the sensitivity of a test selection: re-fill the dataset histograms with the selection,
# build a histogram-based and an envelope model, and fit the envelope model to the Asimov data set
# of the histogram-based model

set -euo pipefail

Help() {
    echo "Estimate the sensitivity of a test selection with the envelope model fit to the MC-template Asimov data"
    echo "Usage: ./selection_test.sh [options] <selection> <tag>"
    echo "  selection                    Cut applied on top of the set selection, using the selection tree"
    echo "                               branch names, e.g. \"trk_t0 > 540\" (\"1\" for no extra cut)"
    echo "  tag                          Test name: histograms go to test_hists/<tag>/, and the models"
    echo "                               use the tags <tag>_hists and <tag>_env"
    echo "  -p, --process <name>         Process name (default: mumep)"
    echo "  -s, --selections <list>      Selection sets, quoted (default: \"45\")"
    echo "  -c, --control-regions <list> Control region set offsets, quoted (default: none), passed to"
    echo "                               full_loop.sh; their sets are re-filled with the same selection"
    echo "      --cl <value>             Confidence level of the limits (default: 0.9)"
    echo "      --rmax <value>           Upper bound of r (default: 1000 for mumep, 100 otherwise)"
    echo "      --loop-args \"<args>\"     Additional full_loop.sh arguments for both models, quoted"
    echo "      --fit-args \"<args>\"      Additional Combine arguments for the limits, quoted"
    echo "      --skip-hists             Use existing test_hists/<tag>/ histograms"
    echo "      --skip-models            Use existing <tag>_hists and <tag>_env models and cards"
    echo "  -h, --help                   Show this help"
    echo ""
    echo "Outputs: test_hists/<tag>/ holds the histograms, logs/, the Combine outputs in combine/,"
    echo "and summary.txt with the limits"
}

PROCESS="mumep"
SELECTIONS="45"
CONTROLREGIONS=""
CL="0.9"
RMAX=""
LOOPARGS=""
FITARGS=""
SKIPHISTS=""
SKIPMODELS=""
CUT=""
TAG=""
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
        -c|--control-regions)
            CONTROLREGIONS="${2:-}"
            shift 2
            ;;
        --cl)
            CL="${2:-}"
            shift 2
            ;;
        --rmax)
            RMAX="${2:-}"
            shift 2
            ;;
        --loop-args)
            LOOPARGS="${2:-}"
            shift 2
            ;;
        --fit-args)
            FITARGS="${2:-}"
            shift 2
            ;;
        --skip-hists)
            SKIPHISTS="1"
            shift
            ;;
        --skip-models)
            SKIPMODELS="1"
            shift
            ;;
        -*)
            echo "Unknown option: $1"
            Help
            exit 1
            ;;
        *)
            if [[ "${CUT}" == "" ]]; then
                CUT="$1"
            elif [[ "${TAG}" == "" ]]; then
                TAG="$1"
            else
                echo "Unknown argument: $1"
                Help
                exit 1
            fi
            shift
            ;;
    esac
done

if [[ "${CUT}" == "" ]] || [[ "${TAG}" == "" ]]; then
    Help
    exit 1
fi
if [[ "${RMAX}" == "" ]]; then
    RMAX="100"
    [[ "${PROCESS}" == "mumep" ]] && RMAX="1000"
fi

cd "$(dirname "$0")"
OUTDIR="test_hists/${TAG}"
LOGDIR="${OUTDIR}/logs"
COMBDIR="${OUTDIR}/combine"
mkdir -p "${LOGDIR}" "${COMBDIR}"

# Sets to re-fill: the selection sets and their control regions
SETS=""
for SELECTION in ${SELECTIONS}; do
    SETS="${SETS} ${SELECTION}"
    for OFFSET in ${CONTROLREGIONS}; do
        SETS="${SETS} $((SELECTION + OFFSET))"
    done
done
SETS="${SETS# }"

# Merged card names, as written by full_loop.sh
merged_card() {
    local card="datacards/combine_total_${PROCESS}"
    for SELECTION in ${SELECTIONS}; do card="${card}_${SELECTION}"; done
    for OFFSET in ${CONTROLREGIONS}; do card="${card}_cr${OFFSET}"; done
    echo "${card}_$1.txt"
}
HISTCARD=$(merged_card "${TAG}_hists")
ENVCARD=$(merged_card "${TAG}_env")

echo ">>> Selection test ${TAG}: \"${CUT}\" on sets ${SETS}"

# 1) Dataset histograms with the test selection
if [[ "${SKIPHISTS}" == "" ]]; then
    echo ">>> Making the test histograms in ${OUTDIR}"
    root.exe -q -b "make_selection_hists.C(\"${SETS}\", \"${CUT//\"/\\\"}\", \"${OUTDIR}\")" > "${LOGDIR}/make_selection_hists.log" 2>&1 \
        || { echo "Histogram creation failed, see ${LOGDIR}/make_selection_hists.log"; exit 1; }
    grep "selected" "${LOGDIR}/make_selection_hists.log" || true
fi

# 2) Histogram-based and envelope models, built in parallel
if [[ "${SKIPMODELS}" == "" ]]; then
    # No systematic fits: the test histograms only have the nominal set histograms
    LOOPOPTS=(--evt-ana --no-sys -p "${PROCESS}" -s "${SELECTIONS}")
    [[ "${CONTROLREGIONS}" != "" ]] && LOOPOPTS+=(-c "${CONTROLREGIONS}")
    echo ">>> Building the models ${TAG}_hists and ${TAG}_env (logs in ${LOGDIR})"
    CONVANA_HIST_PATH="${PWD}/${OUTDIR}" CONVANA_HIST_PDFS=1 CONVANA_USE_ENV=0 \
        ./full_loop.sh "${LOOPOPTS[@]}" -t "${TAG}_hists" ${LOOPARGS} > "${LOGDIR}/model_hists.log" 2>&1 &
    PID_HISTS=$!
    CONVANA_HIST_PATH="${PWD}/${OUTDIR}" CONVANA_HIST_PDFS=0 CONVANA_USE_ENV=1 \
        ./full_loop.sh "${LOOPOPTS[@]}" -t "${TAG}_env" ${LOOPARGS} > "${LOGDIR}/model_env.log" 2>&1 &
    PID_ENV=$!
    STATUS=0
    wait ${PID_HISTS} || { echo "Histogram model failed, see ${LOGDIR}/model_hists.log"; STATUS=1; }
    wait ${PID_ENV}   || { echo "Envelope model failed, see ${LOGDIR}/model_env.log"; STATUS=1; }
    [[ ${STATUS} -ne 0 ]] && exit 1
fi
for CARD in "${HISTCARD}" "${ENVCARD}"; do
    if [[ ! -f "${CARD}" ]]; then
        echo "Card ${CARD} not found"
        exit 1
    fi
done

# 3) Asimov data set of the histogram-based model, and the limits
ARGS="--rMin 0. --rMax ${RMAX} --cl ${CL} --cminDefaultMinimizerStrategy=1 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rAbsAcc 0.001 --rRelAcc 0.001"
ARGS="${ARGS} --cminRunAllDiscreteCombinations --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams"
ARGS="${ARGS} --X-rtd MINIMIZER_multiMin_hideConstants"
if [[ "${FITARGS}" != "" ]]; then
    ARGS="${ARGS} ${FITARGS}"
fi

# Print and run a command, logging its output
run() {
    local log=$1; shift
    echo "$@" | tee "${log}"
    "$@" >> "${log}" 2>&1 || { echo "Command failed, see ${log}"; exit 1; }
}

echo ">>> Generating the Asimov data set of ${HISTCARD}"
run "${LOGDIR}/asimov.log" combine -M GenerateOnly "${HISTCARD}" -t -1 --expectSignal 0 --saveToys -n ".${TAG}_asimov"
ASIMOV="higgsCombine.${TAG}_asimov.GenerateOnly.mH120.123456.root"
mv "${ASIMOV}" "${COMBDIR}/"
ASIMOV="${COMBDIR}/${ASIMOV}"

echo ">>> Fitting the envelope model ${ENVCARD} to the histogram-model Asimov data set"
run "${LOGDIR}/limit_env.log" combine -M AsymptoticLimits -d "${ENVCARD}" ${ARGS} -t -1 --toysFile "${ASIMOV}" -n ".${TAG}_env"
echo ">>> Reference: histogram model on its own Asimov data set"
run "${LOGDIR}/limit_hists.log" combine -M AsymptoticLimits -d "${HISTCARD}" ${ARGS} -t -1 --noFitAsimov -n ".${TAG}_hists"
mv higgsCombine.${TAG}_env.AsymptoticLimits.mH120*.root higgsCombine.${TAG}_hists.AsymptoticLimits.mH120*.root "${COMBDIR}/"

# 4) Summary
limit() { grep "$2" "$1" | tail -n 1 | awk '{print $NF}'; }
SUMMARY="${OUTDIR}/summary.txt"
{
    echo "tag            : ${TAG}"
    echo "selection      : ${CUT}"
    echo "sets           : ${SETS}"
    echo "hists card     : ${HISTCARD}"
    echo "env card       : ${ENVCARD}"
    echo "signal rate    : $(awk '$1 == "rate" && $2 != "group" {print $2}' "${ENVCARD}")"
    echo "observation    : $(grep '^observation' "${HISTCARD}" | awk '{print $2}')"
    echo "${CL} CL limits on r:"
    echo "  envelope fit to the MC-template Asimov : $(limit "${LOGDIR}/limit_env.log" 'Observed Limit')"
    echo "    expected 16% / 50% / 84%             : $(limit "${LOGDIR}/limit_env.log" 'Expected 16.0%') / $(limit "${LOGDIR}/limit_env.log" 'Expected 50.0%') / $(limit "${LOGDIR}/limit_env.log" 'Expected 84.0%')"
    echo "  histogram model on its own Asimov      : $(limit "${LOGDIR}/limit_hists.log" 'Observed Limit')"
} | tee "${SUMMARY}"
echo ">>> Summary written to ${SUMMARY}"
