#! /bin/bash
# Run the RooFit statistical analysis on a model workspace produced by build_model.C
#
# Usage:
#   ./run_stat_analysis.sh [process] [selection] [tag] [analyses_mask] [ntoys] [mu_inject]
#
# Arguments:
#   process        – "mumem" or "mumep"             (default: mumem)
#   selection      – selection integer               (default: 20)
#   tag            – optional tag passed to build    (default: "")
#   analyses_mask  – hex/dec bitmask of analyses:
#                      0x01 = profile likelihood fit
#                      0x02 = CLs asymptotic 90% UL
#                      0x04 = CLs toy-based 90% UL
#                      0x08 = Feldman-Cousins 90% CI
#                      0x10 = toy MC closure / pull study
#                    (default: 0x3 = fit + CLs asymptotic, fast)
#   ntoys          – number of toys for CLs/FC/pull  (default: 500)
#   mu_inject      – injected signal yield for toys  (default: 0.0)
#
# Examples:
#   ./run_stat_analysis.sh mumem 20 ""  0x3          # fit + CLs asymptotic only (fast)
#   ./run_stat_analysis.sh mumem 20 ""  0x1f 500 0   # all analyses, 500 toys, no signal injected
#   ./run_stat_analysis.sh mumem 20 ""  0x1f 1000 1  # all analyses, 1000 toys, 1-event signal
#   ./run_stat_analysis.sh mumem 20 mds1f 0x3        # use the mds1f-tagged workspace

Help() {
    echo "Run the RooFit statistical analysis"
    echo "Usage: ./run_stat_analysis.sh <process> <selection> <tag> <analyses_mask> <ntoys> <mu_inject>"
    echo "  analyses_mask (hex): 0x01=fit, 0x02=CLs_asym, 0x04=CLs_toys, 0x08=FC, 0x10=toys"
}

if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    Help; exit 0
fi

PROCESS="${1:-mumem}"
SELECTION="${2:-20}"
TAG="${3:-}"
ANALYSES="${4:-0x3}"
NTOYS="${5:-500}"
MU_INJECT="${6:-0.}"

echo "============================================================"
echo " Statistical analysis"
echo "   process   : ${PROCESS}"
echo "   selection : ${SELECTION}"
echo "   tag       : ${TAG:-<none>}"
echo "   analyses  : ${ANALYSES}  (0x01=fit,0x02=CLs_asym,0x04=CLs_toys,0x08=FC,0x10=toys)"
echo "   ntoys     : ${NTOYS}"
echo "   mu_inject : ${MU_INJECT}"
echo "============================================================"

root.exe -q -b "stat_analysis.C(\"${PROCESS}\", ${SELECTION}, \"${TAG}\", ${ANALYSES}, ${NTOYS}, ${MU_INJECT})"
STATUS=$?

echo ""
if [[ ${STATUS} -eq 0 ]]; then
    echo "stat_analysis completed successfully."
else
    echo "stat_analysis finished with status ${STATUS}."
fi
exit ${STATUS}
