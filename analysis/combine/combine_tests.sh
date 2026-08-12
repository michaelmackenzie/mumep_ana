# Some standard tests/plots using Combine
#! /bin/bash

Help() {
    echo "Perform Combine fits/tests to an input data card"
    echo "Usage: ./combine_tests.sh <data card> <optional tag> <optional skip fits>"
}

DATACARD=$1
TAG=$2
SKIPFITS=$3

if [[ "${DATACARD}" == "" ]] || [[ "${DATACARD}" == "-h" ]] || [[ "${DATACARD}" == "--help" ]]; then
    Help
    exit
fi

if [ ! -f ${DATACARD} ]; then
    echo "No data card ${DATACARD} found!"
    Help
    exit
fi

FIGDIR="figures/combine_tests"
OUTDIR="outputs/combine_tests"
[ ! -d ${FIGDIR} ] && mkdir -p ${FIGDIR}
[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}

if [[ "${TAG}" == "" ]]; then
    TAG=`basename ${DATACARD} | sed 's/.txt//g' | sed 's/combine_//g'`
fi
echo "Using output naming tag ${TAG}"

#----------------------------------------------------------------
# CLs upper limit (Asimov)
#----------------------------------------------------------------

echo "Asimov CLs upper limit"
if [[ "${SKIPFITS}" == "" ]]; then
    combine -d ${DATACARD} -t -1 -n ".${TAG}"
    FILE="higgsCombine.${TAG}.AsymptoticLimits.mH120.root"
    [ -f ${FILE} ] && mv ${FILE} ${OUTDIR}/
fi

#----------------------------------------------------------------
# Example toy fit
#----------------------------------------------------------------

echo "Single fit, no signal"
if [[ "${SKIPFITS}" == "" ]]; then
    combine -M FitDiagnostics  -d ${DATACARD} -t 1 -n ".${TAG}_Single_B" --rMin -1. --rMax 1. --saveShapes --saveWithUncertainties --plot --out ${FIGDIR}
    FILE="${FIGDIR}/fitDiagnostics.${TAG}_Single_B.root"
    [ -f ${FILE} ] && mv ${FILE} ${OUTDIR}/
    rm higgsCombine.${TAG}_Single_B.FitDiagnostics.mH120.123456.root
fi
FILE="${OUTDIR}/fitDiagnostics.${TAG}_Single_B.root"
if [ -f ${FILE} ]; then
    root.exe -q -b "../tools/plot_combine_fit.C(\"${FILE}\", \"${FIGDIR}/${TAG}_Single_B\")"
else
    echo "File ${FILE} not found!"
fi

echo "Single fit, signal"
if [[ "${SKIPFITS}" == "" ]]; then
    combine -M FitDiagnostics  -d ${DATACARD} -t 1 -n ".${TAG}_Single_S" --rMin -10. --rMax 10. --expectSignal 5.0 --saveShapes --plot --out ${FIGDIR}/
    FILE="${FIGDIR}/fitDiagnostics.${TAG}_Single_S.root"
    [ -f ${FILE} ] && mv ${FILE} ${OUTDIR}/
    rm higgsCombine.${TAG}_Single_S.FitDiagnostics.mH120.123456.root
fi
FILE="${OUTDIR}/fitDiagnostics.${TAG}_Single_S.root"
if [ -f ${FILE} ]; then
    root.exe -q -b "../tools/plot_combine_fit.C(\"${FILE}\", \"${FIGDIR}/${TAG}_Single_S\")"
else
    echo "File ${FILE} not found!"
fi

#----------------------------------------------------------------
# Fit pulls/POI distributions
#----------------------------------------------------------------

echo "Fit diagnostics, no signal"
if [[ "${SKIPFITS}" == "" ]]; then
    combine -M FitDiagnostics  -d ${DATACARD} -t 1000 -n ".${TAG}_B" --rMin -1. --rMax 1.
    FILE="fitDiagnostics.${TAG}_B.root"
    [ -f ${FILE} ] && mv ${FILE} ${OUTDIR}/
    rm higgsCombine.${TAG}_B.FitDiagnostics.mH120.123456.root
fi
FILE="${OUTDIR}/fitDiagnostics.${TAG}_B.root"
if [ -f ${FILE} ]; then
    root.exe -q -b "../tools/plot_combine_fits.C(\"${FILE}\", 0., \"${FIGDIR}/${TAG}_B\")"
else
    echo "File ${FILE} not found!"
fi

echo "Fit diagnostics, signal"
if [[ "${SKIPFITS}" == "" ]]; then
    combine -M FitDiagnostics  -d ${DATACARD} -t 1000 -n ".${TAG}_S" --expectSignal 5.0 --rMin -10. --rMax 10.
    FILE="fitDiagnostics.${TAG}_S.root"
    [ -f ${FILE} ] && mv ${FILE} ${OUTDIR}/
    rm higgsCombine.${TAG}_S.FitDiagnostics.mH120.123456.root
fi
FILE="${OUTDIR}/fitDiagnostics.${TAG}_S.root"
if [ -f ${FILE} ]; then
    root.exe -q -b "../tools/plot_combine_fits.C(\"${FILE}\", 5., \"${FIGDIR}/${TAG}_S\")"
else
    echo "File ${FILE} not found!"
fi

echo "Finished!"
