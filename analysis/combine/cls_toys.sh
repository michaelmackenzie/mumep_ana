#! /bin/bash
# Evaluate the CLs upper limit using toys

CARD=$1
if [[ "${CARD}" == "" ]]; then
    echo "No card given!"
fi

ARGS="-M HybridNew --LHCmode LHC-limits"
ARGS="${ARGS} -t -1 --rMin 0. --rMax 30. --cl 0.9 --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1"
ARGS="${ARGS} --clsAcc 0.005 --confidenceTolerance 0.05 --interpAcc 0.1 --rRelAcc 0.001 --rAbsAcc 0.01"
ARGS="${ARGS} --plot figures/limit.png"
echo "============================================"
echo "Processing card ${CARD}"
echo "============================================"
COMMAND="combine -d ${CARD} ${ARGS}"
echo ${COMMAND}
${COMMAND}
echo "Done."
