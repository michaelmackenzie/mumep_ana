#! /bin/bash
# Compare selection set sensitivities

CARDS=$1
CL=$2
SKIPTOYS=$3
if [[ "${CARDS}" == "" ]]; then
    CARDS="mumem_75_evt_r0101 mumem_77_78_evt_r0101"
fi
if [[ "${CL}" == "" ]]; then
    echo ">>> Using 90% CL CLs"
    CL="0.9"
fi

RANGE="100"
if [[ "${CARDS}" == *"mumep"* ]]; then
    RANGE="1000"
fi

FREEZE="--freezeParameters allConstrainedNuisances"
ARGS="-t -1 --rMin 0. --rMax ${RANGE} --cl ${CL} --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rAbsAcc 0.001 --rRelAcc 0.001"

for CARD in $CARDS; do
    echo "============================================"
    echo "Processing card ${CARD}"
    echo "============================================"
    COMMAND="combine -d datacards/combine_total_${CARD}.txt ${ARGS}"
    echo ${COMMAND}
    ${COMMAND}
    COMMAND="combine -d datacards/combine_total_${CARD}.txt ${ARGS} ${FREEZE}"
    echo ${COMMAND}
    ${COMMAND}
    COMMAND="combine -d datacards/combine_total_${CARD}_cc.txt ${ARGS}"
    echo ${COMMAND}
    ${COMMAND}
    COMMAND="combine -d datacards/combine_total_${CARD}_cc.txt ${ARGS} ${FREEZE}"
    echo ${COMMAND}
    ${COMMAND}

    if [[ "${SKIPTOYS}" == "" ]]; then
        echo "======> Performing toy-based CLs calculations"
        COMMAND="cls_toys.sh datacards/combine_total_${CARD}.txt"
        echo ${COMMAND}
        ${COMMAND}
        COMMAND="cls_toys.sh datacards/combine_total_${CARD}.txt"
        echo ${COMMAND} \"${FREEZE}\"
        ${COMMAND} "${FREEZE}"
    fi
done

echo "Done."
