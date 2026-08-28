#! /bin/bash
# Compare selection set sensitivities

CARDS=$1
CL=$2
SKIPTOYS=$3
if [[ "${CARDS}" == "" ]]; then
    CARDS="20_r0102 24_r0102 25_r0102 24_25_r0102 60_evt_r0101_0d50bins 20_evt_r0100 24_evt_r0100 25_evt_r0100 24_25_evt_r0100"
fi
if [[ "${CL}" == "" ]]; then
    echo ">>> Using 90% CL CLs"
    CL="0.9"
fi

FREEZE="--freezeParameters allConstrainedNuisances"
ARGS="-t -1 --rMin 0. --rMax 100. --cl ${CL} --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rAbsAcc 0.001 --rRelAcc 0.001"
for CARD in $CARDS; do
    echo "============================================"
    echo "Processing card ${CARD}"
    echo "============================================"
    COMMAND="combine -d datacards/combine_total_mumem_${CARD}.txt ${ARGS}"
    echo ${COMMAND}
    ${COMMAND}
    COMMAND="combine -d datacards/combine_total_mumem_${CARD}.txt ${ARGS} ${FREEZE}"
    echo ${COMMAND}
    ${COMMAND}
    COMMAND="combine -d datacards/combine_total_mumem_${CARD}_cc.txt ${ARGS}"
    echo ${COMMAND}
    ${COMMAND}
    COMMAND="combine -d datacards/combine_total_mumem_${CARD}_cc.txt ${ARGS} ${FREEZE}"
    echo ${COMMAND}
    ${COMMAND}

    if [[ "${SKIPTOYS}" == "" ]]; then
        echo "======> Performing toy-based CLs calculations"
        COMMAND="cls_toys.sh datacards/combine_total_mumem_${CARD}.txt"
        echo ${COMMAND}
        ${COMMAND}
        COMMAND="cls_toys.sh datacards/combine_total_mumem_${CARD}.txt"
        echo ${COMMAND} \"${FREEZE}\"
        ${COMMAND} "${FREEZE}"
    fi
done

echo "Done."
