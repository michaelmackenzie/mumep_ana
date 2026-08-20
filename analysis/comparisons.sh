#! /bin/bash

SET="80"
TAG="r0100"

DIR1="/exp/mu2e/data/projects/run1a/mmackenz/combine"
DIR2="/exp/mu2e/data/projects/run1a/sophie-mumem"
FIGDIR="figures/compare/${SET}_${TAG}"
JSON="${DIR2}/cutset-${SET}/normalizations.json"
mkdir -p ${FIGDIR}

do_comp() {
    local name=$1
    local oname=$2
    local hname=$3
    local rname=$4
    if [ -f ${JSON} ]; then
        RATE=`grep ${rname} ${JSON} | tail -n 1 | awk '{print $NF}' | sed 's/,//g'`
    else
        RATE="1"
    fi
    echo ">>> Comparing $name: Extracted rate = ${RATE}"
    python compare_inputs.py \
           --file1 "${DIR1}/comp_mumem_${SET}_evt_${TAG}.root" --hist1 workflow/raw/${name} --rate-hist1 workflow/normalized/${name} \
           --file2 "${DIR2}/cutset-${SET}/hists/postcuts_${oname}_hists.root" --hist2 h_mom_${hname} --rate2 ${RATE} \
           --label1 "Michael" --label2 "Sophie" \
           --title "${name} Raw Comparison" \
           --output ${FIGDIR}/${name}_compare.png
}

do_comp signal ce ce CE
do_comp cosmic onspill cosmic cosmics
do_comp dio dio dio DIO
do_comp rpc_ext extRPC extRPC RPCExt
do_comp rpc_int intRPC intRPC RPCInt
do_comp rmc_ext extRMC0n extRMC0N RMCExt0N
