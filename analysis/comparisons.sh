#! /bin/bash

SET="80"
TAG="r0100"

DIR1="/exp/mu2e/data/projects/run1a/mmackenz/combine"
DIR2="/exp/mu2e/data/projects/run1a/sophie-mumem"
FIGDIR="figures/compare/${SET}_${TAG}"
mkdir -p ${FIGDIR}

do_comp() {
    local name=$1
    local oname=$2
    local hname=$3
    echo ">>> Comparing $name"
    python compare_inputs.py \
           --file1 "${DIR1}/comp_mumem_${SET}_evt_${TAG}.root" --hist1 workflow/raw/${name} --rate1 1 \
           --file2 "${DIR2}/cutset-${SET}/hists/postcuts_${oname}_hists.root" --hist2 h_mom_${hname} --rate2 1 \
           --label1 "Michael" --label2 "Sophie" \
           --title "${name} Raw Comparison" \
           --output ${FIGDIR}/${name}_compare.png
}

do_comp signal ce ce
do_comp cosmic onspill cosmic
do_comp dio dio dio
do_comp rpc_ext extRPC extRPC
do_comp rpc_int intRPC intRPC
do_comp rmc_ext extRMC0n extRMC0N
