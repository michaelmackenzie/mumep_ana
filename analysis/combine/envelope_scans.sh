#! /bin/bash
# Make envelope likelihood scans and plot the overlay of the scans

Help() {
  echo "Usage: envelope_scans.sh [options] <card>
  Options:
        -h, --help      Show this help message and exit
        -n, --npdfs     Number of PDFs to use for the envelope (default: 2)
        -i, --index     Index parameter name to use for the envelope (default: pdf_index)
        -T, --tag       Tag for the output files (default: envelope)
        -r, --range     Range for the scan (default: 100)
        "
}

if [[ $# -lt 1 ]]; then
  Help
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      Help
      exit 0
      ;;
    -i|--index)
      INDEX="$2"
      shift 2
      ;;
    -T|--tag)
      TAG="$2"
      shift 2
      ;;
    -r|--range)
      RANGE="$2"
      shift 2
      ;;
    -n|--npdfs)
      NPDFS="$2"
      shift 2
      ;;
    *)
      if [[ -z "$CARD" ]]; then
        CARD="$1"
      else
        echo "Unknown argument: $1"
        Help
        exit 1
      fi
      shift
      ;;
  esac
done

if [[ -z "$CARD" ]]; then
  echo "Error: No card specified."
  Help
  exit 1
fi
if [[ ! -f "$CARD" ]]; then
  echo "Error: Card file $CARD does not exist."
  exit 1
fi
if [[ -z "$NPDFS" ]]; then
  NPDFS=2
fi
if [[ -z "$INDEX" ]]; then
  INDEX="pdf_index"
fi
if [[ -z "$TAG" ]]; then
  TAG="envelope"
fi
if [[ -z "$RANGE" ]]; then
  RANGE=100
fi

FIT_ARGS="--cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1"
FIT_ARGS="$FIT_ARGS --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
FIT_ARGS="$FIT_ARGS --X-rtd MINIMIZER_freezeDisassociatedParams"
FIT_ARGS="$FIT_ARGS --X-rtd MINIMIZER_multiMin_hideConstants"
FIT_ARGS="$FIT_ARGS --X-rtd MINIMIZER_multiMin_maskConstraints"
FIT_ARGS="$FIT_ARGS --saveNLL"

for((i=0; i<${NPDFS}; i++)); do
  echo "Running combine for PDF $i"
  combine -M MultiDimFit "$CARD" --algo grid --points 40 --rMin -${RANGE} --rMax ${RANGE} --setParameters ${INDEX}=$i -n ".${TAG}_pdf$i" $FIT_ARGS --freezeParameters ${INDEX}
  FIT_FILE="higgsCombine.${TAG}_pdf$i.MultiDimFit.mH120.root"
  if [[ ! -f "$FIT_FILE" ]]; then
    echo "Error: Fit file $FIT_FILE was not created."
    exit 1
  fi
done

echo "Running the total envelope fit"
combine -M MultiDimFit "$CARD" --algo grid --points 40 --rMin -${RANGE} --rMax ${RANGE} -n ".${TAG}_total" $FIT_ARGS --cminRunAllDiscreteCombinations --freezeParameters ${INDEX}
FIT_FILE="higgsCombine.${TAG}_total.MultiDimFit.mH120.root"
if [[ ! -f "$FIT_FILE" ]]; then
  echo "Error: Fit file $FIT_FILE was not created."
  exit 1
fi

# Plot the results
root.exe -q -b '../tools/plot_envelope.C("", "'${TAG}'")'