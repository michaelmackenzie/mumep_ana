#! /bin/bash
# Generate toys from one card and then fit with another

Help() {
  echo "Usage: combine_bias_test.sh [options] <generate_card> <fit_card>
  Options:
    -h, --help      Show this help message and exit
    -F, --fit-args  Additional arguments to pass to the combine tool when fitting the toys
    -G, --gen-args  Additional arguments to pass to the combine tool when generating the toys
    -t, --toys      Number of toys to generate (default: 1000)
    -T, --tag       Tag for the output files (default: bias_test)
    -s, --seed      Random seed for toy generation (default: 123456)
    "
}

if [[ $# -lt 2 ]]; then
  Help
  exit 1
fi

FIT_ARGS=""
GEN_ARGS=""
TOYS=1000
TAG="bias_test"
SEED=123456
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      Help
      exit 0
      ;;
    -F|--fit-args)
      FIT_ARGS="$2"
      shift 2
      ;;
    -G|--gen-args)
      GEN_ARGS="$2"
      shift 2
      ;;
    -t|--toys)
      TOYS="$2"
      shift 2
      ;;
    -T|--tag)
      TAG="$2"
      shift 2
      ;;
    -s|--seed)
      SEED="$2"
      shift 2
      ;;
    *)
      if [[ -z "$CARD1" ]]; then
        CARD1="$1"
      elif [[ -z "$CARD2" ]]; then
        CARD2="$1"
      else
        echo "Unknown argument: $1"
        Help
        exit 1
      fi
      shift
      ;;
  esac
done

if [[ ! -f "$CARD1" ]]; then
  echo "Error: Card $CARD1 does not exist."
  exit 1
fi
if [[ ! -f "$CARD2" ]]; then
  echo "Error: Card $CARD2 does not exist."
  exit 1
fi

echo "Generating $TOYS toys from $CARD1 and fitting with $CARD2 and saving results with tag $TAG"

# Generate the toys
combine -M GenerateOnly "$CARD1" -t "$TOYS" $GEN_ARGS -n ".$TAG" -s "$SEED" --saveToys
TOYFILE="higgsCombine.$TAG.GenerateOnly.mH120.$SEED.root"
if [[ ! -f "$TOYFILE" ]]; then
  echo "Error: Toy file $TOYFILE was not created."
  exit 1
fi

echo "Generated toys saved in $TOYFILE"
FIT_ARGS="$FIT_ARGS --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rMin -500. --rMax 500."
FIT_ARGS="$FIT_ARGS --cminRunAllDiscreteCombinations --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
FIT_ARGS="$FIT_ARGS --X-rtd MINIMIZER_freezeDisassociatedParams"
FIT_ARGS="$FIT_ARGS --X-rtd MINIMIZER_multiMin_hideConstants"
combine -d "$CARD2" -M FitDiagnostics --toysFile="${TOYFILE}" $FIT_ARGS -n ".$TAG" -t "$TOYS" -s "$SEED" --saveWithUncertainties --saveNormalizations

FIT_RESULT_FILE="fitDiagnostics.$TAG.root"
if [[ ! -f "$FIT_RESULT_FILE" ]]; then
  echo "Error: Fit result file $FIT_RESULT_FILE was not created."
  exit 1
fi

echo "Fit results saved in $FIT_RESULT_FILE"
root.exe -q -b -l '../tools/plot_combine_fits.C("'$FIT_RESULT_FILE'", 0., "'$TAG'")'