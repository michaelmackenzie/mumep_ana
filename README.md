# Muon to positron conversion analysis

This repo is intended for muon to positron analysis at Mu2e, focused on Run 1A.

## Best practices

In general, we should try to follow these guidelines:
- We should strive for collaborative work, including all who want to participate.
- All work and results should be scripted and repeatable.
- Figures in the paper should be repeatable and created using scripts added to this repository.
- Studies that evaluate important paper inputs, such as systematic uncertainties, should be documented here.
- Inputs to analysis tools (such as histograms and datasets) should be stored in a central location such as: `/exp/mu2e/data/projects/run1a/mumep_ana/`
- Sensitivity results/analysis workspaces should be documented and similarly stored.
- Tools to evaluate inputs or sensitivity should be available to the collaboration, regularly updated, and version controlled to track changes in results.

## Building

Checkout the relevant repositories:
```bash
git clone https://github.com/michaelmackenzie/mumep_ana.git
git clone --branch dev --single-branch https://github.com/michaelmackenzie/Stntuple.git
git clone https://github.com/michaelmackenzie/Mu2eEvtAna.git
git clone https://github.com/Mu2e/Offline.git
git clone https://github.com/Mu2e/mu2e-trig-config.git
git clone https://github.com/Mu2e/Production.git
git clone https://github.com/Mu2e/EventNtuple.git
git clone https://github.com/Mu2e/ArtAnalysis.git
git clone https://github.com/Mu2e/MLTrain.git
git clone https://github.com/michaelmackenzie/grim.git
```

Compile on mu2ebuild02
```bash
ssh mu2ebuild02
mu2einit
# cd /path/to/area
muse setup
time muse build --mu2eCompactPrint --mu2ePyWrap --mu2eCBD -j20
```

## Combine installation:
Clone the packages
```bash
git clone --branch mu2e_dev --single-branch https://github.com/michaelmackenzie/HiggsAnalysis-CombinedLimit.git combine/HiggsAnalysis/CombinedLimit
git clone https://github.com/drbenmorgan/vdt.git combine/vdt
git clone https://gitlab.com/libeigen/eigen.git --depth 1 --branch 3.4.0 ~/local/eigen-3.4.0
```

Build the packages:
```bash
ssh mu2ebuild02
# cd /path/to/combine

# setup python environment
cd HiggsAnalysis/CombinedLimit/
mu2einit
pyenv rootana 2.5.0
source env_standalone_mu2e.sh
cd ../..

# Compile vdt
cd vdt
[ -e CMakeCache.txt ] && rm CMakeCache.txt
cmake -DCMAKE_INSTALL_PREFIX=build  .
make
make install
cd ..

# Compile Combine
cd HiggsAnalysis/CombinedLimit/
make -j20; make
```

## Setting up Combine environment
Combine is not built with muse/Offline, and requires python libraries not in the normal Offline environment.
Instead, pyenv rootana is used with Combine, without setting up a muse/Offline environment.

```bash
cd combine/HiggsAnalysis/CombinedLimit/
mu2einit
pyenv rootana 2.5.0
source env_standalone_mu2e.sh
cd ../..
```

## Repository organization

## Submitting Stntuple jobs

## Processing Stntuple ntuples

## Processing EventNtuple ntuples

EventNtuple inputs can be processed to make histograms/slim trees that are input to the statistical analysis tools.
An example usage of these with Mu2eEvtAna is:
```bash
# See Mu2eEvtAna/scripts/datasets.C to see the naming convention
# Ensure your rootlogon.C includes: 'gSystem->Load("$MUSE_BUILD_DIR/Mu2eEvtAna/lib/libmu2eevtana.so");'
getToken
root.exe -q -b -l 'Mu2eEvtAna/scripts/make_histograms.C(1/*N(processes)*/, "cpos1b1s5r0100", 1, "cnv_ana", 1 /*N(threads)*/, 1e5/*N(events / thread)*/)'
ls -l 
```

## Making initial plots

Preliminary plots can be made using the [Plotter](analysis/plotter) tool.
An example is:
```bash
cd analysis
root.exe -q -b 'make_plots.C(true, {75})'
```

## Making Combine models

The Combine workflow has several steps:
- Fit input histograms and extract PDFs
- Combine the PDFs and rates into a total model
- Create the Combine data cards/workspaces, merging cards for multiple categories if needed
- Running the statistical tests

An example of the complete workflow:
```bash
cd analysis/combine/
# Single category
time ./full_loop.sh -p mumem -s "75" -t evt_r0100 --evt-ana --no-sys
# Two exclusive categories
time ./full_loop.sh -p mumem -s "77 78" -t evt_r0100 --evt-ana --no-sys
```

## Statistical analysis

Using Combine, many statistical tests can be performed.
See https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit for detailed documentation on Combine.

An example evaluation:
```bash
cd analysis/combine/

# Extract approximate Asimov 90% CL CL_s upper limit
CARD="datacards/combine_total_mumem_75_evt_r0100.txt"
ARGS="-t -1 --rMin -100. --rMax 100. --cl 0.9 --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rAbsAcc 0.0001 --rRelAcc 0.001"
combine -d ${CARD} ${ARGS}
```