# Muon to positron conversion analysis

This repo is intended for muon to positron analysis at Mu2e, focused on Run 1A.


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
git clone https://github.com/michaelmackenzie/Mu2eEvtAna.git
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
cd HiggsAnalysis/CombinedLimit/
mu2einit
pyenv rootana 2.5.0
source env_standalone_mu2e.sh
cd ../..
```

## Repository organization

## Submitting Stntuple jobs

## Processing Stntuple ntuples

## Processing EventNtuple ntuples

## Making initial plots


## Making Combine models

## Statistical analysis