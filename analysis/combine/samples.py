"""Mu2eEvtAna (EventNtuple) inputs of the mu- --> e-/e+ search: datasets, rates, and exposures.

This is the single dataset list of the analysis: the Combine models (make_samples) and the ROOT plotting and
selection macros (export_to_root, through ../run_macro.py) both use it, with the exposures (Exposure) and the
histogram directories (hist_paths).

    from samples import Exposure, make_samples
    exposure = Exposure("run1a", two_bb=False)
    samples = make_samples(exposure)
    h = samples["dio"].hist("trk_20/obs")        # TH1D at the expected yield for the exposure

The shared constants come from Mu2eEvtAna's plotter/physics/Mu2ePhysics.C through mu2ecombine.mu2e_physics.
"""

import importlib.util
import os
import sys
from dataclasses import dataclass


def _mu2ecombine_dir():
    """The directory holding the mu2ecombine package, from PYTHONPATH, $MU2ECOMBINE_DIR, or a Mu2eEvtAna clone next
    to this repository (the package is not imported)."""
    spec = importlib.util.find_spec("mu2ecombine")
    if spec is not None and spec.submodule_search_locations:
        return os.path.dirname(list(spec.submodule_search_locations)[0])
    here = os.path.dirname(os.path.realpath(__file__))
    for path in (os.environ.get("MU2ECOMBINE_DIR", ""),
                 os.path.normpath(os.path.join(here, "..", "..", "..", "Mu2eEvtAna", "combine"))):
        if path and os.path.isdir(os.path.join(path, "mu2ecombine")):
            return path
    raise ImportError("mu2ecombine not found: add <Mu2eEvtAna>/combine to PYTHONPATH or set MU2ECOMBINE_DIR")


def _import_mu2ecombine():
    """The full mu2ecombine package (with its fitting dependencies, e.g. scipy): only make_samples needs it."""
    path = _mu2ecombine_dir()
    if path not in sys.path:
        sys.path.insert(0, path)
    import mu2ecombine
    return mu2ecombine


def _load_physics():
    """mu2ecombine.mu2e_physics (Mu2eEvtAna's Mu2ePhysics.C constants). Its module needs only ROOT, so it is loaded on
    its own when the package isn't already imported: the plotting and selection macros (run_macro.py) don't need the
    package's fitting dependencies."""
    if "mu2ecombine" in sys.modules:
        return sys.modules["mu2ecombine"].mu2e_physics
    path = os.path.join(_mu2ecombine_dir(), "mu2ecombine", "physics.py")
    spec = importlib.util.spec_from_file_location("_mu2ecombine_physics", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.mu2e_physics


phys = _load_physics()

DEFAULT_HIST_PATH = "/exp/mu2e/data/projects/run1a/mumep_ana/histograms"
DEFAULT_REFERENCE = 1e-15   # branching fraction the signal is normalized to: r is in units of it
REDUCED_TRACKER   = True
SIGNALS = ("mumem", "mumep")

# Simulated samples, per tracker configuration (True: reduced tracker for Run 1, False: nominal tracker)
# name: dataset, N(generated) (the generated livetime in s for the cosmic rays), events in the full dataset, and
#       the rate key of rates_per_pot
COSMIC_LIVETIME = 4437713.0 * 5.72 / 4.438   # generated livetime x stored-livetime bug correction
DATASETS = {
    True: {
        "mumem":      ("cele1b1s5r0102", 1e7            , 4126584, "signal"),
        "mumep":      ("cpos1b1s5r0102", 1e7            , 3214876, "signal"),
        "cosmic":     ("cry4ab1s5r0102", COSMIC_LIVETIME, 4141826, "cosmic"),
        "dio":        ("dio00b1s5r0102", 25e6           , 9323791, "dio_95"),
        "rmc_ext_0n": ("rmce0b1s5r0102", 7e9            , 4962417, "rmc_0n_80"),
        "rmc_ext_1n": ("rmce1b1s5r0102", 7e9            , 2972496, "rmc_1n_80"),
        "rmc_int_0n": ("rmci0b1s5r0102", 5e7            , 1211523, "rmc_int_0n_80"),
        "rmc_int_1n": ("rmci1b1s5r0102", 5e7            ,  507998, "rmc_int_1n_80"),
        "rpc_ext":    ("rpce1b1s5r0102", 5e9            ,  477976, "rpc"),
        "rpc_int":    ("rpci1b1s5r0102", 1.25e8         , 1910423, "rpc_int"),
        "pbar":       ("pbar1b1s5r0102", 3e7            , 6306269, "pbar"),
    },
    False: {
        "mumem":      ("cele1b1s5r0101", 1e7            , 4141125, "signal"),
        "mumep":      ("cpos1b1s5r0101", 1e7            , 3235878, "signal"),
        "cosmic":     ("cry4ab1s5r0101", COSMIC_LIVETIME, 4155435, "cosmic"),
        "dio":        ("dio00b1s5r0101", 25e6           , 9368976, "dio_95"),
        "rmc_ext_0n": ("rmce0b1s5r0101", 7e9            , 4967393, "rmc_0n_80"),
        "rmc_ext_1n": ("rmce1b1s5r0101", 7e9            , 2974188, "rmc_1n_80"),
        "rmc_int_0n": ("rmci0b1s5r0101", 5e7            , 1229006, "rmc_int_0n_80"),
        "rmc_int_1n": ("rmci1b1s5r0101", 5e7            ,  517029, "rmc_int_1n_80"),
        "rpc_ext":    ("rpce1b1s5r0101", 5e9            ,  458818, "rpc"),
        "rpc_int":    ("rpci1b1s5r0101", 1.25e8         , 1899806, "rpc_int"),
        "pbar":       ("pbar1b1s5r0101", 3e7            , 6461314, "pbar"),
    },
}
# mu- --> e- keeps a single component per RMC source, taken from the 0 neutron knockout sample
ALIASES = {"rmc_ext": "rmc_ext_0n", "rmc_int": "rmc_int_0n"}
# MDS ensembles (1BB, unweighted): dataset and events in the full dataset
DATA = {
    "data_mds3c": ("mds3cb1s5r0100", 4921433),
    "data_mds3d": ("mds3db0s5r0100", 9837962),
}

# Stopping-target designs: OnSpill (unmixed) signal and DIO samples that replace the nominal Mix1BB mu- --> e- signal
# and DIO, with the target's stopped muons per POT. The other processes keep their nominal samples; the RMC follows the
# stops per POT, the RPC (pion stops) does not. All are MDC2025ba CeEndpoint and DIOtail95 samples; "nominal" is the
# nominal target, with the nominal stops per POT.
# Stops per POT = MuBeamCat events per POT (SimEfficiencies2, Sim_best v1_1, run 1430) x mu- stops per resampled
# MuBeamCat event (MuminusStopsCat events / generated events). The nominal TargetStops (MDC2025ac) kept 1 in 1000 stops,
# the alternate ones (from the same MuBeamCat files) every stop.
TARGET_HIST_PATH = "/exp/mu2e/data/projects/run1a/mumep_ana/histograms_alt_target_2026-10-31"
MUBEAM_PER_POT = 213816 / 1e8
TARGETS = {
    # name: signal and DIO (dataset, N(generated), events in the full dataset), mu- stops per resampled MuBeamCat event
    "nominal":        (("cele0b0s5r0101", 1e7, 4247008), ("dio00b0s5r0101", 8e6, 2995328), 1435092 * 1000 / 4e9),
    "st37plus6up6dn": (("cele0b0s5r0110", 1e7, 4002040), ("dio00b0s5r0110", 8e6, 2764431), 5856806 / 1e7),
    "st49hole1to5":   (("cele0b0s5r0120", 1e7, 3967395), ("dio00b0s5r0120", 8e6, 2717026), 6527950 / 1e7),
    "st49hole9to41":  (("cele0b0s5r0130", 1e7, 3961096), ("dio00b0s5r0130", 8e6, 2700249), 5826669 / 1e7),
}


class Exposure:
    """N(POT), livetime, and stops per POT of a data set.

    name:   run1a (28 days of running), run2, mds3c, or mds3d (MDS ensembles, 1BB)
    two_bb: the Run 1A/2 exposure at the 2BB duty cycle and intensity, using the inputs re-weighted to 2BB
            (Mu2eEvtAna histogram mode 2, ConvAna.*.m2.root)
    target: a stopping-target design of TARGETS (Run 1A at 1BB only): its stops per POT, and its OnSpill signal and
            DIO samples in make_samples
    """

    NAMES = ("run1a", "run2", "mds3c", "mds3d")

    def __init__(self, name="run1a", two_bb=False, target=None):
        if name not in self.NAMES:
            raise ValueError(f"unknown exposure '{name}' (available: {', '.join(self.NAMES)})")
        if two_bb and name.startswith("mds"):
            raise ValueError(f"the {name} ensemble is 1BB data, it has no 2BB exposure")
        if target is not None:
            if target not in TARGETS:
                raise ValueError(f"unknown target '{target}' (available: {', '.join(TARGETS)})")
            if name != "run1a" or two_bb:
                raise ValueError("the stopping-target samples are only defined for the Run 1A 1BB exposure")
        self.name, self.two_bb, self.target = name, two_bb, target
        self.nmuons_per_pot = phys.run1a.nmuons_per_pot if target is None else MUBEAM_PER_POT * TARGETS[target][2]
        self.reference = DEFAULT_REFERENCE
        if name == "run1a":
            # the same wall time at the 2BB duty cycle and N(POT) rate when using the 2BB inputs
            duty_cycle = phys.duty_cycle_2bb if two_bb else phys.duty_cycle_1bb
            self.livetime = phys.run1a.livetime / phys.duty_cycle_1bb * duty_cycle
            self.npot = self.livetime * (phys.npot_rate_2bb if two_bb else phys.npot_rate_1bb)
        elif name == "run2":
            self.npot, self.livetime = phys.run2.npot, phys.run2.livetime
            self.nmuons_per_pot = phys.run2.nmuons_per_pot
        else:
            # MDS3c/d EventNtuple ensembles; MDS3c is corrected for the events of its earlier version
            self.livetime = 4401645.8 * (4921433. / 4863968. if name == "mds3c" else 1.)
            self.npot = self.livetime * phys.npot_rate_1bb
            self.reference = 1e-13
        self.nmuons = self.npot * self.nmuons_per_pot

    @property
    def hist_mode(self):
        return 2 if self.two_bb else 1

    def __repr__(self):
        return (f"Exposure({self.name}{', 2BB' if self.two_bb else ''}"
                f"{f', target {self.target}' if self.target else ''}: N(POT) = {self.npot:.3e}, "
                f"livetime = {self.livetime:.3e} s, N(muons) = {self.nmuons:.3e})")


def rates_per_pot(nmuons_per_pot):
    """{rate key: rate per POT} of the DATASETS rate keys at the given stopped muons per POT. The signal rate is per
    unit branching fraction (R_mue is relative to muon capture); the cosmic ray "rate" is 1 (per second of livetime)."""
    f_cap = phys.muon_capture_fraction_al
    rate_rmc = nmuons_per_pot * f_cap * phys.rmc_br_above_57
    rate_rmc_0n_80 = rate_rmc * phys.rmc_ps_0n_r_above_57 / phys.rmc_ps_0n_frac_above_57 * phys.rmc_ps_0n_frac_above_80
    rate_rmc_1n_80 = rate_rmc * phys.rmc_ps_1n_r_above_57 / phys.rmc_ps_1n_frac_above_57 * phys.rmc_ps_1n_frac_above_80
    rate_rpc = phys.run1a.pion_stops_per_pot * phys.rpc_br * phys.run1a.pion_survive_frac
    return {
        "signal":        f_cap * nmuons_per_pot,
        "cosmic":        1.,
        "dio_95":        (1. - f_cap) * nmuons_per_pot * phys.dio_frac_above_95,
        "rmc_0n_80":     rate_rmc_0n_80,
        "rmc_1n_80":     rate_rmc_1n_80,
        "rmc_int_0n_80": rate_rmc_0n_80 * phys.rmc_internal_ratio,
        "rmc_int_1n_80": rate_rmc_1n_80 * phys.rmc_internal_ratio,
        "rpc":           rate_rpc,
        "rpc_int":       rate_rpc * phys.rpc_internal_ratio,
        "pbar":          phys.run1a.pbar_stops_per_pot,
    }


@dataclass
class DatasetInfo:
    """One input of the analysis, independent of the exposure's N(POT) and livetime.

    dataset: dataset name of the ConvAna histogram file (ConvAna.cnv_ana.<dataset>.m<mode>.root)
    ngen:    N(generated events), or the generated livetime (s) for the cosmic rays
    nevents: events in the full dataset (for missing-file corrections), None if unknown
    rate:    rate per POT (per unit branching fraction for the signal), or 1 for the cosmic rays and data
    path:    histogram directory
    kind:    "signal", "background", or "data"
    per_livetime: the yield scales with the livetime (cosmic rays) instead of N(POT)
    """
    dataset: str
    ngen: float
    nevents: int
    rate: float
    path: str
    kind: str = "background"
    per_livetime: bool = False

    def norm(self):
        """Expected yield per generated event per POT (per s of livetime for the cosmic rays), as DatasetInfo_t::norm."""
        return self.rate / self.ngen


def hist_paths(hist_path=None, target_hist_path=None):
    """The histogram directory (default: CONVANA_HIST_PATH, or DEFAULT_HIST_PATH) and the stopping-target histogram
    directory (default: CONVANA_TARGET_HIST_PATH, or TARGET_HIST_PATH)."""
    return (hist_path or os.environ.get("CONVANA_HIST_PATH") or DEFAULT_HIST_PATH,
            target_hist_path or os.environ.get("CONVANA_TARGET_HIST_PATH") or TARGET_HIST_PATH)


def dataset_info(exposure, reduced_tracker=None, hist_path=None, target_hist_path=None):
    """{name: DatasetInfo} of the signal, background, and data inputs at the exposure's stops per POT. With
    exposure.target, the mu- --> e- signal and the DIO are the target's samples. reduced_tracker: default
    REDUCED_TRACKER. hist_path: default CONVANA_HIST_PATH, or DEFAULT_HIST_PATH. target_hist_path: default
    CONVANA_TARGET_HIST_PATH, or TARGET_HIST_PATH."""
    reduced_tracker = REDUCED_TRACKER if reduced_tracker is None else reduced_tracker
    hist_path, target_path = hist_paths(hist_path, target_hist_path)
    rates = rates_per_pot(exposure.nmuons_per_pot)
    infos = {}
    for name, (dataset, ngen, nevents, rate) in DATASETS[reduced_tracker].items():
        infos[name] = DatasetInfo(dataset, ngen, nevents, rates[rate], hist_path,
                                  "signal" if name in SIGNALS else "background", per_livetime=(rate == "cosmic"))
    if exposure.target:
        (signal, signal_gen, signal_events), (dio, dio_gen, dio_events), _ = TARGETS[exposure.target]
        infos["mumem"] = DatasetInfo(signal, signal_gen, signal_events, rates["signal"], target_path, "signal")
        infos["dio"] = DatasetInfo(dio, dio_gen, dio_events, rates["dio_95"], target_path)
    for alias, name in ALIASES.items():
        infos[alias] = infos[name]
    for name, (dataset, nevents) in DATA.items():
        infos[name] = DatasetInfo(dataset, 1., nevents, 1., hist_path, "data")
    return infos


def make_samples(exposure, hist_path=None, target_hist_path=None, reduced_tracker=None):
    """{name: Sample} of the signal, background, and data inputs (dataset_info), normalized to the exposure. The
    signal is normalized to exposure.reference. The data are unweighted 1BB inputs."""
    Sample = _import_mu2ecombine().Sample
    samples = {}
    by_info = {}  # aliases share one Sample
    for name, info in dataset_info(exposure, reduced_tracker, hist_path, target_hist_path).items():
        if id(info) in by_info:
            samples[name] = by_info[id(info)]
            continue
        if info.kind == "data":
            weight, expected_events, mode = 1., None, 1   # no missing-file scaling of the data in the fits
        else:
            weight = (exposure.livetime if info.per_livetime else exposure.npot) * info.norm()
            if info.kind == "signal":
                weight *= exposure.reference
            expected_events, mode = info.nevents, exposure.hist_mode
        path = os.path.join(info.path, f"ConvAna.cnv_ana.{info.dataset}.m{mode}.root")
        samples[name] = by_info[id(info)] = Sample(path, weight=weight, expected_events=expected_events,
                                                   name=info.dataset)
    return samples


def export_to_root(exposure, reduced_tracker=None, hist_path=None, target_hist_path=None):
    """Set up the ROOT analysis macros (loaded beforehand) for the exposure: the exposure of physics.C, the histogram
    directories of defaults.C (hist_paths), and, if datasets.C is loaded, its datasets_ map (dataset_info). The macros
    apply N(POT) or the livetime, and the signal branching fraction, themselves; their init_physics(tag) checks that the
    tag names this exposure. Run the macros through ../run_macro.py, which calls this."""
    import ROOT
    if not hasattr(ROOT, "init_physics"):
        raise RuntimeError("export_to_root: load physics.C (or a macro including it) into ROOT first")
    hist_path, target_path = hist_paths(hist_path, target_hist_path)
    cpp = lambda v: ("true" if v else "false") if isinstance(v, bool) else (f'"{v}"' if isinstance(v, str) else repr(v))
    values = {
        "livetime_": float(exposure.livetime), "npot_": float(exposure.npot), "nmuons_": float(exposure.nmuons),
        "nmuons_per_pot_": float(exposure.nmuons_per_pot), "exposure_reference_": float(exposure.reference),
        "exposure_name_": exposure.name, "exposure_two_bb_": exposure.two_bb, "use_2bb_": exposure.two_bb,
        "exposure_target_": exposure.target or "", "muon_capture_fraction_": float(phys.muon_capture_fraction_al),
        "hist_path_str_": os.path.join(hist_path, ""), "target_hist_path_": os.path.join(target_path, ""),
    }
    ROOT.gROOT.ProcessLine(" ".join(f"{k} = {cpp(v)};" for k, v in values.items()) + " hist_path_ = hist_path_str_.Data();")
    if hasattr(ROOT, "datasets_"):
        ROOT.datasets_.clear()
        for name, info in dataset_info(exposure, reduced_tracker, hist_path, target_path).items():
            ROOT.datasets_[ROOT.TString(name)] = ROOT.DatasetInfo_t(info.ngen, info.nevents or 0, 0., 1., info.rate,
                                                                    info.dataset)


def rate_uncertainties(model, rmc_floating=False):
    """The normalization uncertainties of the mu- --> e-/e+ models, correlated across channels. The internal
    conversion fraction uncertainty (intN) covers both RPC and RMC, which use the same ratio. rmc_floating: the RMC
    rate is a free parameter fit to the control region (rmcNorm), so lumi and rmcN do not apply to the RMC."""
    model.lnN("lumi", 1.1, processes=["!cosmic*"] + (["!rmc*"] if rmc_floating else []), group="rate")
    model.lnN("csmN", 1.2, processes="cosmic*", group="rate")
    model.lnN("dioN", 1.025, processes="dio*", group="rate")
    model.lnN("rpcN", 1.27, processes="rpc*", group="rate")
    model.lnN("pbrN", 2.0, processes="pbar*", group="rate")
    if rmc_floating:
        model.rate_param("rmcNorm", "rmc*", range=(0., 10.))
    else:
        model.lnN("rmcN", 1.079, processes="rmc*", group="rate")   # 1.40 +- 0.11 (TRIUMF, 1999)
    model.lnN("intN", 1.045, processes="*_int*", group="rate")
