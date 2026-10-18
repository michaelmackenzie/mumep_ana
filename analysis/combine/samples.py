"""Mu2eEvtAna (EventNtuple) inputs of the mu- --> e-/e+ search for mu2ecombine: datasets, rates, and exposures.

    from samples import Exposure, make_samples
    exposure = Exposure("run1a", two_bb=False)
    samples = make_samples(exposure)
    h = samples["dio"].hist("trk_20/obs")        # TH1D at the expected yield for the exposure

The rates follow datasets.C (EventNtuple inputs) and physics.C; the shared constants come from Mu2eEvtAna's
plotter/physics/Mu2ePhysics.C through mu2ecombine.mu2e_physics.
"""

import os

from mu2ecombine import Sample, mu2e_physics as phys

DEFAULT_HIST_PATH = "/exp/mu2e/data/projects/run1a/mumep_ana/histograms"
DEFAULT_REFERENCE = 1e-15   # branching fraction the signal is normalized to: r is in units of it


class Exposure:
    """N(POT), livetime, and stops per POT of a data set.

    name:   run1a (28 days of running), run2, mds3c, or mds3d (MDS ensembles, 1BB)
    two_bb: the Run 1A/2 exposure at the 2BB duty cycle and intensity, using the inputs re-weighted to 2BB
            (Mu2eEvtAna histogram mode 2, ConvAna.*.m2.root)
    """

    NAMES = ("run1a", "run2", "mds3c", "mds3d")

    def __init__(self, name="run1a", two_bb=False):
        if name not in self.NAMES:
            raise ValueError(f"unknown exposure '{name}' (available: {', '.join(self.NAMES)})")
        if two_bb and name.startswith("mds"):
            raise ValueError(f"the {name} ensemble is 1BB data, it has no 2BB exposure")
        self.name, self.two_bb = name, two_bb
        self.nmuons_per_pot = phys.run1a.nmuons_per_pot
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
        return (f"Exposure({self.name}{', 2BB' if self.two_bb else ''}: N(POT) = {self.npot:.3e}, "
                f"livetime = {self.livetime:.3e} s, N(muons) = {self.nmuons:.3e})")


def make_samples(exposure, hist_path=None):
    """{name: Sample} of the signal, background, and data inputs, normalized to the exposure. The signal is
    normalized to exposure.reference. hist_path: the histogram directory (default: CONVANA_HIST_PATH, or the
    analysis histogram directory)."""
    hist_path = hist_path or os.environ.get("CONVANA_HIST_PATH") or DEFAULT_HIST_PATH
    npot = exposure.npot
    nmu_per_pot = exposure.nmuons_per_pot
    f_cap = phys.muon_capture_fraction_al

    rate_signal = f_cap * nmu_per_pot * exposure.reference   # R_mue is relative to muon capture
    rate_dio_95 = (1. - f_cap) * nmu_per_pot * phys.dio_frac_above_95
    rate_rmc = nmu_per_pot * f_cap * phys.rmc_br_above_57
    rate_rmc_0n_80 = rate_rmc * phys.rmc_ps_0n_r_above_57 / phys.rmc_ps_0n_frac_above_57 * phys.rmc_ps_0n_frac_above_80
    rate_rmc_1n_80 = rate_rmc * phys.rmc_ps_1n_r_above_57 / phys.rmc_ps_1n_frac_above_57 * phys.rmc_ps_1n_frac_above_80
    rate_rpc = phys.run1a.pion_stops_per_pot * phys.rpc_br * phys.run1a.pion_survive_frac
    rmc_int, rpc_int = phys.rmc_internal_ratio, phys.rpc_internal_ratio
    cosmic_livetime = 4437713.0 * 5.72 / 4.438   # generated livetime x stored-livetime bug correction

    def sample(dataset, weight, expected_events=None, mode=exposure.hist_mode):
        path = os.path.join(hist_path, f"ConvAna.cnv_ana.{dataset}.m{mode}.root")
        return Sample(path, weight=weight, expected_events=expected_events, name=dataset)

    # name: dataset, per-event weight (exposure x rate / N(generated)), and events in the full dataset
    samples = {
        "mumem":      sample("cele1b1s5r0101", npot * rate_signal / 1e7, 4141125),
        "mumep":      sample("cpos1b1s5r0101", npot * rate_signal / 1e7, 3235878),
        "cosmic":     sample("cry4ab1s5r0101", exposure.livetime / cosmic_livetime, 4155435),
        "dio":        sample("dio00b1s5r0101", npot * rate_dio_95 / 25e6, 9368976),
        "rmc_ext_0n": sample("rmce0b1s5r0101", npot * rate_rmc_0n_80 / 7e9, 4967393),
        "rmc_ext_1n": sample("rmce1b1s5r0101", npot * rate_rmc_1n_80 / 7e9, 2974188),
        "rmc_int_0n": sample("rmci0b1s5r0101", npot * rate_rmc_0n_80 * rmc_int / 5e7, 1229006),
        "rmc_int_1n": sample("rmci1b1s5r0101", npot * rate_rmc_1n_80 * rmc_int / 5e7, 517029),
        "rpc_ext":    sample("rpce1b1s5r0101", npot * rate_rpc / 5e9, 458818),
        "rpc_int":    sample("rpci1b1s5r0101", npot * rate_rpc * rpc_int / 1.25e8, 1899806),
        "pbar":       sample("pbar1b1s5r0101", npot * phys.run1a.pbar_stops_per_pot / 3e7, 6461314),
    }
    # mu- --> e- keeps a single component per RMC source, taken from the 0 neutron knockout sample
    samples["rmc_ext"] = samples["rmc_ext_0n"]
    samples["rmc_int"] = samples["rmc_int_0n"]
    # data (1BB ensembles), unweighted
    samples["data_mds3c"] = sample("mds3cb1s5r0100", 1., mode=1)
    samples["data_mds3d"] = sample("mds3db0s5r0100", 1., mode=1)
    return samples


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
