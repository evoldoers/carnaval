"""Shared helpers: run the simulator, parse its fold log, score structures."""

import os
import subprocess
import statistics
import math

LEFT = "<[{(abcdefghijklmnopqrstuvwxyz"
RIGHT = ">]})ABCDEFGHIJKLMNOPQRSTUVWXYZ"

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.environ.get("CARNAVAL", os.path.join(HERE, os.pardir, "bin", "carnaval"))


def parse_pairs(fold):
    """Base pairs encoded in a CARNAVAL fold string, across all pseudoknot levels."""
    pairs = set()
    for level in range(len(LEFT)):
        stack = []
        for i, c in enumerate(fold):
            if c == LEFT[level]:
                stack.append(i)
            elif c == RIGHT[level] and stack:
                pairs.add((stack.pop(), i))
    return pairs


class Sample:
    __slots__ = ("pairs", "energy", "sweeps", "accepted")

    def __init__(self, pairs, energy, sweeps, accepted):
        self.pairs = pairs
        self.energy = energy
        self.sweeps = sweeps
        self.accepted = accepted


def run(construct, seed, unit_moves, period, extra=()):
    """Run one simulation and return (samples, stderr_text)."""
    x, y, z = construct["board"]
    cmd = [BINARY, "-x", str(x), "-y", str(y), "-z", str(z),
           "--init", construct["seq"],
           "--unit-moves", str(unit_moves),
           "--rnd", str(seed),
           "--folds", "--period", str(period), "--monochrome"] + list(extra)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"carnaval failed: {proc.stderr.strip()}")
    alphabet = set(LEFT + RIGHT + ".")
    samples = []
    for line in proc.stdout.splitlines():
        tok = line.split()
        if len(tok) < 5 or set(tok[2]) - alphabet:
            continue
        sweeps = next((float(t.split("=", 1)[1]) for t in tok if t.startswith("sweeps=")), None)
        if sweeps is None:
            continue
        samples.append(Sample(parse_pairs(tok[2]), float(tok[3]), sweeps, int(tok[0])))
    return samples, proc.stderr


def helix_formed(pairs, helix, tolerate_fraying=True):
    """A helix counts as formed when all its pairs are present, or all but one
    terminal pair if fraying is tolerated and the helix is longer than 2 bp."""
    present = sum(1 for p in helix if p in pairs)
    need = len(helix)
    if tolerate_fraying and len(helix) > 2:
        need = len(helix) - 1
    return present >= need


def core_unpaired(pairs, core):
    if not core:
        return True
    bound = {i for pair in pairs for i in pair}
    return not any(c in bound for c in core)


def folded(sample, construct, require_core=False):
    if not all(helix_formed(sample.pairs, h) for h in construct["helices"].values()):
        return False
    return core_unpaired(sample.pairs, construct["core"]) if require_core else True


def first_passage(samples, predicate):
    """Sweeps at which `predicate` first holds, or None if it never does."""
    for s in samples:
        if predicate(s):
            return s.sweeps
    return None


def quantile(sorted_vals, q):
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (pos - lo)


def censored_summary(times, n_total, budget):
    """Summarize first-passage times where some runs never reached the target.

    Returns a dict with the observed quantiles and the censoring fraction.  No
    mean is reported: with censored observations a mean is not defined without
    assuming a distribution, and the spread here is wide enough that the median
    is the honest summary.
    """
    obs = sorted(t for t in times if t is not None)
    n_cens = n_total - len(obs)
    return {
        "n": n_total,
        "n_reached": len(obs),
        "n_censored": n_cens,
        "frac_reached": len(obs) / n_total if n_total else 0.0,
        "min": obs[0] if obs else None,
        "q25": quantile(obs, 0.25),
        "median": quantile(obs, 0.50),
        "q75": quantile(obs, 0.75),
        "max": obs[-1] if obs else None,
        "budget": budget,
    }


def fmt_sweeps(v):
    return "never" if v is None else f"{v:,.0f}"


def mean_sem(vals):
    if not vals:
        return (None, None)
    if len(vals) == 1:
        return (vals[0], 0.0)
    return (statistics.mean(vals), statistics.stdev(vals) / math.sqrt(len(vals)))
