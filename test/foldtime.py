#!/usr/bin/env python3
"""How long does a construct take to fold?

Reports two different things, which are easy to confuse:

  * First-passage time to the target structure.  How many sweeps until each
    reference helix, and then all of them together, first appear.  This is a
    first-passage distribution over independent seeds, and some seeds never
    arrive within the budget, so it is summarized by median and quartiles over
    the seeds that did, alongside the fraction that did.  A mean would not be
    defined without assuming a distribution for the censored runs.

  * Equilibration of the energy.  How many sweeps until the seed-averaged fold
    energy reaches a plateau.  The energy settles long before the reference
    structure is found, so this is the weaker and more optimistic notion of
    "equilibrated"; it says the simulation has stopped drifting, not that it
    has found the native fold.

One sweep is one attempted move per mobile unit, so moves = sweeps x length.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import carnavallib as lib
from constructs import CONSTRUCTS


def equilibration(traces, frac=0.95, tail=0.2):
    """Sweep at which the seed-averaged energy first reaches `frac` of plateau.

    traces: list of per-seed lists of (sweeps, energy), all on the same grid.
    """
    n = min(len(t) for t in traces)
    if n < 5:
        return None, None, None
    grid = [traces[0][k][0] for k in range(n)]
    mean = [sum(t[k][1] for t in traces) / len(traces) for k in range(n)]
    start = max(1, int(n * (1 - tail)))
    plateau = sum(mean[start:]) / len(mean[start:])
    base = mean[0]
    if plateau <= base:
        return None, plateau, mean
    target = base + frac * (plateau - base)
    for k in range(n):
        if mean[k] >= target:
            return grid[k], plateau, mean
    return None, plateau, mean


def report(construct, seeds, unit_moves, samples, require_core):
    period = max(1, (unit_moves * len(construct["seq"])) // samples)
    nt = len(construct["seq"])
    runs = []
    traces = []
    for seed in seeds:
        s, _ = lib.run(construct, seed, unit_moves, period)
        runs.append(s)
        traces.append([(x.sweeps, x.energy) for x in s])

    budget = max(r[-1].sweeps for r in runs if r)
    print(f"\n=== {construct['name']}  ({nt} nt, {len(seeds)} seeds, "
          f"budget {budget:,.0f} sweeps = {budget*nt:,.0f} moves) ===")
    print(f"{'target':<22} {'reached':>8} {'median':>13} {'q25':>13} {'q75':>13} "
          f"{'min':>13} {'max':>13}")

    def line(label, pred):
        times = [lib.first_passage(r, pred) for r in runs]
        st = lib.censored_summary(times, len(runs), budget)
        print(f"{label:<22} {st['n_reached']}/{st['n']:<6} "
              f"{lib.fmt_sweeps(st['median']):>13} {lib.fmt_sweeps(st['q25']):>13} "
              f"{lib.fmt_sweeps(st['q75']):>13} {lib.fmt_sweeps(st['min']):>13} "
              f"{lib.fmt_sweeps(st['max']):>13}")
        return st

    for name, helix in construct["helices"].items():
        line(name, lambda s, h=helix: lib.helix_formed(s.pairs, h))
    allst = line("ALL helices", lambda s: lib.folded(s, construct))
    compst = None
    if construct["core"]:
        compst = line("ALL + core unpaired",
                      lambda s: lib.folded(s, construct, require_core=True))
    if construct["tertiary"]:
        line("tertiary pair", lambda s: all(p in s.pairs for p in construct["tertiary"]))

    eq, plateau, mean = equilibration(traces)
    e0 = mean[0] if mean else 0
    print(f"\n  energy: starts {e0:.1f}, plateau {plateau:.1f}" if plateau is not None
          else "\n  energy: no plateau detected")
    if eq is not None:
        print(f"  energy reaches 95% of plateau at {eq:,.0f} sweeps "
              f"= {eq*nt:,.0f} moves")
    occupancy = [sum(1 for x in r if lib.folded(x, construct)) / len(r) for r in runs]
    om, osem = lib.mean_sem(occupancy)
    print(f"  fraction of sampled time in the target structure: "
          f"{om:.3f} +- {osem:.3f}  (per seed: "
          f"{', '.join(f'{o:.2f}' for o in occupancy)})")
    return allst, compst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--construct", default="all", choices=["all"] + list(CONSTRUCTS))
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--unit-moves", type=int, default=1000000)
    ap.add_argument("--samples", type=int, default=400)
    ap.add_argument("--long", action="store_true",
                    help="10x the budget, matching the run length used in the paper")
    ap.add_argument("--require-core", action="store_true")
    args = ap.parse_args()
    if args.long:
        args.unit_moves *= 10
    names = list(CONSTRUCTS) if args.construct == "all" else [args.construct]
    seeds = list(range(1, args.seeds + 1))
    for name in names:
        report(CONSTRUCTS[name], seeds, args.unit_moves, args.samples, args.require_core)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
