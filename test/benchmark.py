#!/usr/bin/env python3
"""Codified throughput benchmark: fold the 1ehz tRNA and measure move rates.

Records the measurement against test/baseline.json and fails if throughput has
regressed by more than a tolerance factor, so this doubles as a performance
regression test.  Pass --update to record the current machine as the baseline.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import carnavallib as lib
from constructs import TRNA

BASELINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "baseline.json")


def one_run(unit_moves):
    x, y, z = TRNA["board"]
    cmd = [lib.BINARY, "-x", str(x), "-y", str(y), "-z", str(z),
           "--init", TRNA["seq"], "--unit-moves", str(unit_moves),
           "--rnd", "1", "--folds", "--period", str(10 ** 12), "--monochrome"]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip())
    tried = int(re.search(r"Tried (\d+) moves", proc.stderr).group(1))
    accepted = int(re.search(r"moves, (\d+) succeeded", proc.stderr).group(1))
    sweeps = float(re.search(r"Simulated ([0-9.e+]+) sweeps", proc.stderr).group(1))
    return elapsed, tried, accepted, sweeps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unit-moves", type=int, default=1000000,
                    help="moves per mobile unit; 1e6 is 76M moves for this 76-mer")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--tolerance", type=float, default=2.0,
                    help="fail if accepted moves/s is worse than baseline by this factor")
    ap.add_argument("--update", action="store_true", help="record this run as the baseline")
    args = ap.parse_args()

    runs = [one_run(args.unit_moves) for _ in range(args.repeats)]
    best = min(r[0] for r in runs)
    elapsed, tried, accepted, sweeps = min(runs, key=lambda r: r[0])

    print(f"tRNA folding benchmark ({TRNA['name']}, {len(TRNA['seq'])} nt, "
          f"board {TRNA['board'][0]}x{TRNA['board'][1]})")
    print(f"  attempted moves      {tried:,}")
    print(f"  accepted moves       {accepted:,}  ({100.0*accepted/tried:.2f}% of attempts)")
    print(f"  sweeps               {sweeps:,.0f}")
    print(f"  wall time, best of {args.repeats}  {best:.2f}s   "
          f"(all: {', '.join(f'{r[0]:.2f}' for r in runs)})")
    print(f"  attempted moves/s    {tried/best/1e6:.2f}M")
    print(f"  accepted moves/s     {accepted/best/1e6:.2f}M")
    print(f"  sweeps/s             {sweeps/best:,.0f}")

    record = {"unit_moves": args.unit_moves, "attempted": tried, "accepted": accepted,
              "seconds": round(best, 4), "accepted_per_s": accepted / best}
    if args.update:
        with open(BASELINE, "w") as fh:
            json.dump(record, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"\nbaseline written to {BASELINE}")
        return 0

    if not os.path.exists(BASELINE):
        print(f"\nno baseline at {BASELINE}; run with --update to record one")
        return 0
    with open(BASELINE) as fh:
        base = json.load(fh)
    if base["unit_moves"] != args.unit_moves:
        print(f"\nbaseline used --unit-moves {base['unit_moves']}; not comparing")
        return 0
    ratio = record["accepted_per_s"] / base["accepted_per_s"]
    print(f"\nbaseline {base['accepted_per_s']/1e6:.2f}M accepted/s, "
          f"this run {record['accepted_per_s']/1e6:.2f}M, ratio {ratio:.2f}x")
    if ratio < 1.0 / args.tolerance:
        print(f"REGRESSION: slower than baseline by more than {args.tolerance}x")
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
