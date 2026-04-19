"""Offline search driver — runs the solver across many seeds in parallel
and prints the best POINTS configuration found.

Run: python3 orbits/01-algebraic-sa/_search.py --seeds 16 --budget 90
Then hardcode the printed POINTS into solution.py.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
SOLUTION_PY = HERE / "solution.py"


def worker(seed: int, budget: float):
    """Run solve() with given seed and budget via subprocess."""
    env = os.environ.copy()
    env["NO3_SEED"] = str(seed)
    env["NO3_BUDGET"] = str(budget)
    t0 = time.time()
    result = subprocess.run(
        [sys.executable, "-c", (
            "import os, sys\n"
            f"sys.path.insert(0, {str(HERE.parent.parent)!r})\n"
            f"import importlib.util\n"
            f"spec = importlib.util.spec_from_file_location('s', {str(SOLUTION_PY)!r})\n"
            "m = importlib.util.module_from_spec(spec)\n"
            "spec.loader.exec_module(m)\n"
            "pts = m.POINTS\n"
            "ok = m._verify_no_three(pts)\n"
            "print('LEN', len(pts))\n"
            "print('OK', ok)\n"
            "print('PTS', pts)\n"
        )],
        env=env, capture_output=True, text=True, timeout=budget + 60,
    )
    dt = time.time() - t0
    if result.returncode != 0:
        print(f"  seed={seed} FAILED:\n{result.stderr}", file=sys.stderr)
        return seed, [], False, dt
    out = result.stdout
    pts = []
    ok = False
    for line in out.splitlines():
        if line.startswith("PTS "):
            pts = eval(line[4:])
        elif line.startswith("OK "):
            ok = (line[3:].strip() == "True")
    return seed, pts, ok, dt


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, default=8)
    p.add_argument("--budget", type=float, default=60.0)
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()

    seeds = list(range(1, args.seeds + 1))
    print(f"Launching {args.seeds} trials, budget={args.budget}s each, workers={args.workers}...",
          file=sys.stderr)

    best_size = 0
    best_pts = None
    best_seed = None

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(worker, s, args.budget): s for s in seeds}
        for fut in concurrent.futures.as_completed(futures):
            seed, pts, ok, dt = fut.result()
            print(f"  seed={seed}: {len(pts)} pts, valid={ok}, {dt:.1f}s",
                  file=sys.stderr)
            if ok and len(pts) > best_size:
                best_size = len(pts)
                best_pts = pts
                best_seed = seed

    print(f"\nBEST: seed={best_seed}, size={best_size}", file=sys.stderr)
    print(f"BEST_POINTS = {best_pts}")


if __name__ == "__main__":
    main()
