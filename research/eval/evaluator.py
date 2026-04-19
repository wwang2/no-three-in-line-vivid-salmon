#!/usr/bin/env python3
"""
Evaluator for the no-three-in-line benchmark (Erdős extremal, 1917 / Dudeney 1917).

Task: place as many points as possible on an N×N integer grid such that
no three of them are collinear. The solution file must export:

    N = 10          # grid size (must equal 10 here — the agent cannot change it)
    POINTS = [...]  # list of (row, col) integer tuples, 0 <= row, col < N

Metric (maximize):
    n_valid = len(POINTS) if no three are collinear else 0

The evaluator:
  1. Validates points are distinct integers in [0, N).
  2. Checks every triple (i, j, k) for collinearity via integer cross product
     (p_j - p_i) × (p_k - p_i) == 0.
  3. Returns -len(POINTS) as the metric (negative because direction:minimize
     in the existing /autorun harness). Agents should try to make this as
     negative as possible — i.e., pack as many points in as they can.

Background:
  - Erdős-Szekeres conjecture (1951): on any N × N grid, max is 2N points.
  - Proven for N ≤ 46 (numerically). The conjecture itself is still open
    for general N.
  - For N=10: best known = 20 (meets 2N).
  - For N=12: best known = 22 (below 2N=24; the conjecture may fail at
    larger N).

  Multi-orbit motivation: agents exploring naive lattices, diagonals, or
  affine constructions will easily hit ~12-15 points before running into
  collinearity floors. Reaching 18-20 requires algebraic constructions
  (e.g., (r, r² mod p) curves modulo a prime p close to N) or careful
  local search. A single-shot agent rarely hits the 2N ceiling — the
  problem is designed to exercise 3-4 orbits of iteration.

Usage:
    python evaluator.py --solution <path> --seed <int>
    Output: METRIC=<float>  (NEGATIVE point count; direction:minimize)
"""

import argparse
import importlib.util
import os
import sys


def load_solution(path):
    spec = importlib.util.spec_from_file_location("solution", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cross(o, a, b):
    """Integer cross product (a - o) × (b - o). Zero iff collinear."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def validate_points(points, N):
    """Check points are distinct, integer, and within the N × N grid."""
    if not isinstance(points, (list, tuple)):
        return "POINTS must be a list or tuple"
    seen = set()
    for i, p in enumerate(points):
        if not (isinstance(p, (list, tuple)) and len(p) == 2):
            return f"points[{i}] is not a length-2 (row, col) tuple"
        r, c = p
        if not (isinstance(r, int) and isinstance(c, int)):
            return f"points[{i}] has non-integer coords: {p}"
        if not (0 <= r < N and 0 <= c < N):
            return f"points[{i}] = {p} is outside the {N}×{N} grid"
        if (r, c) in seen:
            return f"points[{i}] = {p} is a duplicate"
        seen.add((r, c))
    return None


def check_no_collinear(points):
    """Returns (i, j, k) indices of first collinear triple, or None if none."""
    n = len(points)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if cross(points[i], points[j], points[k]) == 0:
                    return (i, j, k)
    return None


def evaluate(solution_path, seed=1):
    try:
        module = load_solution(solution_path)
    except Exception as e:
        print(f"ERROR: Could not load solution: {e}", file=sys.stderr)
        sys.exit(1)

    N = getattr(module, "N", None)
    if N != 10:
        print(
            f"ERROR: this benchmark requires N=10 (got {N}); do not change it.",
            file=sys.stderr,
        )
        sys.exit(1)

    points = list(getattr(module, "POINTS", []))
    err = validate_points(points, N)
    if err:
        print(f"ERROR: invalid POINTS — {err}", file=sys.stderr)
        # Invalid submission scored as 0 — agent gets penalized (metric=0 on
        # a minimize axis = worst possible; real answers give metric < 0).
        print("METRIC=0.0")
        return

    triple = check_no_collinear(points)
    if triple is not None:
        i, j, k = triple
        print(
            f"INVALID: points[{i}]={points[i]}, points[{j}]={points[j]}, "
            f"points[{k}]={points[k]} are collinear",
            file=sys.stderr,
        )
        print("METRIC=0.0")
        return

    # Score = -count, so minimize = pack more points.
    metric = -float(len(points))
    print(f"METRIC={metric:.6f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--solution", required=True)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    evaluate(args.solution, args.seed)
