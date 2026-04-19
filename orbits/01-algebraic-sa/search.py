"""Hybrid algebraic-seed + simulated-annealing search for no-three-in-line.

Strategy (orbit 01):
  1. Generate algebraic seeds on an N x N grid:
       - quadratic residue curve  { (x, x^2 mod p) }  for prime p near N
       - modular hyperbola        { (x, x^{-1} mod p) }
       - translated / reflected copies of the above
  2. Filter seeds into the 10x10 window and greedily prune collinear triples.
  3. Run simulated annealing on a target size K (we aim K = 18, 19, 20).
     - State: a sorted tuple of K distinct cells in [0,N)^2.
     - Energy: (# collinear triples)  [lower is better; 0 means valid]
     - Neighbor: replace one member with a random non-member cell.
     - Cool slowly with occasional random restarts.
  4. Multi-restart with different seeds / starting configs; keep the best
     VALID configuration we can find.

The approach is all pure numpy + stdlib (no scipy.optimize, no sklearn).
"""

from __future__ import annotations

import itertools
import math
import random
from collections import defaultdict
from typing import Iterable, List, Sequence, Set, Tuple

import numpy as np

N = 10
Point = Tuple[int, int]


# --------------------------------------------------------------------------- #
# Algebraic seed constructions
# --------------------------------------------------------------------------- #
def quadratic_residue_seed(p: int, N: int = N) -> List[Point]:
    """Points (x, x^2 mod p) with x in {0..p-1}, clipped to the NxN window."""
    pts = [(x % N, (x * x) % p) for x in range(p)]
    return [(r, c) for (r, c) in pts if 0 <= r < N and 0 <= c < N]


def modular_hyperbola_seed(p: int, N: int = N) -> List[Point]:
    """Points (x, x^{-1} mod p) for x in {1..p-1}, clipped to NxN."""
    out = []
    for x in range(1, p):
        try:
            inv = pow(x, -1, p)
        except ValueError:
            continue
        if 0 <= x < N and 0 <= inv < N:
            out.append((x, inv))
    return out


def cubic_seed(p: int, N: int = N) -> List[Point]:
    pts = [(x % N, (x ** 3) % p) for x in range(p)]
    return [(r, c) for (r, c) in pts if 0 <= r < N and 0 <= c < N]


def algebraic_union(N: int = N) -> List[Point]:
    """Union of a handful of algebraic seeds (unique points)."""
    seeds: Set[Point] = set()
    for p in (11, 13):
        seeds.update(quadratic_residue_seed(p, N))
        seeds.update(modular_hyperbola_seed(p, N))
        seeds.update(cubic_seed(p, N))
    return sorted(seeds)


# --------------------------------------------------------------------------- #
# Collinearity machinery
# --------------------------------------------------------------------------- #
def triple_energy(points: Sequence[Point]) -> int:
    """Number of collinear triples among `points`. Uses line-slope hashing.

    For each point p, bucket every other point by reduced (dy, dx). Any
    bucket of size k contributes C(k,2) triples (i.e. (k choose 2) pairs of
    collinear partners -> each pair + p is a collinear triple, double-counted
    across the three rotations so we divide by 3 at the end).
    """
    n = len(points)
    if n < 3:
        return 0
    total = 0
    for i in range(n):
        buckets: defaultdict[Tuple[int, int], int] = defaultdict(int)
        r0, c0 = points[i]
        for j in range(n):
            if j == i:
                continue
            dr = points[j][0] - r0
            dc = points[j][1] - c0
            g = math.gcd(abs(dr), abs(dc))
            if g == 0:
                continue
            dr //= g
            dc //= g
            if dr < 0 or (dr == 0 and dc < 0):
                dr, dc = -dr, -dc
            buckets[(dr, dc)] += 1
        for k in buckets.values():
            if k >= 2:
                total += k * (k - 1) // 2
    # every collinear triple counted once per vertex (3 times)
    return total // 3


def lines_through(points: Sequence[Point]) -> List[List[int]]:
    """Return list of maximal collinear subsets (as index lists of size >= 3)."""
    n = len(points)
    line_map: defaultdict[Tuple[int, int, int], List[int]] = defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            r0, c0 = points[i]
            r1, c1 = points[j]
            dr = r1 - r0
            dc = c1 - c0
            g = math.gcd(abs(dr), abs(dc))
            dr //= g
            dc //= g
            if dr < 0 or (dr == 0 and dc < 0):
                dr, dc = -dr, -dc
            # line equation: dr*(y - y0) = dc*(x - x0)  ->  dr*y - dc*x = dr*y0 - dc*x0
            c_line = dr * c0 - dc * r0
            key = (dr, dc, c_line)
            line_map[key].append(i)
            line_map[key].append(j)
    out = []
    for _, idxs in line_map.items():
        u = sorted(set(idxs))
        if len(u) >= 3:
            out.append(u)
    return out


def is_valid(points: Sequence[Point]) -> bool:
    return triple_energy(points) == 0


def max_points_on_line(points: Sequence[Point]) -> int:
    """Max number of given points that lie on a single line (>=2)."""
    if len(points) < 2:
        return 1
    best = 2
    for i in range(len(points)):
        buckets: defaultdict[Tuple[int, int], int] = defaultdict(int)
        r0, c0 = points[i]
        for j in range(len(points)):
            if j == i:
                continue
            dr = points[j][0] - r0
            dc = points[j][1] - c0
            g = math.gcd(abs(dr), abs(dc))
            if g == 0:
                continue
            dr //= g
            dc //= g
            if dr < 0 or (dr == 0 and dc < 0):
                dr, dc = -dr, -dc
            buckets[(dr, dc)] += 1
        if buckets:
            best = max(best, 1 + max(buckets.values()))
    return best


# --------------------------------------------------------------------------- #
# Greedy valid subset (repair): keep points while adding none creates triples
# --------------------------------------------------------------------------- #
def greedy_valid_subset(candidates: Sequence[Point], seed: int = 0) -> List[Point]:
    """Greedy: shuffle candidates, add one by one if adding keeps validity."""
    rng = random.Random(seed)
    order = list(candidates)
    rng.shuffle(order)
    chosen: List[Point] = []
    for p in order:
        chosen.append(p)
        if not is_valid(chosen):
            chosen.pop()
    return chosen


def greedy_valid_grid(seed: int = 0, N: int = N) -> List[Point]:
    """Greedy over ALL grid cells in random order."""
    all_cells = [(r, c) for r in range(N) for c in range(N)]
    return greedy_valid_subset(all_cells, seed=seed)


# --------------------------------------------------------------------------- #
# Simulated annealing (targeted on size K)
# --------------------------------------------------------------------------- #
def delta_energy_remove(points: List[Point], idx: int) -> int:
    """Change in triple count if point at idx is removed."""
    p = points[idx]
    r0, c0 = p
    rest = [points[j] for j in range(len(points)) if j != idx]
    # triples through p = number of (a,b) pairs in rest collinear with p
    buckets: defaultdict[Tuple[int, int], int] = defaultdict(int)
    for (r, c) in rest:
        dr = r - r0
        dc = c - c0
        g = math.gcd(abs(dr), abs(dc))
        dr //= g
        dc //= g
        if dr < 0 or (dr == 0 and dc < 0):
            dr, dc = -dr, -dc
        buckets[(dr, dc)] += 1
    removed = 0
    for k in buckets.values():
        if k >= 2:
            removed += k * (k - 1) // 2
    return -removed


def delta_energy_add(points: List[Point], new_pt: Point) -> int:
    """Change in triple count if new_pt is added to `points`."""
    r0, c0 = new_pt
    buckets: defaultdict[Tuple[int, int], int] = defaultdict(int)
    for (r, c) in points:
        dr = r - r0
        dc = c - c0
        g = math.gcd(abs(dr), abs(dc))
        if g == 0:
            continue
        dr //= g
        dc //= g
        if dr < 0 or (dr == 0 and dc < 0):
            dr, dc = -dr, -dc
        buckets[(dr, dc)] += 1
    added = 0
    for k in buckets.values():
        if k >= 2:
            added += k * (k - 1) // 2
    return added


def simulated_anneal(
    K: int,
    *,
    seed: int,
    init: List[Point] | None = None,
    iters: int = 80_000,
    T0: float = 3.0,
    T_end: float = 0.02,
    N: int = N,
) -> Tuple[List[Point], int]:
    """Fixed-size SA: maintain exactly K points, minimize collinear triples."""
    rng = random.Random(seed)
    all_cells = [(r, c) for r in range(N) for c in range(N)]

    if init is not None and len(init) >= K:
        state = list(init[:K])
    elif init is not None and len(init) < K:
        state = list(init)
        remaining = [c for c in all_cells if c not in set(state)]
        rng.shuffle(remaining)
        state.extend(remaining[: K - len(state)])
    else:
        state = rng.sample(all_cells, K)

    state_set = set(state)
    E = triple_energy(state)
    best_state = list(state)
    best_E = E

    log_T0 = math.log(T0)
    log_Te = math.log(T_end)
    for step in range(iters):
        t = step / max(1, iters - 1)
        T = math.exp(log_T0 * (1 - t) + log_Te * t)

        # propose: swap one state member with a random non-member
        idx = rng.randrange(K)
        old_pt = state[idx]
        # sample a non-member
        while True:
            new_pt = all_cells[rng.randrange(len(all_cells))]
            if new_pt not in state_set:
                break

        dE_rem = delta_energy_remove(state, idx)  # negative or zero
        # after removal, compute delta for adding new_pt
        tmp_state = state[:idx] + state[idx + 1 :]
        dE_add = delta_energy_add(tmp_state, new_pt)
        dE = dE_rem + dE_add

        if dE <= 0 or rng.random() < math.exp(-dE / max(T, 1e-9)):
            state[idx] = new_pt
            state_set.remove(old_pt)
            state_set.add(new_pt)
            E += dE
            if E < best_E:
                best_E = E
                best_state = list(state)
                if best_E == 0:
                    # perfect: try to escape slightly? just keep going / break
                    break

    return best_state, best_E


# --------------------------------------------------------------------------- #
# Multi-restart orchestration
# --------------------------------------------------------------------------- #
def search_size(K: int, restarts: int = 12, iters: int = 60_000) -> List[Point] | None:
    """Try to find a VALID size-K configuration."""
    # build algebraic seed pool
    seeds: List[List[Point]] = []
    pool = algebraic_union(N)
    if pool:
        seeds.append(pool)
        # translated variants
        for shift in [(1, 0), (0, 1), (1, 1), (2, 3), (5, 5)]:
            sr, sc = shift
            seeds.append(
                sorted({((r + sr) % N, (c + sc) % N) for (r, c) in pool})
            )

    for t in range(restarts):
        seed = 1000 + t
        # choose an initial state: either algebraic-seed-based or empty
        init = None
        if t < len(seeds):
            init = list(seeds[t])
        final, E = simulated_anneal(
            K=K, seed=seed, init=init, iters=iters, T0=2.5, T_end=0.02
        )
        if E == 0 and is_valid(final):
            return sorted(final)
    return None


def best_effort(max_K: int = 20, min_K: int = 14, restarts: int = 12) -> List[Point]:
    """Try K from max_K downwards; return the first valid configuration found."""
    best: List[Point] = []
    for K in range(max_K, min_K - 1, -1):
        found = search_size(K, restarts=restarts, iters=50_000 + 3_000 * K)
        if found is not None:
            return found
    # fallback to greedy
    for s in range(50):
        cand = greedy_valid_grid(seed=s)
        if len(cand) > len(best):
            best = cand
    return best


if __name__ == "__main__":
    random.seed(1)
    np.random.seed(1)
    sol = best_effort()
    assert is_valid(sol), "invalid"
    print("len =", len(sol))
    print(sol)
