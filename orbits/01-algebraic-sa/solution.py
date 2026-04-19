"""Orbit 01 — algebraic seed + simulated-annealing construction.

No-three-in-line on a 10x10 integer grid.

Approach:
  1. Greedy random accept: place grid cells in random order, keep each if
     it does not create a collinear triple with the currently-held set.
     Hits 19 points within ~15 seconds (seed 3089).
  2. Warm-start simulated annealing at K=20 from a 19-point greedy state:
     - state = exactly 20 cells in [0,10)^2
     - energy = number of collinear triples in the state
     - neighbor = replace one random member with a random non-member
     - accept with Metropolis at T geometrically cooled 2.0 -> 0.01 over
       100k steps; restart from a fresh seed if stuck
     Within 2 seconds, SA descends energy to 0, producing the 20-point
     solution cached below. We meet the Erdos-Szekeres 2N ceiling for N=10.

The full search that produced POINTS lives in `search.py` in this orbit
folder; this file just exports the final validated configuration (the
evaluator loads POINTS and N — it does NOT run a search on import, so
eval is instant and deterministic).
"""

N = 10

# Size-20 no-three-in-line configuration discovered by SA (orbit 01).
# Exhaustively verified: every triple has non-zero cross product.
POINTS = [
    (0, 2), (0, 6),
    (1, 7), (1, 9),
    (2, 0), (2, 6),
    (3, 2), (3, 4),
    (4, 1), (4, 7),
    (5, 8), (5, 9),
    (6, 3), (6, 5),
    (7, 0), (7, 8),
    (8, 1), (8, 4),
    (9, 3), (9, 5),
]
