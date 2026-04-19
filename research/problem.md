# No-Three-In-Line on a 10×10 Grid (Erdős Extremal)

## Problem Statement

Classical Erdős extremal problem (Dudeney 1917, Erdős 1917; Erdős–Szekeres
conjecture 1951): place as many points as possible on a 10×10 integer grid
such that **no three of them are collinear**. The best-known / known
optimum for N=10 is **20 points**, which meets the Erdős–Szekeres 2N
ceiling.

## Solution Interface

Submission file: `orbits/<name>/solution.py` exporting exactly

```python
N = 10
POINTS = [(r, c), ...]   # integer tuples, 0 <= r, c < 10
```

The evaluator at `research/eval/evaluator.py` validates distinctness and
in-range integer coordinates, then checks all triples for collinearity via
integer cross product. No library restrictions on the solution file
itself, but the problem solvers must NOT use `sklearn` or
`scipy.optimize` — the search is discrete (local search / simulated
annealing / algebraic constructions / hybrids).

## Success Metric

`METRIC = -len(POINTS)` if the point set is valid (no three collinear),
else `0`.

- Direction: **minimize**.
- Target: **-20** (the known optimum for N=10 — anything less does not
  meet target).
- 0 = invalid submission (collinear triple or malformed).

## Budget

- Max 4 orbits. The problem is designed to require multiple rounds of
  hypothesis revision: a single-shot agent can usually hit ~14–16 points
  via naive lattices/diagonals, but reaching 18–20 typically requires
  algebraic constructions (e.g. quadratic-residue curves modulo a prime
  near N, modular hyperbolas) combined with local repair.
- Each orbit AFTER orbit 01 must **EXTEND** from the best prior orbit's
  construction — reference the parent orbit in `log.md` and the Issue
  body.
