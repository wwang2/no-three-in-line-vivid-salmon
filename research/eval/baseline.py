"""Sanity baseline: diagonal-only placement.

10 points along a single diagonal is collinear (three on one line),
so this returns metric=0 (invalid). Used ONLY to prove the evaluator
runs end-to-end. Real orbit solutions live in orbits/<name>/solution.py.
"""

N = 10
POINTS = [(i, i) for i in range(10)]
