---
issue: 2
parents: []
eval_version: eval-v1
metric: null
---

# Research Notes

## Strategy

Hybrid algebraic + simulated-annealing approach for the no-three-in-line
problem on the 10x10 grid (target: 20 points, i.e. metric = -20).

### Families to try
- Algebraic seeds: quadratic-residue curves y = x^2 mod p for primes near N
  (p = 11 is the natural choice since 10 < 11); modular hyperbolas
  y = x^{-1} mod p; random-translated copies of these.
- Simulated annealing with objective = point count − λ · (# collinear triples),
  restarted from multiple algebraic/random seeds.
- Local repair: given a high-count candidate, remove/add single points to
  eliminate remaining collinear triples.

