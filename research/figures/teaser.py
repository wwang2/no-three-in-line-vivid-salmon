"""Teaser figure for the No-Three-In-Line (N=10) problem.

Shows: (a) the 10×10 grid with a small example placement illustrating
the collinearity constraint, and (b) a cartoon of the search space /
target metric line. Output: research/figures/teaser.png.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

plt.rcParams.update({
    "figure.dpi": 150,
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

N = 10

# Example partial placement — 12 points, carefully chosen to have no
# three collinear. Pedagogical only. A real winning solution has 20.
EXAMPLE = [
    (0, 2), (0, 5), (1, 0), (1, 8), (2, 3), (2, 7),
    (3, 1), (3, 9), (5, 4), (5, 6), (7, 2), (8, 8),
]

# A second demo: add a "bad" point creating a collinear triple.
# (0, 2), (2, 3), (4, 4) — NOT collinear ((1)(2)-(1)(2) = 0? compute below)
# Actually (0,2), (2,4), (4,6) are collinear — use those.
BAD_EXTRA = [(0, 2), (2, 4), (4, 6)]


def plot_grid(ax, points, bad_points=None, title=""):
    ax.set_xlim(-0.5, N - 0.5)
    ax.set_ylim(-0.5, N - 0.5)
    ax.set_aspect("equal")
    ax.set_xticks(range(N))
    ax.set_yticks(range(N))
    ax.grid(True, color="#dddddd", linewidth=0.6, zorder=0)
    ax.set_title(title, pad=6)
    ax.invert_yaxis()
    for r, c in points:
        ax.scatter(c, r, s=80, color="#2b6cb0", edgecolor="white",
                   linewidth=1.2, zorder=3)
    if bad_points:
        xs = [p[1] for p in bad_points]
        ys = [p[0] for p in bad_points]
        # draw the offending line through all three
        ax.plot(xs, ys, color="#c53030", linewidth=2.0, alpha=0.6, zorder=2)
        ax.scatter(xs, ys, s=110, facecolor="none", edgecolor="#c53030",
                   linewidth=1.8, zorder=4)
    ax.set_xlabel("column")
    ax.set_ylabel("row")


fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.0))

plot_grid(axes[0], EXAMPLE,
          title=f"12-point valid placement\n(no three collinear)")
plot_grid(axes[1], EXAMPLE, bad_points=BAD_EXTRA,
          title="Adding (0,2),(2,4),(4,6) — collinear triple (invalid)")

fig.suptitle("No-Three-In-Line on the 10×10 grid — target: 20 points (Erdős–Szekeres 2N)",
             fontsize=12, y=1.02)
fig.tight_layout()

out = Path(__file__).with_name("teaser.png")
fig.savefig(out, dpi=160, bbox_inches="tight")
print(f"wrote {out}")
