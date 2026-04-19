"""Generate narrative.png (qualitative) and results.png (quantitative)
for the no-three-in-line 10x10 orbit.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

HERE = Path(__file__).parent
FIGS = HERE / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

# Style from research/style.md
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.titleweight": "medium",
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "axes.grid": True,
    "grid.alpha": 0.15,
    "grid.linewidth": 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlepad": 10.0,
    "axes.labelpad": 6.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "legend.borderpad": 0.3,
    "legend.handletextpad": 0.5,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "figure.constrained_layout.use": True,
})

COLOR_POINT = "#4C72B0"
COLOR_VIOLATION = "#C44E52"
COLOR_BASELINE = "#888888"
COLOR_METHOD = "#55A868"

N = 10

SOLUTION = [
    (0, 4), (0, 6),
    (1, 1), (1, 2),
    (2, 1), (2, 7),
    (3, 6), (3, 9),
    (4, 0), (4, 4),
    (5, 0), (5, 8),
    (6, 2), (6, 5),
    (7, 8), (7, 9),
    (8, 3), (8, 5),
    (9, 3), (9, 7),
]


def collinear_triples(points):
    n = len(points)
    trips = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                (x1, y1), (x2, y2), (x3, y3) = points[i], points[j], points[k]
                if (x2 - x1) * (y3 - y1) == (x3 - x1) * (y2 - y1):
                    trips.append((i, j, k))
    return trips


# ---------- baseline: diagonal (invalid) ----------
BASELINE = [(i, i) for i in range(N)]  # invalid, 10 collinear


# ---------- narrative.png: baseline vs our method ----------
def draw_grid(ax, points, title, color_pts, color_violations=None, show_lines=True):
    """Draw the 10x10 grid with points + any collinear triples highlighted."""
    ax.set_xlim(-0.5, N - 0.5)
    ax.set_ylim(-0.5, N - 0.5)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_xticks(range(N))
    ax.set_yticks(range(N))
    ax.set_xticklabels(range(N))
    ax.set_yticklabels(range(N))
    ax.invert_yaxis()

    # shade violating collinear lines (for the baseline)
    viol_points = set()
    if show_lines and len(points) >= 3:
        trips = collinear_triples(points)
        if trips:
            # Draw lines through all collinear triples
            drawn_lines = set()
            for i, j, k in trips:
                p1, p2 = points[i], points[k]
                # canonical line id via direction
                dx, dy = p2[0] - p1[0], p2[1] - p1[1]
                g = math.gcd(abs(dx), abs(dy))
                if g > 0:
                    dx, dy = dx // g, dy // g
                # cross: all pairs with same line
                key = (dx, dy, p1[0] * dy - p1[1] * dx)  # cross-prod offset
                if key in drawn_lines:
                    pass
                drawn_lines.add(key)
                # draw segment from left extreme to right extreme on this line
                ax.plot(
                    [p1[1], p2[1]], [p1[0], p2[0]],
                    color=COLOR_VIOLATION, alpha=0.25, linewidth=1.4, zorder=1,
                )
                viol_points.update([points[i], points[j], points[k]])

    # points: plot non-violating in primary color, violating in red
    for (r, c) in points:
        if (r, c) in viol_points:
            ax.plot(c, r, "o", markersize=13, markerfacecolor=COLOR_VIOLATION,
                    markeredgecolor="white", markeredgewidth=1.2, zorder=3)
        else:
            ax.plot(c, r, "o", markersize=13, markerfacecolor=color_pts,
                    markeredgecolor="white", markeredgewidth=1.2, zorder=3)

    ax.set_xlabel("column")
    ax.set_ylabel("row")
    ax.set_title(title)


def make_narrative():
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.6), sharey=True)

    # LEFT: naive diagonal baseline (10 points — all collinear)
    draw_grid(axes[0], BASELINE,
              "Baseline: naive diagonal\n10 points — ALL collinear (invalid)",
              COLOR_BASELINE, show_lines=True)
    axes[0].text(
        -0.12, 1.05, "(a)", transform=axes[0].transAxes,
        fontsize=14, fontweight="bold",
    )
    axes[0].annotate(
        "metric = 0\n(invalid)",
        xy=(N - 1, N - 1), xycoords="data",
        xytext=(N + 1.1, N - 1.5), fontsize=10, color=COLOR_VIOLATION,
        ha="left",
    )

    # RIGHT: our 20-point solution
    draw_grid(axes[1], SOLUTION,
              "SA hybrid: 20-point solution\nno 3 collinear — meets target (2N ceiling)",
              COLOR_METHOD, show_lines=False)
    axes[1].text(
        -0.12, 1.05, "(b)", transform=axes[1].transAxes,
        fontsize=14, fontweight="bold",
    )
    axes[1].annotate(
        "metric = -20\n(optimal)",
        xy=(0, 0), xycoords="data",
        xytext=(N + 0.6, 0.3), fontsize=10, color=COLOR_METHOD,
        ha="left",
    )

    fig.suptitle(
        "No-three-in-line on 10×10: baseline failure vs. SA solution",
        fontsize=14, fontweight="medium", y=1.02,
    )
    out = FIGS / "narrative.png"
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


# ---------- results.png: quantitative panel ----------
# Simulated-annealing trial data (from _search.py run above)
SEARCH_TRIALS = [
    # (seed, n_points, time_sec)
    (8, 20, 8.5),
    (2, 20, 13.6),
    (4, 19, 105.0),
    (1, 19, 105.0),
    (6, 19, 105.1),
    (5, 19, 105.3),
    (7, 19, 105.3),
    (3, 19, 105.3),
]

# Evaluator results — 3 seeds, all on the same POINTS list → deterministic.
EVAL_SEEDS = [1, 2, 3]
EVAL_METRICS = [-20.0, -20.0, -20.0]

# Reference: baselines
BASELINE_METRIC = 0.0  # diagonal baseline = invalid
TARGET_METRIC = -20.0


def make_results():
    fig = plt.figure(figsize=(13, 5.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.2, 1.1], wspace=0.28)

    # Panel (a): histogram of trial sizes
    ax = fig.add_subplot(gs[0])
    sizes = [t[1] for t in SEARCH_TRIALS]
    bin_edges = np.arange(16.5, 22.5)
    counts, _ = np.histogram(sizes, bins=bin_edges)
    xs = np.arange(17, 22)
    colors = [COLOR_POINT if x < 20 else COLOR_METHOD for x in xs]
    ax.bar(xs, counts, color=colors, edgecolor="white", width=0.85)
    for x, c in zip(xs, counts):
        if c > 0:
            ax.text(x, c + 0.12, str(int(c)), ha="center",
                    fontsize=12, fontweight="medium",
                    color=COLOR_METHOD if x == 20 else COLOR_POINT)
    ax.text(20, max(counts) + 0.9, "TARGET",
            ha="center", fontsize=9, color=COLOR_METHOD, fontweight="medium")
    ax.set_xlabel("solution size (# points)")
    ax.set_ylabel("# of seeds")
    ax.set_title("SA trial outcomes (8 seeds)")
    ax.set_xticks(xs)
    ax.set_ylim(0, max(counts) + 1.7)
    ax.grid(True, axis="y", alpha=0.2)
    ax.text(-0.14, 1.05, "(a)", transform=ax.transAxes,
            fontsize=14, fontweight="bold")

    # Panel (b): time-to-solution
    ax = fig.add_subplot(gs[1])
    ts = np.array([t[2] for t in SEARCH_TRIALS])
    sz = np.array([t[1] for t in SEARCH_TRIALS])
    seeds_arr = [t[0] for t in SEARCH_TRIALS]
    colors = [COLOR_METHOD if s == 20 else COLOR_POINT for s in sz]
    ax.scatter(ts, sz, c=colors, s=130, edgecolor="white", linewidth=1.4, zorder=3,
               alpha=0.9)
    # Collapse clustered labels: 6 seeds hit 19 at ~105s; 2 hit 20 at ~10s.
    ax.text(11, 20.12,
            "seeds 8, 2 — 20 pts in < 14 s",
            fontsize=10, color=COLOR_METHOD, ha="left", fontweight="medium")
    ax.text(105, 18.82,
            "seeds 1, 3, 4, 5, 6, 7 — 19 pts at 105 s budget",
            fontsize=10, color=COLOR_POINT, ha="right")
    ax.axhline(20, color=COLOR_METHOD, linestyle="--", linewidth=1.0, alpha=0.5)
    ax.axhline(19, color=COLOR_POINT, linestyle=":", linewidth=0.8, alpha=0.3)
    ax.set_xlabel("wall-clock time (s)")
    ax.set_ylabel("solution size (# points)")
    ax.set_title("Time to reach each solution")
    ax.set_ylim(18.4, 20.6)
    ax.set_xlim(-5, 120)
    ax.text(-0.14, 1.05, "(b)", transform=ax.transAxes,
            fontsize=14, fontweight="bold")

    # Panel (c): final evaluator metric
    ax = fig.add_subplot(gs[2])
    xs = np.arange(len(EVAL_SEEDS))
    ax.bar(xs, EVAL_METRICS, color=COLOR_METHOD, edgecolor="white", width=0.6)
    ax.axhline(BASELINE_METRIC, color=COLOR_BASELINE,
               linestyle="--", linewidth=1.2, label="baseline (invalid, 0.0)")
    ax.axhline(TARGET_METRIC, color=COLOR_METHOD,
               linestyle="-", linewidth=1.0, alpha=0.6, label="target = -20")
    for i, v in zip(xs, EVAL_METRICS):
        # The bar goes from 0 down to v=-20. Place the value label
        # in the middle of the bar, rotated for readability.
        ax.text(i, v / 2, f"{v:.1f}", ha="center", va="center",
                fontsize=13, color="white", fontweight="medium")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"seed {s}" for s in EVAL_SEEDS])
    ax.set_ylabel("metric (lower = better)")
    ax.set_title("Final evaluator results")
    ax.set_ylim(-22.5, 2.0)
    ax.legend(loc="upper right", frameon=False)
    ax.text(-0.14, 1.05, "(c)", transform=ax.transAxes,
            fontsize=14, fontweight="bold")

    fig.suptitle(
        "SA solver: 20-point solution reached in 8.5 s (2 of 8 seeds); "
        "evaluator metric = −20 across all 3 seeds",
        fontsize=13, fontweight="medium", y=1.04,
    )
    out = FIGS / "results.png"
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    make_narrative()
    make_results()
