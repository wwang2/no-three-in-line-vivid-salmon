"""Generate figures/narrative.png and figures/results.png for orbit 01.

narrative.png — the final 20-point placement on the 10x10 grid, shown
                beside a naive baseline (diagonal, invalid) and a
                16-point greedy, so the reader SEES the improvement.

results.png   — quantitative summary: metric across strategies, and a
                bar chart of greedy outcome distribution.
"""
from __future__ import annotations

import os
import random
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import search  # noqa: E402
from solution import POINTS as FINAL_POINTS  # noqa: E402

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

N = 10
FIG_DIR = HERE / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {
    "final": "#4C72B0",      # deep blue — our 20-point solution
    "greedy": "#DD8452",     # warm orange — intermediate greedy
    "baseline": "#888888",   # gray — diagonal baseline (invalid)
    "grid": "#cccccc",
    "line": "#C44E52",       # collinearity highlight red
    "target": "#55A868",     # green target line
}


# --------------------------------------------------------------------------- #
def draw_grid(ax, points, title, color, marker_size=180, show_lines=None):
    """Draw an NxN grid with scatter of `points`."""
    ax.set_xlim(-0.6, N - 0.4)
    ax.set_ylim(-0.6, N - 0.4)
    ax.set_aspect("equal")
    ax.set_xticks(range(N))
    ax.set_yticks(range(N))
    ax.grid(True, color=COLORS["grid"], linewidth=0.6, alpha=0.7)
    ax.tick_params(length=0)
    ax.set_title(title)

    # empty cells as faint dots
    all_cells = [(r, c) for r in range(N) for c in range(N)]
    occupied = set(points)
    empty = [p for p in all_cells if p not in occupied]
    if empty:
        er = [p[1] for p in empty]
        ec = [p[0] for p in empty]
        ax.scatter(er, ec, s=12, c=COLORS["grid"], marker="o", zorder=1)

    # collinear lines (if any) drawn first so points sit on top
    if show_lines:
        for line_idxs in show_lines:
            xs = [points[i][1] for i in line_idxs]
            ys = [points[i][0] for i in line_idxs]
            # draw across the grid along the line
            if len(xs) >= 2:
                ax.plot(xs, ys, color=COLORS["line"], linewidth=1.2,
                        alpha=0.8, zorder=2)

    # chosen points
    if points:
        r = [p[1] for p in points]
        c = [p[0] for p in points]
        ax.scatter(r, c, s=marker_size, c=color, edgecolors="white",
                   linewidth=1.2, zorder=3)

    ax.set_xlabel("column (c)")
    ax.set_ylabel("row (r)")
    ax.invert_yaxis()  # row 0 on top, like a matrix


def make_narrative():
    # grow a 16-point greedy for contrast, and a 10-diagonal baseline
    baseline = [(i, i) for i in range(N)]   # invalid (collinear)
    # grow greedy: find a 16-point valid subset
    greedy16 = search.greedy_valid_subset(
        [(r, c) for r in range(N) for c in range(N)], seed=0
    )
    # trim to exactly 16 if longer
    greedy16 = greedy16[:16]

    final = list(FINAL_POINTS)

    # also regenerate a fresh greedy with more luck to fill to ~16-17
    best_mid = greedy16
    for s in range(200):
        cand = search.greedy_valid_subset(
            [(r, c) for r in range(N) for c in range(N)], seed=s
        )
        if len(cand) > len(best_mid):
            best_mid = cand
        if len(best_mid) >= 17:
            break
    greedy16 = best_mid

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.2), sharex=True, sharey=True)

    # (a) baseline: diagonal
    draw_grid(
        axes[0], baseline,
        title=f"Diagonal baseline — {len(baseline)} pts  (INVALID)",
        color=COLORS["baseline"], marker_size=150,
        show_lines=[list(range(len(baseline)))],
    )
    # (b) greedy
    draw_grid(
        axes[1], greedy16,
        title=f"Random greedy — {len(greedy16)} pts  (valid, sub-optimal)",
        color=COLORS["greedy"], marker_size=160,
    )
    # (c) final 20
    draw_grid(
        axes[2], final,
        title=f"Algebraic + SA — {len(final)} pts  (Erdős–Szekeres 2N optimum)",
        color=COLORS["final"], marker_size=180,
    )

    # panel labels ABOVE the title area on the left margin
    for ax, lab in zip(axes, ["(a)", "(b)", "(c)"]):
        ax.text(0.02, 1.10, lab, transform=ax.transAxes,
                fontsize=15, fontweight="bold", color="#333")

    fig.suptitle(
        "No-three-in-line on a 10×10 grid — diagonal baseline -> greedy -> optimal",
        fontsize=15, fontweight="medium", y=1.05,
    )

    fig.savefig(FIG_DIR / "narrative.png", dpi=200, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print("wrote", FIG_DIR / "narrative.png")


def _collect_sa_traces():
    """Run SA a few times with different init conditions, capture energy traces.

    Returns a list of (label, color, trace, final_E) tuples where
    trace is a list of (step, best_E_so_far) points.
    """
    traces = []

    # (1) Warm-start from greedy-19 (seed 3089), SA seed 0 -> reaches E=0 fast
    greedy19 = search.greedy_valid_grid(seed=3089)[:19]
    _, E, trace = search.simulated_anneal(
        K=20, seed=0, init=greedy19, iters=100_000,
        T0=2.0, T_end=0.01, record_trace=True,
    )
    traces.append(("warm-start (greedy-19, SA seed 0)  -> E=0",
                   COLORS["final"], trace, E))

    # (2) Cold K=20 seed=5 -> gets stuck at E=1
    _, E, trace = search.simulated_anneal(
        K=20, seed=5, init=None, iters=100_000,
        T0=2.0, T_end=0.01, record_trace=True,
    )
    traces.append(("cold start (SA seed 5)  -> stuck at E=1",
                   COLORS["greedy"], trace, E))

    # (3) Cold K=20 seed=12 -> also stuck
    _, E, trace = search.simulated_anneal(
        K=20, seed=12, init=None, iters=100_000,
        T0=2.0, T_end=0.01, record_trace=True,
    )
    traces.append(("cold start (SA seed 12)  -> stuck at E=1",
                   COLORS["baseline"], trace, E))

    return traces


def make_results():
    """Quantitative panel: strategy bars + greedy-size distribution + SA convergence."""
    fig = plt.figure(figsize=(18.0, 6.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.2, 1.2])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    # --- (a) bar chart of strategies vs point count
    strategies = [
        ("Diagonal\nbaseline", 0, COLORS["baseline"], "invalid"),
        ("Algebraic\nunion (raw)", 11, COLORS["baseline"], ""),
        ("Random greedy\n(single seed)", 16, COLORS["greedy"], ""),
        ("Random greedy\n(5000 restarts)", 16, COLORS["greedy"], "mode"),
        ("Algebraic + SA\n(warm-start K=20)", 20, COLORS["final"], ""),
    ]
    labels = [s[0] for s in strategies]
    counts = [s[1] for s in strategies]
    bar_colors = [s[2] for s in strategies]
    bars = ax1.bar(range(len(strategies)), counts, color=bar_colors,
                   edgecolor="white", linewidth=1.5, zorder=2)
    for i, (bar, s) in enumerate(zip(bars, strategies)):
        v = s[1]
        note = s[3]
        if note:
            lab = f"{v} ({note})"
        else:
            lab = str(v)
        ax1.text(bar.get_x() + bar.get_width() / 2, v + 0.6, lab,
                 ha="center", va="bottom", fontsize=11,
                 color=s[2] if v > 0 else COLORS["baseline"])

    # mark the max hit by 5000-restart greedy with a small cap line above its bar
    ax1.plot([3 - 0.28, 3 + 0.28], [19, 19],
             color="#333", linewidth=1.6, zorder=3)
    ax1.text(3, 19.4, "max=19", ha="center", va="bottom", fontsize=9,
             color="#333")

    ax1.axhline(20, color=COLORS["target"], linestyle="--", linewidth=1.2,
                zorder=1)
    # put target text on the upper-left side, clear of bars
    ax1.text(0.02, 22.6, "target = 20 (Erdős–Szekeres 2N optimum)",
             color=COLORS["target"], fontsize=10, ha="left")
    ax1.set_xticks(range(len(strategies)))
    ax1.set_xticklabels(labels, fontsize=9)
    ax1.set_ylabel("# valid points placed")
    ax1.set_ylim(0, 25)
    ax1.set_title("Strategy comparison on 10×10 grid", pad=28)
    ax1.grid(True, axis="y", alpha=0.2)
    ax1.text(0.0, 1.12, "(a)", transform=ax1.transAxes,
             fontsize=15, fontweight="bold", color="#333")

    # --- (b) greedy size distribution (5000 restarts)
    hist = {12: 2, 13: 38, 14: 454, 15: 1653, 16: 1980, 17: 801, 18: 70, 19: 2}
    xs = sorted(hist.keys())
    ys = [hist[k] for k in xs]
    ax2.bar(xs, ys, color=COLORS["greedy"], alpha=0.85,
            edgecolor="white", linewidth=1.0)
    ax2.axvline(20, color=COLORS["final"], linestyle="-", linewidth=2.2,
                zorder=3, label="SA warm-start -> 20")
    ax2.axvline(19, color=COLORS["target"], linestyle=":", linewidth=1.5,
                zorder=3, label="greedy max (2 of 5000)")
    for x, y in zip(xs, ys):
        ax2.text(x, y + 45, str(y), ha="center", fontsize=9,
                 color=COLORS["greedy"])
    ax2.set_xlabel("# valid points found by random greedy")
    ax2.set_ylabel("frequency over 5000 restarts")
    ax2.set_xticks(list(range(12, 21)))
    ax2.set_ylim(0, 2300)
    ax2.set_title("Greedy outcome distribution", pad=28)
    ax2.legend(loc="upper right")
    ax2.grid(True, axis="y", alpha=0.2)
    ax2.text(0.0, 1.12, "(b)", transform=ax2.transAxes,
             fontsize=15, fontweight="bold", color="#333")

    # --- (c) SA convergence curves (energy vs step)
    traces = _collect_sa_traces()
    for (label, color, trace, final_E) in traces:
        steps = [t[0] for t in trace]
        # +0.5 shift so E=0 is visible on log-y
        energies = [t[1] + 0.5 for t in trace]
        lw = 2.4 if final_E == 0 else 1.6
        alpha = 1.0 if final_E == 0 else 0.75
        ax3.plot(steps, energies, color=color, linewidth=lw,
                 alpha=alpha, label=label)

    # horizontal guide at E=0 (shown at 0.5 because of the shift)
    ax3.axhline(0.5, color=COLORS["target"], linestyle="--",
                linewidth=1.2, zorder=1)
    ax3.text(5e4, 0.55, "E = 0  (no three collinear)",
             color=COLORS["target"], fontsize=9, va="bottom", ha="center")

    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlim(1, 1.2e5)
    ax3.set_ylim(0.4, 80)
    ax3.set_xlabel("SA step")
    ax3.set_ylabel("best triple-count energy E  (shift +0.5)")
    ax3.set_title("SA energy descent — warm-start vs cold-start", pad=28)
    ax3.grid(True, which="both", alpha=0.2)
    ax3.legend(loc="upper right", fontsize=8.5)
    ax3.text(0.0, 1.12, "(c)", transform=ax3.transAxes,
             fontsize=15, fontweight="bold", color="#333")

    fig.suptitle(
        "Orbit 01 results — algebraic-seed + SA hits the 2N ceiling (metric = -20)",
        fontsize=14, fontweight="medium", y=1.03,
    )
    fig.savefig(FIG_DIR / "results.png", dpi=200, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print("wrote", FIG_DIR / "results.png")


if __name__ == "__main__":
    make_narrative()
    make_results()
