#!/usr/bin/env python3
"""
FIG-1: k-resolved dissociation map for Paper 5a.
Plots the fixed-axis Leggett-Garg witness K3(Q=z) against the stabilizer
2-Renyi magic M2 (log2 basis) across k=2..10, from the verified dissociation
table (p5a_r1_engine.py + p5a_r2_dissociation_table.py, cross-checked
against p5a_theorem_b_spot_checks_gram.py). Values hard-coded here from the
already-verified table (Table I of the paper) -- this script only renders
the figure, it does not recompute the physics.
Strips matplotlib metadata on save (no author/host leakage), >=300 dpi.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

k = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])
K3_fixedaxis = np.array([1.000, 1.498, 1.000, 1.479, 1.495, 1.494, 1.342, 1.497, 1.491])
K3_optimized = np.array([1.000, 1.500, 1.000, 1.500, 1.500, 1.498, 1.427, 1.500, 1.500])
M2 = np.array([0.000, 0.585, 0.5585, 0.585, 0.584, 0.585, 0.5530, 0.585, 0.585])

# Drawn at the width it is printed at (revtex \columnwidth is about 3.4 in), so that the font
# sizes below are the sizes on paper. The earlier 4.6 in version was scaled down by
# \includegraphics[width=\columnwidth] by a factor 0.739, which put EVERY label below 6 pt
# (8.0 -> 5.91, 6.5 -> 4.80). LaTeX reports nothing, because nothing overflows.
fig, ax1 = plt.subplots(figsize=(3.4, 2.6))

ax1.plot(k, K3_fixedaxis, "o-", color="#1f5fa8", label=r"$K_3(Q{=}\hat z)$", zorder=3)
ax1.scatter([8], [K3_optimized[6]], marker="^", s=40, color="#7fa8d9",
            label=r"$K_3^{\mathrm{opt}}$ ($k{=}8$ only)", zorder=4)
ax1.axhline(1.0, color="gray", lw=0.6, ls=":")
ax1.axhline(1.5, color="gray", lw=0.6, ls=":")
ax1.set_ylabel(r"Leggett–Garg witness $K_3$", color="#1f5fa8", fontsize=8)
ax1.tick_params(axis="y", labelcolor="#1f5fa8", labelsize=7)
ax1.tick_params(axis="x", labelsize=7)
ax1.set_xlabel(r"$k$ (SU(2)$_k$ level)", fontsize=8)
ax1.set_xticks(k)
ax1.set_ylim(0.9, 1.6)

ax2 = ax1.twinx()
ax2.plot(k, M2, "s--", color="#c9622a", label=r"$M_2$ (state magic, $\log_2$)", zorder=3)
ax2.axhline(np.log2(1.5), color="#c9622a", lw=0.6, ls=":")
ax2.set_ylabel(r"stabilizer 2-Rényi magic $M_2$", color="#c9622a", fontsize=8)
ax2.tick_params(axis="y", labelcolor="#c9622a", labelsize=7)
ax2.set_ylim(-0.03, 0.62)

# k=4 blind-spot annotation
ax1.scatter([4], [1.000], s=90, facecolors="none", edgecolors="black", linewidths=1.3, zorder=5)
ax1.annotate("k=4\nblind", xy=(4, 1.0), xytext=(4.5, 1.12),
             fontsize=7, ha="left",
             arrowprops=dict(arrowstyle="-", lw=0.6))

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=6.5, loc="lower right", framealpha=0.9)

fig.tight_layout()
if __name__ == "__main__":
    fig.savefig(Path(__file__).resolve().parent / "p5a_fig1_dissociation_map.png",
                dpi=300, metadata={"Author": None, "Software": None,
                                   "Title": None, "Description": None})
    print("wrote p5a_fig1_dissociation_map.png")
