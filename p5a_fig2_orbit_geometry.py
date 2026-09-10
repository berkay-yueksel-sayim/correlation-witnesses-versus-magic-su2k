#!/usr/bin/env python3
"""
FIG-2: what does NOT decide blindness of the fixed-axis Leggett-Garg witness.

Every plotted number is READ from p5a_ergebnis_t4.json (blocks `map` and `curve`);
nothing is typed here. That is the difference to fig1, which renders the already
verified Table I. The Bloch dot products quoted in the annotation are exact
constants taken from p5a_THEOREM_NOTE_k4_blindness.md (Theorem B, checks G1/G2/G3):
tetrahedron {-1/3}, octahedron {0,+-1}, icosahedron {+-1,+-1/sqrt5}.

Statement of the figure: k=2 (2O) and k=4 (2T) are blind, k=8 (2I) is not, although
all three have a finite braid image; the group orders are not ordered along k
(48 -> 24 -> 120); and k=4 and k=8 carry nearly the same magic while sitting on
opposite sides of the witness. So neither finiteness, nor group order, nor the
amount of magic decides. What does decide -- the orbit Gram structure -- is stated
in the text, not drawn here: no column of the JSON carries it as a scalar.

The two dense levels k=3 and k=5 are shown on purpose. They are the evidence for
their own exclusion: among the three finite levels alone, the gate-magic column
would look like a clean discriminator (0.000/0.765 blind against 0.815 seeing),
and it is exactly k=3 (0.562) that breaks it.

Conventions: group orders are printed in the BINARY convention used by Table I
(48/24/120). The JSON stores the PROJECTIVE order (24/12/60); the doubling is
performed below as an explicit expression, not typed in.

Strips matplotlib metadata on save (no author/host leakage), >=300 dpi.
"""
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if __name__ == "__main__" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
with open(HERE / "p5a_ergebnis_t4.json", encoding="utf-8") as fh:
    data = json.load(fh)

# `map` carries K3/M2/group, `curve` carries order/dense -- joined on k.
map_by_k = {e["k"]: e for e in data["map"]}
curve_by_k = {e["k"]: e for e in data["curve"] if e["k"] in map_by_k}

TOL = 1e-9
points = []
for k in sorted(map_by_k):
    m, c = map_by_k[k], curve_by_k[k]
    ordnung_projektiv = c["order"]                     # 24 / 12 / 60, None for dense
    # Table I prints the BINARY order; the JSON stores the projective one.
    ordnung_binaer = 2 * ordnung_projektiv if ordnung_projektiv is not None else None
    # Cross-check: the order is in the JSON twice -- as `curve.order` and inside
    # the `map.group` string. If the two ever drift apart, stop instead of drawing.
    if ordnung_projektiv is not None:
        im_label = re.search(r"\((\d+)\)", m["group"])
        assert im_label and int(im_label.group(1)) == ordnung_projektiv, (
            f"k={k}: curve.order={ordnung_projektiv} vs map.group={m['group']}")
    if c["dense"]:
        zustand = "dense"                              # sampling, not conclusive
    elif m["K3"] <= 1.0 + TOL:
        zustand = "blind"
    else:
        zustand = "sees"
    points.append({"k": k, "K3": m["K3"], "M2": m["M2"], "zustand": zustand,
                   "name": m["group"].split("(")[0], "ordnung": ordnung_binaer})

STIL = {                       # verdict lives in marker+colour, K3 lives in position
    "blind": dict(marker="o", color="#1f5fa8", label="blind (established)"),
    "sees":  dict(marker="^", color="#c9622a", label="witness fires (established)"),
    "dense": dict(marker="o", color="#8a8a8a", label="dense image (not conclusive)"),
}

# Drawn at the width it is printed at (revtex \columnwidth is about 3.4 in), so that the
# font sizes below are the sizes on paper. Rendering wider and letting \includegraphics
# scale down silently shrinks every label -- LaTeX reports nothing, because nothing overflows.
fig, ax = plt.subplots(figsize=(3.4, 2.75))
ax.axhline(1.0, color="gray", lw=0.7, ls=":")
# left of the data: the lower right corner is taken by the legend (checked by looking
# at the rendered PNG -- the first layout hid this label underneath it)
ax.text(1.62, 1.023, "macrorealist bound", fontsize=6.2, color="gray", ha="left")

gezeigt = set()
for p in points:
    s = STIL[p["zustand"]]
    ax.scatter(p["k"], p["K3"], marker=s["marker"], s=42, zorder=4,
               color=s["color"] if p["zustand"] != "dense" else "none",
               edgecolors=s["color"], linewidths=1.1,
               label=s["label"] if p["zustand"] not in gezeigt else None)
    gezeigt.add(p["zustand"])
    if p["ordnung"] is not None:                        # finite: name the group
        ax.annotate(f"$|{p['name']}|$ = {p['ordnung']}", xy=(p["k"], p["K3"]),
                    xytext=(0, -11), textcoords="offset points",
                    fontsize=6.0, ha="center", color="#333333")

# The twin pair: nearly the same magic, opposite sides of the witness.
vier = next(p for p in points if p["k"] == 4)
acht = next(p for p in points if p["k"] == 8)
ax.annotate("", xy=(acht["k"], acht["K3"]), xytext=(vier["k"], vier["K3"]),
            arrowprops=dict(arrowstyle="<->", lw=0.8, color="#444444",
                            connectionstyle="arc3,rad=-0.32"))
# below the arc, not across it (the first layout let the arrow run through the text)
ax.text(6.15, 1.115,
        f"same magic\n$M_2$ = {vier['M2']:.4f} vs {acht['M2']:.4f}",
        fontsize=6.0, ha="center", va="bottom", color="#444444")

# lower left, with the legend moved to the upper right: at >=6 pt the legend box grew and
# covered the |2T| = 24 label and the end of "macrorealist bound" (seen in the rendered PNG)
ax.text(1.55, 1.40,
        "a finite braid image\n"
        "does not imply blindness:\n"
        "neither finiteness,\n"
        "nor group order,\n"
        "nor magic decides",
        fontsize=6.2, ha="left", va="top", color="#222222", linespacing=1.35)

ax.set_xlabel(r"$k$ (SU(2)$_k$ level)", fontsize=8)
ax.set_ylabel(r"Leggett–Garg witness $K_3(Q{=}\hat z)$", fontsize=8)
ax.set_xticks(sorted(map_by_k))
ax.tick_params(labelsize=7)
ax.set_xlim(1.3, 8.85)      # left margin: the |2O| = 48 label at k=2 was clipped otherwise
ax.set_ylim(0.93, 1.70)
ax.legend(fontsize=6.2, loc="upper right", framealpha=0.95,
          borderpad=0.4, handletextpad=0.5, labelspacing=0.35)
fig.tight_layout()

if __name__ == "__main__":
    target = HERE / "p5a_fig2_orbit_geometry.png"
    fig.savefig(target, dpi=300, metadata={"Author": None, "Software": None,
                                         "Title": None, "Description": None})
    for p in points:
        print(f"k={p['k']}  {p['name']:5s} order={p['ordnung']}  "
              f"K3={p['K3']:.6f}  M2={p['M2']:.6f}  -> {p['zustand']}")
    print(f"wrote {target.name}")
