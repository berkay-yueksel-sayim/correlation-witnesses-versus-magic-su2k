"""
p5a_korbany_lattice_route.py -- genuinely independent second route for the
doubled-Fibonacci Korbany mutual-information witness H (Sec. "Three certificates
of genuine nonstabilizerness", item 1).

The primary (direct) route computes H = -sum_a p_a log2 p_a from the analytic
quantum dimensions d = (1, phi, phi, phi^2), p_a = d_a^2 / D^2. The closed-form
rearrangement H = log2(D^2) - sum_a p_a log2(d_a^2) is algebraically identical to
that sum (substitute p_a = d_a^2/D^2) and is therefore NOT an independent route --
it only re-measures floating-point noise of the same computation.

The independent route used here reads the topological-sector weights |S_{0a}|^2
from the VACUUM COLUMN of the doubled-Fibonacci modular S-matrix that was
extracted from a 3x3 lattice (Wilson-loop + C3-rotation channel; see
p5a_korbany_lattice_vacuum.json for provenance and input md5s). No analytic quantum
dimensions and no F/R symbols enter the lattice extraction, so its vacuum column
is an independent measurement of the same weights. The two entropies agree to the
lattice-extraction residual (|Delta| ~ 1.8e-11), the value quoted in the paper.

Deterministic; reads only the shipped JSON. No lattice re-run required.
"""
import json
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def shannon_bits(weights):
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    return float(-np.sum(w * np.log2(w)))


def main():
    # --- analytic (direct) route ---
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    d = np.array([1.0, phi, phi, phi ** 2])
    p = d ** 2 / np.sum(d ** 2)
    H_direct = shannon_bits(p)

    # --- independent route: lattice-extracted S vacuum column ---
    with open(os.path.join(HERE, "p5a_korbany_lattice_vacuum.json"), encoding="utf-8") as f:
        rec = json.load(f)
    w_lat = rec["lattice_vacuum_column_abs2"]
    H_lat = shannon_bits(w_lat)

    delta = abs(H_lat - H_direct)
    print("analytic direct route : H = %.12f" % H_direct)
    print("lattice S-column route: H = %.12f" % H_lat)
    print("|Delta|               = %.3e  (lattice-extraction residual)" % delta)
    print("provenance run        : %s (center %s)"
          % (rec["provenance"].get("quantity", "?")[:60], rec["provenance"]["center"]))

    # sanity: the lattice weights must match the analytic weights to the residual
    assert abs(H_direct - 1.700979) < 1e-5, "analytic H off"
    assert delta < 1e-8, "lattice route does not reproduce H within the expected residual"
    print("OK: independent lattice route reproduces H to the extraction residual.")


if __name__ == "__main__":
    main()
