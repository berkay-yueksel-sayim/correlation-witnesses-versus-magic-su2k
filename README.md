# Correlation Witnesses versus Magic: A k-Resolved Nonstabilizerness Map for SU(2)_k Anyon Fusion Spaces

**Author:** Berkay Yüksel Sayim
**ORCID:** [0009-0004-4993-7352](https://orcid.org/0009-0004-4993-7352)
**DOI (all versions):** [10.5281/zenodo.21134454](https://doi.org/10.5281/zenodo.21134454)

## Abstract

Standard correlation witnesses — spatial Bell/CHSH, temporal Leggett-Garg
K3, and KCBS contextuality — are known to detect nonstabilizerness
("magic") in some settings. We report a worked example, a k-resolved
nonstabilizerness map of the SU(2)_k anyon braid-representation family, in
which a standard temporal witness is instead systematically blind: at k=4
the three-time Leggett-Garg witness saturates its macrorealistic bound
exactly (K3=1.000000 over every state, braid element, and measurement
axis, in the equal-time-step protocol with V1=V2) while the single-qubit
fusion channel carries near-maximal nonstabilizerness (M2=0.5585, 95% of
the finite-dimensional ceiling log2(3/2)). We prove this blindness as a
structural theorem: the witness depends only on the Bloch-sphere Gram
geometry of the braid orbit, not on nonstabilizerness, and k=4 happens to
align the fixed measurement axis with a threefold orbit symmetry that
caps the witness at 1; a finite, Clifford-generating group — the chiral
octahedral group, the k=2 braid image — with a misaligned axis reaches
K3=3/2 under a two-propagator protocol. We complement this with three
independent certificates of genuine nonstabilizerness — two of them
long-range, the third a complementary gate-based non-Cliffordness
measure, an honest negative control, a hardware-oriented prediction,
and an interface with KCBS contextuality at d=3.

## Contents

This record is compiled from `main_v1.1.tex` (RevTeX 4-2). `p5a_fig1_dissociation_map.png`
is Fig. 1 of the paper.

**Reproduction code and deposited data** (all deterministic, NumPy only
unless noted; each script was re-run standalone from this folder to confirm
it reproduces the value cited in the paper text):

- `p5a_korbany_abstract_crosscheck.py` + `p5a_korbany_witness.json` — the
  doubled-Fibonacci mutual-information witness H=1.700979 (Sec. "Three
  certificates of genuine nonstabilizerness", item 1), reproduced from the
  quantum dimensions in closed form; two independent routes agree to
  |Δ|=2.2e-16 (see script output). Includes the Z2-toric-code null control
  (H=2.000000 exactly) and the toric-T-state anchor (H=0.600876).
- `p5a_r1_engine.py`, `p5a_r2_dissociation_table.py`,
  `p5a_r1r2_AUSWERTUNG.md` — the full k-resolved dissociation table
  (Table I / Fig. 1): three independent magic monotones (M2, RoM, geometric
  distance) and the Leggett-Garg witness across k=2..10, plus the McKay
  binary-polyhedral-group identification.
- `p5a_fig1_dissociation_map.py` — renders `p5a_fig1_dissociation_map.png` from
  the verified table values (metadata-stripped, 300 dpi).
- `p5a_theorem_a_spot_checks.py`, `p5a_theorem_b_spot_checks_gram.py`
  (+ dependency `p5a_theorem_find_counterexample.py`), `p5a_THEOREM_NOTE_k4_blindness.md`
  — the structural theorem of Sec. "Structural theorem": Theorem A
  (Pauli-aligned-Clifford special case) and Theorem B (general Bloch-Gram
  geometry, covers the real k=4 braid image), with the alignment
  counterexample (misaligned octahedral group reaching K3=3/2) and the
  operational two-qubit-unitary confirmation (Check G5).
- `p5a_n2_longrange_tsre.py` (+ dependency `p5a_n2_tsre_base.py`) +
  `p5a_n2_longrange_tsre.json` — the gauge-invariant minimum ground-space
  stabilizer Rényi entropy certificate (item 2): M2(Fib)=6.393 (Monte Carlo
  lower bound), M2(Toric)=0 exact null control, {0,±1}-fraction
  discriminator.
- `p5a_tsre_robust_m2_range.json` — the gauge-invariant ground-space M2 range
  [6.5, 8.9] over the full four-dimensional ground space (values-and-provenance
  record; the exact-diagonalization construction behind it is not part of this
  record, see Data and Code Availability).
- `p5a_gate_magic_t3.py` + `p5a_ergebnis_t3.json` — the amortized gate-based
  non-Cliffordness measure (item 3) across k=2,3,4,5,8, with the closed-form
  T-gate anchor M2^A(T)=2−log2(3)=0.4150 (Sec. "Hardware-oriented
  prediction").
- `p5a_n1_hardening.py` + `p5a_ergebnis_hardening_b.json` — the honest negative
  control (Sec. "Honest negative control"): leakage-free strict additivity,
  and the d=5 Wigner-Mana sanity values across k=3,4,5,8 (log2 basis).
- `p5a_t5_kcbs_magic.py` + `p5a_ergebnis_t5.json` — the KCBS-contextuality/magic
  interface at d=3 (Sec. "Interface with contextuality at d=3").
- `p5a_mckay_four_axis_map.py` + `p5a_ergebnis_t4.json` — a fourth independent
  cross-check of the k=4 dissociation (LGI, magic, gate-magic, McKay
  correspondence together).
- `p5a_fig2_orbit_geometry.py` + `p5a_fig2_orbit_geometry.png` — renders Fig. 2
  (orbit geometry: group order, K3 and M2 per k) from the verified values in
  `p5a_ergebnis_t4.json` (metadata-stripped, 300 dpi).
- `p5a_korbany_lattice_route.py` + `p5a_korbany_lattice_vacuum.json` — the
  genuinely independent second route for the Korbany witness H (item 1 above):
  the vacuum column of the modular S-matrix extracted from a 3x3 lattice via the
  rotation channel, with no analytic quantum dimensions entering. The two lattice
  inputs of that extraction live in the companion record; the JSON names them.
- `p5a_n2_hex_cycle_generator.py` — geometric prerequisite of
  `p5a_n2_tsre_base.py`: the ordered hexagon cycle (A-B-A-B-A-B, six distinct
  edges) on an N x M honeycomb torus.
- `p5a_V_gs_2x2.npy` — the doubled-Fibonacci ground space of the 2x2 honeycomb
  torus, (175, 4) float64, built by `fib_groundstate` in
  `p5a_n2_tsre_base.py`; see the note below.

All scripts were verified to reproduce the values cited in the paper text to
the stated precision when re-run from this folder (2026-07-13).

## The deposited ground space

This file carries the ground space of the restricted sector in the basis the
generator produces; the export is deterministic — two separate processes give
byte-identical output. It spans the same subspace as `p6_V_gs_2x2.npy` in the
companion record (doi:10.5281/zenodo.21362245) to machine precision — the
projectors agree to 9.4e-17, while the arrays themselves differ
(max |A−B| = 2.6e-1) — because the basis within that subspace differs: this
file carries the raw output, the companion record a vacuum-adapted basis. Taken
together the two files are the check rather than a claim: the five distinguished
rows have rank 1 in both — that is the physics — while the occupancy pattern
differs, [175,175,175,175] here against [170,170,175,170] there — that is the
convention. The structure of those rows, and the lattice geometry behind it, are
documented in the companion record.

## License
- Paper, figures, and data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — see `LICENSE`
- Source code (`*.py`): [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0) — see `LICENSE-CODE`
