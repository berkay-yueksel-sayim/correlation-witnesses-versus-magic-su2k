# Proof note: blindness of the dichotomic Lüders LGI witness under finite braid images — alignment geometry (orbit Gram), with a Pauli-aligned Clifford special case

**Run:** 2026-07-07_B_k4blind_r1 (task B2), repair 2026-07-08 · **Status:** two theorems, every lemma numerically spot-checked.
**Theorem A** (Pauli-aligned Clifford special case): `spot_checks.py` (seed=42, deterministic). **Theorem B** (orbit Gram geometry, covers the real k=4): `spot_checks_gram.py` (deterministic, reproduces the A-p07 reference reimplementation). Raw output of both at the end of this note.
**Wording:** What is sayable is "this witness is structurally blind because the Q axis and the orbit are aligned" — no more, no priority claim. Results by others are attributed (Jones; Tuba–Wenzl; Freedman–Larsen–Wang; Budroni–Emary).

---

## Setup and definitions

- Qubit, fixed dichotomic Q = n₀·σ (standard frame: Q = Pauli Z, n₀ = e_z). Three measurement times t1 < t2 < t3; propagation between the times by elements of a group G ⊂ SU(2): V1 (t1→t2), V2 (t2→t3). Heisenberg observables:
  Q1 = Q, Q2 = V1† Q V1, Q3 = (V2V1)† Q (V2V1). Their Bloch directions n1 = n₀, n2, n3 lie in the orbit of n₀ under the SO(3) image of G; **n2 and n3 can be chosen freely and independently, because V1 and V2 are independent.**
- **Lüders correlator** for dichotomic projective measurements (standard convention, cf. Fritz 2010; Budroni–Emary arXiv:1309.3678):
  C_ij = Σ_{a,b=±1} a·b · P(a at t_j, b at t_i) = (1/2) tr(ρ {Q_i, Q_j}).
- **Witness:** K3 = C21 + C32 − C31; macrorealist bound K3 ≤ 1.

---

# THEOREM A (Pauli-aligned-Clifford-Spezialfall)

Let **G be a subgroup of the single-qubit Clifford group in the frame of Q** — i.e. G normalizes the Pauli group spanned by Q, equivalently: Q lies in the signed Pauli axis grid ±{X, Y, Z} and the SO(3) image of G permutes that grid. Then, for the dichotomic Lüders LGI witness with fixed Q = Pauli Z and propagators from G:

**K3 ≤ 1 for ALL ρ.** The witness is structurally blind, independently of the state.

> **⚠ WARNING (core of the 2026-07-08 repair):** The **real** k=4 braid image of SU(2)_k=4 does **NOT** satisfy this hypothesis **in the fusion-Z frame.** There the SO(3) image of σ₁ is a projective **3-fold rotation about z** (cos θ = −1/2), so z is a 3-fold axis of the tetrahedral group T, not a Pauli 2-fold axis; orbit(Z) consists of **4 tetrahedron directions** (z + 3 directions with z component −1/3, pairwise Bloch dots exactly −1/3), NOT of ±{X, Y, Z}. In particular σ₂† Z σ₂ = Bloch(−√2/3, −√6/3, −1/3) is **not** a ±Pauli. The k=4 group is indeed conjugate to a Clifford subgroup, but the conjugation does not fix Z — Lemma 1 breaks. **The k=4 blindness nevertheless holds, but via the orbit Gram geometry (THEOREM B), not via THEOREM A.** (The first version of this note wrongly claimed k=4 ⊂ Clifford in the Z frame; see the honesty block.)

## Lemma A1 (Clifford conjugation keeps Z in the signed Pauli orbit)

**Claim:** U Clifford (in the Q frame) ⇒ U†ZU = ±P for some P ∈ {X, Y, Z}.

**Proof:** By definition the Clifford group is the normalizer of the Pauli group: U†ZU lies in the Pauli group, hence U†ZU = ω·P' with a phase ω, P' ∈ {I, X, Y, Z}. Hermiticity of U†ZU forces ω = ±1; tracelessness under conjugation (tr(U†ZU) = tr Z = 0) rules out P' = I. Hence U†ZU = ±P, P ∈ {X, Y, Z}. Since G is a subgroup, V1 and V2V1 both lie in G ⊂ Clifford ⇒ **all three Heisenberg observables Z1, Z2, Z3 are signed Paulis.** ∎

**Spot check 1:** All 24 elements of 2T (8 Lipschitz + 16 Hurwitz unit quaternions → SU(2)) built explicitly, unitarity + det=1 + group closure (24×24 products) verified; U†ZU ∈ {±X, ±Y, ±Z} for **24/24** elements, distribution uniform at 4 per target. **Note:** 2T appears here as an **abstract** Pauli-aligned Clifford subgroup (axes in the Pauli grid) — that is a DIFFERENT embedding from the real k=4 fusion image (see the warning / THEOREM B). See raw output CHECK 1.

## Lemma A2 (Lüders correlators are discrete and ρ-independent)

**Claim:** For signed Paulis A, B one has (1/2){A, B} = ±I if B = ±A, and 0 otherwise. Consequently C_ij = (1/2) tr(ρ {Z_i, Z_j}) is either ±1 (**ρ-independent**) or 0.

**Proof:** In d=2 there are only two cases for two signed Paulis: (i) same Pauli axis, B = ±A ⇒ (1/2){A, B} = ±A² = ±I; (ii) different axes ⇒ the Paulis anticommute, {A, B} = 0. Substituting: C_ij = ±tr(ρ) = ±1 resp. C_ij = 0 — ρ drops out completely. ∎

**Spot check 2:** All 36 signed Pauli pairs: 12× (1/2){A,B} = ±I (same axis), 24× = 0 (different axis), tolerance 1e-12. Correlator over 206 random states: maximal spread 4.44e-16 (ρ-independent), all values in {−1, 0, +1}. See raw output CHECK 2.

## Lemma A3 (sign consistency: all cases end at K3 ≤ 1)

Write Z2 = s2·P2, Z3 = s3·P3 (Lemma A1). From Lemma A2 it follows exactly that:
C21 = s2·δ(P2=Z), C32 = s2·s3·δ(P3=P2), C31 = s3·δ(P3=Z). Complete case list (36 combinations):

| Case | Condition | C21 | C32 | C31 | K3 = C21+C32−C31 | # |
|---|---|---|---|---|---|---|
| A | P2 = Z, P3 = Z | e1 := s2 | e1e2 (e2 := s3s2) | e1e2·e1 = s3 | e1 + e2 − e1e2 ∈ {1, −3} | 4 |
| B | P2 = Z, P3 ≠ Z | ±1 | 0 | 0 | ±1 | 8 |
| C | P2 ≠ Z, P3 = P2 | 0 | ±1 | 0 (since P3 = P2 ≠ Z) | ±1 | 8 |
| D | P2 ≠ Z, P3 = Z | 0 | 0 | ±1 | ∓1 | 8 |
| E | P2 ≠ Z, P3 ∉ {Z, P2} | 0 | 0 | 0 | 0 | 8 |

**Case A written out** (the only one with three nontrivial correlators): Z2 = e1·Z, Z3 = e2·Z2 ⇒ Z3 = e1e2·Z ⇒ C31 = C21·C32, hence K3 = e1 + e2 − e1e2. Values: (+,+) → 1, (+,−) → 1, (−,+) → 1, (−,−) → −3. Maximum 1.
**Cases B–E:** As soon as at least one correlator is 0, the (anti)commutation relations propagate and at most one term of K3 is ±1. **None of the 36 combinations reaches K3 > 1.** ∎

**Spot check 3:** All 36 combinations: case formula against brute force over 56 random ρ, deviation < 1e-10; K3 value sets per case exactly as in the table; max K3 = 1. See raw output CHECK 3.

## Proof of Theorem A

Lemma A1 (with subgroup closure: V1, V2V1 ∈ G) ⇒ Z1, Z2, Z3 are signed Paulis. Lemma A2 ⇒ each of the three Lüders correlators lies, independently of ρ, in {−1, 0, +1}. Lemma A3 ⇒ each of the 36 configurations gives K3 ≤ 1. Since the correlators are ρ-independent, the bound holds for **all** ρ. ∎

**Spot check 4 (end to end):** All 576 propagator pairs (V1, V2) ∈ 2T × 2T, 14 random ρ each: Z2, Z3 always signed Paulis; **max K3 = 1.000000000000** — the bound is attained exactly (e.g. V1 = V2 = I), never exceeded. See raw output CHECK 4.

**Corollary A-a (k=2, abstractly Pauli-aligned):** The braid image at k=2 is binary octahedral 2O; **as an abstract group** 2O is a Pauli-aligned Clifford subgroup (SO(3) image = full octahedral group, permutes ±{X,Y,Z}) — Theorem A applies verbatim. **Spot check 1b:** all 48 elements of 2O conjugate Z to a ±Pauli (48/48), group closure 48×48 verified. (Whether the real k=2 fusion image carries this axis alignment is a matter of convention; THEOREM B G2 computes the k=2 orbit directly and likewise finds K3max = 1.)

**Spot check 5 (contrast, non-Clifford):** The non-Clifford precession U = exp(−iπX/6) violates Lemma A1 (U†ZU is not a ±Pauli) and yields K3 = 1.500000000000 for all ρ tested = the Lüders bound 3/2 (Budroni–Emary). See raw output CHECK 5.

---

# THEOREM B (orbit Gram geometry — covers the real k=4)

Let G ⊂ SU(2) be arbitrary and Q = n₀·σ fixed. Then the Lüders correlators are **Bloch inner products of the Heisenberg directions, and ρ-independent:**

**C_ij = n_i · n_j for EVERY ρ**, where n_i ∈ orbit(n₀) under the SO(3) image of G, with n1 = n₀ fixed and n2, n3 free and independent (V1, V2 independent).

**Proof:** For single-qubit dichotomics {n·σ, m·σ} = 2(n·m)·I (standard Pauli algebra), hence C_ij = (1/2) tr(ρ · 2(n_i·n_j) I) = (n_i·n_j)·tr(ρ) = n_i·n_j — ρ drops out exactly (no Clifford, no Pauli grid needed). ∎

**Consequence:** K3max = max over orbit pairs (n2, n3) of (n₀·n2 + n2·n3 − n₀·n3). Blindness (K3max = 1) ⇔ the Q axis n₀ is **aligned** with a symmetry axis of the orbit; finiteness or the Clifford property alone are NOT sufficient (see application iv).

**Applications (exact values, `spot_checks_gram.py`; reproduces the A-p07 reference reimplementation):**

| # | System | Orbit | Bloch dots (offdiag) | K3max | Reading |
|---|---|---|---|---|---|
| (i) | **k=4 / fusion Z** | 4 tetrahedron directions | {−1/3} | **1 exact** | **blind** (alignment: z = 3-fold axis) |
| (ii) | k=2 (O, z = 4-fold axis) | 6 (±{x,y,z}) | {0, ±1} | 1 exact | blind |
| (iii) | k=8 (I) | 12 icosahedron vertices | {±1, ±1/√5} | 3/√5 ≈ 1.3416 | **sees** (> 1) |
| (iv) | **alignment corollary:** octahedral group, rotated Q=(1,1,0)/√2 | 12 edge midpoints | {0, ±½, ±1} | **3/2 exact** | **not blind** |

**Case list for (i)** (tetrahedron Gram, dot = 1 for equal index, otherwise −1/3): over all 4³ = 64 orbit triples, max(d12 + d23 − d13) = 1 exactly (attained e.g. at n1 = n2, n3 ≠ n1: 1 − 1/3 + 1/3 = 1). **Spot check G1/G2** (below): real k=4 image → 4 tetrahedron directions, all pair dots −1/3, σ₂†Zσ₂ not a ±Pauli; K3max = 1.

**Application (iv) — core of the second repair:** The same abstract group (chiral octahedral group O, 24 elements; as 2O even the **full Clifford image**), but with the Q axis **rotated** by 45° to (1,1,0)/√2, yields K3max = **3/2 exactly** — the Lüders bound, far above 1. **As a slogan, "finite group ⇒ blind" is WRONG.** Blindness is a property of the **aligned pair (group, Q axis)**, not of finiteness and not of the Clifford property per se. **Spot check G4** (Gram) + **G5** (operational quantum cross-check with two real SU(2) qubit unitaries): K3 = 3/2, ρ-independent, == the Gram formula.

**Scope note on (iii):** In the **fixed Z frame** the witness sees K3max = 3/√5 ≈ 1.3416 at k=8. The reference value K3opt = 1.427 = 2cos72° − cos144° comes from the **max-Q protocol** (Q direction optimized ⊥ to the 5-fold axis); in the fixed Z frame no 5-fold axis lies ⊥ z, so it is unreachable there. Not a contradiction but a difference of protocol scope — consistent with application (iv): at k=8 too, K3max depends on the alignment of the Q axis.

---

## Corollaries (across both theorems)

**(a) Dense braid images (contrast).** For k=3 and k≥5 the braid image is dense in the projective unitary group (Freedman–Larsen–Wang universality; finiteness/density classification Jones 1986, Tuba–Wenzl math/9912013 — results by others). The orbit of Q is then dense on the Bloch sphere, C_ij becomes continuous, and K3 can rise to the Lüders bound **3/2** (Budroni–Emary, arXiv:1309.3678). **Spot check 5:** U = exp(−iπX/6) yields K3 = 3/2.

**(b) The driver is the ALIGNMENT, not Clifford/finiteness (sharpened).** The witness is blind exactly when the fixed Q axis is aligned with a symmetry axis of the G orbit whose Gram structure forces max(d12+d23−d13) = 1 (k=4: 3-fold axis → tetrahedron Gram {1,−1/3}; k=2: 4-fold axis → {0,±1}). Finiteness is NOT sufficient (k=8 is finite, yet K3max = 3/√5 > 1); the Clifford property is NOT sufficient (octahedral group with rotated Q: K3max = 3/2, application iv). The Pauli-aligned Clifford case (Theorem A) is the **special case** in which the alignment automatically hits the ±Pauli grid.

**(c) Connection to the magic finding (corrected).** At k=4 the reference measures magic M₂ = 0.5585 (run r1, deposited via `p5a_r1_engine.py`) while K3 = 1.000 (inert). The theorem explains the dissociation structurally: the witness sees only the **orbit Gram structure** of the Heisenberg directions — at k=4 the tetrahedron Gram with dots {1, −1/3} (NOT a "discrete ±Pauli grid {0,±1}"; that holds only at k=2), which caps K3 at 1 — not the magic resource. Blindness despite M₂ = 0.5585 is therefore state-independent alignment structure — **this witness is blind to the k=4 magic**; it says nothing about other witnesses or other Q axes.

## Delimitation (honest): Gottesman–Knill is NOT the argument

Gottesman–Knill gives classical simulability of Clifford circuits on stabilizer states — a statement about computational complexity, not about temporal correlation bounds. Our argument nowhere uses stabilizer simulability; it uses the **orbit Gram geometry** of the Heisenberg directions (Theorem B) and, in the special case, the discreteness of the Pauli orbit (Theorem A). A bridge "simulable ⇒ LGI-inert" would have to be built separately, and is neither needed nor claimed here.

---

## Honesty block (history of the repair, transparent)

- **First version (2026-07-07):** claimed the theorem via "G ⊂ Clifford" alone and gave as an example "the finite k=4 braid image 2T ⊂ Clifford in the Z frame". The corollary phrased blindness via a "discrete ±Pauli grid ⇒ correlators {0,±1}".
- **Two independent refuter passes (2026-07-08)** caught two errors (error class **R6.2** — a confidently plausible but factually wrong intermediate step):
  1. The **real** k=4 fusion-Z image is NOT a Clifford subgroup in the Z frame: σ₁ acts projectively as a 3-fold rotation about z, orbit(Z) = 4 tetrahedron directions with dots −1/3, σ₂†Zσ₂ not a ±Pauli. The correlators are {1, −1/3}, not {0, ±1}. (Refuter B1, `find_counterexample.py`/`ergebnis_b1.json`.)
  2. **Misalignment:** the same abstract 2T/2O group with a rotated Q = (X+Y)/√2 yields K3max = 3/2 — "finite group ⇒ blind" resp. "the driver is Clifford" is wrongly nuanced; the driver is the **alignment**. (Refuter A-p07 reimplementation, `reimpl_ergebnis.json['gram_checks']`.)
- **Repair = this update (2026-07-08):** Theorem A with an explicit Pauli-aligned hypothesis + warning; a new Theorem B (orbit Gram) covers the real k=4 correctly; corollaries (b)/(c) recast onto alignment; spot checks `spot_checks_gram.py` (G1–G5), reproducing the independent A-p07 reimplementation (tet=1, octa=1, ico=3/√5, octa_misaligned=3/2) to machine precision. `spot_checks.py` (Theorem A) remained unchanged and green.

---

## Verification status

| Component | Status |
|---|---|
| Theorem A · Lemma A1 (24/24 elements of 2T abstract; group closure) | verified (CHECK 1) |
| Theorem A · Corollary A-a (48/48 elements of 2O abstract) | verified (CHECK 1b) |
| Theorem A · Lemma A2 (36/36 Pauli pairs; ρ spread 4.44e-16) | verified (CHECK 2) |
| Theorem A · Lemma A3 (36/36 cases, formula == brute force, 56 ρ) | verified (CHECK 3) |
| Theorem A end to end (576 pairs × 14 ρ, max K3 = 1.0 exact) | verified (CHECK 4) |
| Theorem A · contrast non-Clifford (K3 = 3/2) | verified (CHECK 5) |
| **Theorem B · C_ij = n_i·n_j ρ-independent (eigenvalue spread ≤ 6.7e-16)** | **verified (CHECK G1/G5)** |
| **Theorem B (i) k=4 tetrahedron orbit, dots −1/3, K3max = 1 exact** | **verified (CHECK G1/G2), real braid image** |
| **Theorem B (ii) k=2 octahedron aligned, dots {0,±1}, K3max = 1** | **verified (CHECK G2)** |
| **Theorem B (iii) k=8 icosahedron, dots {±1,±1/√5}, K3max = 3/√5** | **verified (CHECK G3), real braid image** |
| **Theorem B (iv) alignment: octahedral group, rotated Q → K3max = 3/2** | **verified (CHECK G4 Gram + G5 operational)** |
| **A-p07 reference reimplementation (tet=1, octa=1, ico=3/√5, octa_mis=3/2) reproduced** | **verified (|diff| = 0 resp. ≤ 1e-15)** |
| Finiteness of 2T/2O/2I as the k=4/2/8 braid image | literature (Jones 1986; Tuba–Wenzl math/9912013), attributed |
| Lüders bound 3/2 | literature (Budroni–Emary 1309.3678), attributed |
| K3opt(k=8, max-Q) = 1.427 · M₂(k=4) = 0.5585 | taken from the reference (r1/r2), not recomputed here |

**Scope limits (honest):** The statement holds for the dichotomic Lüders LGI witness K3 with **fixed** Q and propagators from a group G ⊂ SU(2). Blindness is tied to the **alignment** of Q with an orbit symmetry axis (Theorem B) — not to finiteness or Clifford per se; with a rotated Q the same finite/Clifford G can reach K3 = 3/2. No statement about other measurement schemes (weak/POVM readout), other witnesses, or multi-qubit extensions.

---

**Note on the two raw-output blocks below.** They are the verbatim console log of the original
verification runs (2026-07-07/08) and are reproduced unchanged, in the original German, as a record
of those runs. The scripts were translated afterwards, so their present output is in English; the
glossary below uses the scripts' own wording. `abgeschlossen` = closed (under multiplication) ·
`Spannweite` = spread · `Lueders-Schranke` = Lüders bound · `Abw.` = deviation.

## Raw output Theorem A (spot_checks.py, seed=42, PYTHONIOENCODING=utf-8)

```
== CHECK 1 (Lemma 1): 2T konjugiert Z auf +-Pauli ==
2T: 24 Elemente, alle unitaer mit det=1; abgeschlossen unter Produkt (24x24): True
U^dag Z U in {+-X,+-Y,+-Z} fuer 24/24 Elemente; Verteilung: [('+X', 4), ('+Y', 4), ('+Z', 4), ('-X', 4), ('-Y', 4), ('-Z', 4)]

== CHECK 1b (Korollar a): 2O (k=2-Bild) konjugiert Z auf +-Pauli ==
2O: 48 Elemente, abgeschlossen (48x48): True; U^dag Z U +-Pauli fuer 48/48 Elemente

== CHECK 2 (Lemma 2): alle 36 signierten Pauli-Paare ==
36 Paare: 12x (1/2){A,B} = +-I (gleiche Achse), 24x = 0 (verschiedene Achse)
C = (1/2)tr(rho{A,B}) ueber 206 Zufalls-rho: max. Spannweite = 4.44e-16 (rho-unabhaengig); alle Werte in {-1,0,+1}: True

== CHECK 3 (Lemma 3): Fall-Liste A-E vs. Brute-Force ueber Zufalls-rho ==
  Fall A: 4 Kombinationen, K3-Werte [-3, 1]
  Fall B: 8 Kombinationen, K3-Werte [-1, 1]
  Fall C: 8 Kombinationen, K3-Werte [-1, 1]
  Fall D: 8 Kombinationen, K3-Werte [-1, 1]
  Fall E: 8 Kombinationen, K3-Werte [0]
36/36 Kombinationen: Fall-Formel == Brute-Force (56 rho, Toleranz 1e-10); max K3 = 1

== CHECK 4 (Theorem end-zu-end): alle 576 Propagator-Paare (V1,V2) aus 2T ==
576 Paare, je 14 rho (Ginibre+rein+gemischt): Z2,Z3 stets +-Pauli: True; max K3 = 1.000000000000 <= 1 (erreicht, z.B. V1=V2=I)

== CHECK 5 (Korollar b, Kontrast): Nicht-Clifford-Praezession theta=pi/3 ==
U = exp(-i pi/6 X): U^dag Z U ist +-Pauli? None  -> Lemma 1 verletzt
K3 = 1.500000000000 fuer alle 26 rho (Spannweite 8.88e-16) = 3/2 = Lueders-Schranke (Budroni-Emary, arXiv:1309.3678) > 1

ALLE CHECKS BESTANDEN (seed=42, deterministisch)
```

## Raw output Theorem B (spot_checks_gram.py, deterministic, PYTHONIOENCODING=utf-8, repair 2026-07-08)

```
== CHECK G1 (THEOREM B): Tetraeder-Orbit -- Dots -1/3, K3max = 1 (Blindheit) ==
  4 Einheitsvektoren (1,1,1)-Typ/sqrt3; alle 6 Paar-Dots = -1/3 (max|Abw.| = 1.1e-16); K3max ueber 64 Tripel = 1.000000000000000 = 1 EXAKT (Ref reimpl tet_aligned=1.0000000000000002, |diff|=0.0e+00)
  reales k=4-Bild: proj. Ordnung 12 (T); sigma1 = proj. 3-fach-Drehung um z (Ordnung 3, cos(theta)=-0.500000=-1/2); Orbit(Z) = 4 Tetraeder-Richtungen, Paar-Dots = -1/3 (max|Abw.|=2.2e-16)
  sigma2^dag Z sigma2 -> Bloch(-0.471405,-0.816497,-0.333333) = (-sqrt2/3,-sqrt6/3,-1/3), KEIN +-Pauli; K3max (Q1=Z fest) = 1.000000000000000 = 1

== CHECK G2 (THEOREM B, k=2): Oktaeder aligned z -- Dots {0,+-1}, K3max = 1 ==
  6 Achsen +-{x,y,z}; Paar-Dots (offdiag) = [-1.0, 0.0] in {0,+-1}; K3max (frei) = 1.000000000000000, (Q1=z fest) = 1.000000000000000 = 1 (Ref reimpl octa_aligned=1.0, |diff|=0.0e+00)

== CHECK G3 (THEOREM B, k=8): Ikosaeder -- Dots {+-1,+-1/sqrt5}, K3max = 3/sqrt5 ==
  12 Ecken (0,+-1,+-phi)&zykl.; Paar-Dots (offdiag) = [-1.0, -0.447214, 0.447214] in {+-1,+-1/sqrt5}; K3max = 1.341640786499874 = 3/sqrt5 = 1.341640786499874 > 1 (Ref reimpl ico=1.341640786499874, |diff|=0.0e+00)
  reales k=8-Bild: proj. Ordnung 60 (I); Orbit(Z) = 12 Ikosaeder-Richtungen; K3max (Q1=Z fest) = 1.341640786499873 = 3/sqrt5 -> Zeuge SIEHT k=8

== CHECK G4 (Alignment-Korollar): Oktaeder-Rotationsgruppe, Q=(1,1,0)/sqrt2 -> K3max = 3/2 ==
  24 Vorzeichen-Permutationsmatrizen (det=1) = chirale Oktaedergruppe O; Orbit(Q) = 12 Kantenmitten-Richtungen; Dots (offdiag) = [-1.0, -0.5, -0.0, 0.5]
  K3max (Q1=Q=n0 fest) = 1.500000000000000 = 3/2 EXAKT (Ref reimpl octa_misaligned_110=1.4999999999999996, |diff|=0.0e+00) -> ENDLICHE Gruppe, aber NICHT blind

== CHECK G5 (operational): 2 echte Qubit-Unitaries realisieren das G4-Argmax, K3 = Gram ==
  V1,V2 in SU(2) (unitaer, konstruiert ueber n2,n3); Bloch(Q2)=n2, Bloch(Q3)=n3 verifiziert
  operationaler Lueders-K3 = 1.500000000000000 fuer alle 42 rho (Spannweite 6.7e-16, rho-unabhaengig) == Gram-Formel 1.500000000000000 == 3/2

ALLE GRAM-CHECKS BESTANDEN (deterministisch; Geometrie exakt, Referenz A-p07 reproduziert)
```
