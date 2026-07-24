# -*- coding: utf-8 -*-
"""spot_checks_gram.py -- Spot-Checks fuer THEOREM B (Orbit-Gram-Geometrie) + Alignment-Korollar.

Reparatur-Lauf 2026-07-08. ERGAENZT spot_checks.py (Theorem A) -- ersetzt/aendert es NICHT.
Deterministisch (seed=7 nur fuer rho-Sampling im operationalen Check; alle Geometrie exakt,
kein Sampling). PYTHONIOENCODING=utf-8.

Gram-Kern: fuer Single-Qubit-Dichotome ist {n.sigma, m.sigma} = 2(n.m) I, also der Lueders-
Korrelator C_ij = (1/2)tr(rho{Qi,Qj}) = n_i.n_j fuer JEDES rho (rho-unabhaengig, exakt). K3max
folgt allein aus der Orbit-Gram-Struktur der Heisenberg-Richtungen. Blindheit <=> Alignment der
Q-Achse mit einer Symmetrieachse des Orbits, NICHT Endlichkeit / Clifford per se.

Checks (Nummerierung wie in THEOREM_NOTE.md, THEOREM B):
  G1  Tetraeder-Orbit-Geometrie: 4 Einheitsvektoren (1,1,1)-Typ/sqrt3, alle Paar-Dots = -1/3,
      K3max ueber alle Tripel = 1 EXAKT (Blindheit). PLUS: reales k=4-Braid-Bild im Fusions-Z-
      Rahmen liefert genau diesen Tetraeder-Orbit (proj. Ordnung 12 = T; sigma2^dag Z sigma2
      = Bloch(-sqrt2/3,-sqrt6/3,-1/3), KEIN +-Pauli).
  G2  Oktaeder aligned z (= k=2-Orbit +-{X,Y,Z}): Dots {0,+-1}, K3max = 1.
  G3  Ikosaeder (12 Ecken, goldener Schnitt): Dots {+-1,+-1/sqrt5}, K3max = 3/sqrt5 (> 1).
      PLUS: reales k=8-Braid-Bild (proj. Ordnung 60 = I) liefert genau diesen Orbit.
  G4  Alignment-Korollar: Oktaeder-ROTATIONSGRUPPE (24 Vorzeichen-Permutationsmatrizen, det=1)
      auf gedrehtem Q=(1,1,0)/sqrt2: Orbit 12, K3max = 3/2 EXAKT -> NICHT blind.
  G5  operationale Quanten-Gegenprobe zu G4: baue 2 Qubit-Unitaries (SU(2)-Rotationen) zu zwei
      Orbit-Richtungen, K3 aus echten Lueders-Korrelatoren == Gram-Formel == 3/2.

Referenz-Gegenrechnung (unabhaengige Reimpl A-p07, reimpl_ergebnis.json['gram_checks']):
  tet_aligned=1.0, octa_aligned=1.0, ico=1.341640786499874, octa_misaligned_110=1.5.
"""
import os
import sys
from itertools import combinations, permutations, product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from p5a_theorem_find_counterexample import braid_data  # NUR die Gate-Konvention; Rest eigenstaendig

TOL = 1e-12
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIG = [X, Y, Z]
ez = np.array([0.0, 0.0, 1.0])

# Referenzwerte der unabhaengigen A-p07-Reimpl (reimpl_ergebnis.json['gram_checks'])
REF = {"tet": 1.0000000000000002, "octa": 1.0, "ico": 1.341640786499874, "octa_mis": 1.4999999999999996}


# ---------- eigenstaendige Hilfen ----------
def rot(U):
    """SO(3)-Bild: R[i,j] = (1/2) Re tr(sigma_i U sigma_j U^dag). U (m.sigma) U^dag = (R m).sigma."""
    Ud = U.conj().T
    return np.array([[0.5 * np.real(np.trace(SIG[i] @ U @ SIG[j] @ Ud))
                      for j in range(3)] for i in range(3)])


def canon(U):
    return tuple(np.round(rot(U).flatten(), 6))


def close_group(gens, cap=400):
    """BFS-Gruppenschluss projektiv (mod Phase) ueber das SO(3)-Bild."""
    elems = {canon(I2): I2.copy()}
    frontier = [I2.copy()]
    gg = gens + [g.conj().T for g in gens]
    while frontier:
        U = frontier.pop()
        for g in gg:
            V = g @ U
            c = canon(V)
            if c not in elems:
                assert len(elems) < cap, "Gruppe schliesst nicht (cap erreicht)"
                elems[c] = V
                frontier.append(V)
    return list(elems.values())


def proj_order(U, nmax=100):
    P = U.copy()
    c0 = canon(I2)
    for n in range(1, nmax + 1):
        if canon(P) == c0:
            return n
        P = U @ P
    return None


def bloch(Q):
    return np.array([0.5 * np.real(np.trace(SIG[i] @ Q)) for i in range(3)])


def sigvec(n):
    """n.sigma."""
    return n[0] * X + n[1] * Y + n[2] * Z


def orbit_bloch(elems, Q0):
    """Orbit-Richtungen n(g) mit g^dag Q0 g = n.sigma, dedupliziert."""
    vecs = []
    for g in elems:
        n = bloch(g.conj().T @ Q0 @ g)
        if not any(np.linalg.norm(n - v) < 1e-8 for v in vecs):
            vecs.append(n)
    return vecs


def k3max_orbit(orb, n1_fixed=None):
    """max (n1.n2 + n2.n3 - n1.n3). n1 fest (=Q) falls gegeben, sonst frei ueber Orbit."""
    firsts = [n1_fixed] if n1_fixed is not None else orb
    best, arg = -np.inf, None
    for a in firsts:
        for b in orb:
            for c in orb:
                v = float(np.dot(a, b) + np.dot(b, c) - np.dot(a, c))
                if v > best:
                    best, arg = v, (a, b, c)
    return best, arg


def is_pm_pauli(n):
    axes = [s * e for e in np.eye(3) for s in (1.0, -1.0)]
    return min(np.linalg.norm(n - a) for a in axes) < 1e-9


def su2_rot_map(a, b):
    """SU(2)-Element U, dessen SO(3)-Rotation den Einheitsvektor a auf b abbildet (R_U a = b)."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    c = float(np.dot(a, b))
    if c > 1 - 1e-12:
        return I2.copy()
    if c < -1 + 1e-12:                       # Drehung um pi um eine Achse senkrecht zu a
        perp = np.cross(a, np.array([1.0, 0.0, 0.0]))
        if np.linalg.norm(perp) < 1e-6:
            perp = np.cross(a, np.array([0.0, 1.0, 0.0]))
        k = perp / np.linalg.norm(perp)
        return -1j * sigvec(k)
    k = np.cross(a, b)
    k = k / np.linalg.norm(k)
    th = np.arccos(c)
    return np.cos(th / 2) * I2 - 1j * np.sin(th / 2) * sigvec(k)


# =====================================================================================
# CHECK G1: Tetraeder-Orbit-Geometrie + reales k=4-Braid-Bild
# =====================================================================================
print("== CHECK G1 (THEOREM B): Tetraeder-Orbit -- Dots -1/3, K3max = 1 (Blindheit) ==")
tet = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float) / np.sqrt(3.0)
assert all(abs(np.linalg.norm(v) - 1.0) < TOL for v in tet)
worst = max(abs(float(np.dot(a, b)) + 1.0 / 3.0) for a, b in combinations(tet, 2))
assert worst < TOL, f"Tetraeder-Dot weicht von -1/3 ab: {worst:.2e}"
k3_tet, _ = k3max_orbit(list(tet))
assert abs(k3_tet - 1.0) < TOL, f"K3max Tetraeder = {k3_tet}"
print(f"  4 Einheitsvektoren (1,1,1)-Typ/sqrt3; alle 6 Paar-Dots = -1/3 (max|Abw.| = {worst:.1e}); "
      f"K3max ueber 64 Tripel = {k3_tet:.15f} = 1 EXAKT (Ref reimpl tet_aligned={REF['tet']}, "
      f"|diff|={abs(k3_tet - REF['tet']):.1e})")

# reales k=4-Braid-Bild im Fusions-Z-Rahmen
s1, s2, _F = braid_data(4)
assert np.max(np.abs(s1 @ s2 @ s1 - s2 @ s1 @ s2)) < TOL, "YBE verletzt"
Gk4 = close_group([s1, s2])
assert len(Gk4) == 12, f"proj. Ordnung {len(Gk4)} != 12"
o1 = proj_order(s1)
R1 = rot(s1)
cos_th = 0.5 * (np.trace(R1) - 1.0)
assert o1 == 3 and np.linalg.norm(R1 @ ez - ez) < TOL and abs(cos_th + 0.5) < TOL
orb4 = orbit_bloch(Gk4, Z)
assert len(orb4) == 4 and any(np.linalg.norm(v - ez) < 1e-9 for v in orb4)
worst4 = max(abs(float(np.dot(a, b)) + 1.0 / 3.0) for a, b in combinations(orb4, 2))
assert worst4 < TOL
n_s2 = bloch(s2.conj().T @ Z @ s2)
assert not is_pm_pauli(n_s2) and abs(n_s2[2] + 1.0 / 3.0) < 1e-9
k3_k4, _ = k3max_orbit(orb4, n1_fixed=ez)
assert abs(k3_k4 - 1.0) < TOL
print(f"  reales k=4-Bild: proj. Ordnung {len(Gk4)} (T); sigma1 = proj. 3-fach-Drehung um z "
      f"(Ordnung {o1}, cos(theta)={cos_th:+.6f}=-1/2); Orbit(Z) = 4 Tetraeder-Richtungen, "
      f"Paar-Dots = -1/3 (max|Abw.|={worst4:.1e})")
print(f"  sigma2^dag Z sigma2 -> Bloch({n_s2[0]:+.6f},{n_s2[1]:+.6f},{n_s2[2]:+.6f}) "
      f"= (-sqrt2/3,-sqrt6/3,-1/3), KEIN +-Pauli; K3max (Q1=Z fest) = {k3_k4:.15f} = 1")

# =====================================================================================
# CHECK G2: Oktaeder aligned z (k=2-Orbit)
# =====================================================================================
print("\n== CHECK G2 (THEOREM B, k=2): Oktaeder aligned z -- Dots {0,+-1}, K3max = 1 ==")
octa = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float)
dots_o = sorted({round(float(np.dot(a, b)), 9) for a, b in combinations(octa, 2)})
assert set(dots_o) <= {-1.0, 0.0, 1.0}
k3_octa_free, _ = k3max_orbit(list(octa))
k3_octa_fix, _ = k3max_orbit(list(octa), n1_fixed=ez)
assert abs(k3_octa_free - 1.0) < TOL and abs(k3_octa_fix - 1.0) < TOL
print(f"  6 Achsen +-{{x,y,z}}; Paar-Dots (offdiag) = {dots_o} in {{0,+-1}}; "
      f"K3max (frei) = {k3_octa_free:.15f}, (Q1=z fest) = {k3_octa_fix:.15f} = 1 "
      f"(Ref reimpl octa_aligned={REF['octa']}, |diff|={abs(k3_octa_free - REF['octa']):.1e})")

# =====================================================================================
# CHECK G3: Ikosaeder (goldener Schnitt) + reales k=8-Braid-Bild
# =====================================================================================
print("\n== CHECK G3 (THEOREM B, k=8): Ikosaeder -- Dots {+-1,+-1/sqrt5}, K3max = 3/sqrt5 ==")
phi = (1.0 + np.sqrt(5.0)) / 2.0
raw = []
for s1_, s2_ in product((1.0, -1.0), repeat=2):
    raw.append((0.0, s1_ * 1.0, s2_ * phi))
    raw.append((s1_ * 1.0, s2_ * phi, 0.0))
    raw.append((s1_ * phi, 0.0, s2_ * 1.0))
ico = np.array(raw) / np.sqrt(1.0 + phi * phi)
assert len(ico) == 12 and all(abs(np.linalg.norm(v) - 1.0) < TOL for v in ico)
allowed = [-1.0, -1.0 / np.sqrt(5.0), 1.0 / np.sqrt(5.0), 1.0]
dots_i = sorted({round(float(np.dot(a, b)), 9) for a, b in combinations(ico, 2)})
assert all(min(abs(d - t) for t in allowed) < 1e-9 for d in dots_i)
k3_ico, _ = k3max_orbit(list(ico))
t35 = 3.0 / np.sqrt(5.0)
assert abs(k3_ico - t35) < TOL, f"K3max Ikosaeder = {k3_ico}"
print(f"  12 Ecken (0,+-1,+-phi)&zykl.; Paar-Dots (offdiag) = {[round(d, 6) for d in dots_i]} "
      f"in {{+-1,+-1/sqrt5}}; K3max = {k3_ico:.15f} = 3/sqrt5 = {t35:.15f} > 1 "
      f"(Ref reimpl ico={REF['ico']}, |diff|={abs(k3_ico - REF['ico']):.1e})")
# reales k=8-Braid-Bild
t1, t2, _F8 = braid_data(8)
Gk8 = close_group([t1, t2])
assert len(Gk8) == 60
orb8 = orbit_bloch(Gk8, Z)
assert len(orb8) == 12
k3_k8, _ = k3max_orbit(orb8, n1_fixed=ez)
assert abs(k3_k8 - t35) < TOL
print(f"  reales k=8-Bild: proj. Ordnung {len(Gk8)} (I); Orbit(Z) = 12 Ikosaeder-Richtungen; "
      f"K3max (Q1=Z fest) = {k3_k8:.15f} = 3/sqrt5 -> Zeuge SIEHT k=8")

# =====================================================================================
# CHECK G4: Alignment-Korollar -- Oktaeder-ROTATIONSGRUPPE, gedrehtes Q
# =====================================================================================
print("\n== CHECK G4 (Alignment-Korollar): Oktaeder-Rotationsgruppe, Q=(1,1,0)/sqrt2 -> K3max = 3/2 ==")
rots = []
for perm in permutations(range(3)):
    P = np.zeros((3, 3))
    for i, p in enumerate(perm):
        P[i, p] = 1.0
    for signs in product((1.0, -1.0), repeat=3):
        M = P * np.array(signs)[:, None]
        if abs(np.linalg.det(M) - 1.0) < 1e-9:
            rots.append(M)
assert len(rots) == 24, f"Rotationsgruppe hat {len(rots)} Elemente != 24"
n0 = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
orbO = []
for R in rots:
    v = R @ n0
    if not any(np.linalg.norm(v - w) < 1e-9 for w in orbO):
        orbO.append(v)
assert len(orbO) == 12, f"Orbit-Groesse {len(orbO)} != 12"
dots_O = sorted({round(float(np.dot(a, b)), 9) for a, b in combinations(orbO, 2)})
k3_O, arg_O = k3max_orbit(orbO, n1_fixed=n0)
assert abs(k3_O - 1.5) < TOL, f"K3max = {k3_O}"
print(f"  24 Vorzeichen-Permutationsmatrizen (det=1) = chirale Oktaedergruppe O; "
      f"Orbit(Q) = {len(orbO)} Kantenmitten-Richtungen; Dots (offdiag) = {[round(d, 4) for d in dots_O]}")
print(f"  K3max (Q1=Q=n0 fest) = {k3_O:.15f} = 3/2 EXAKT (Ref reimpl octa_misaligned_110={REF['octa_mis']}, "
      f"|diff|={abs(k3_O - REF['octa_mis']):.1e}) -> ENDLICHE Gruppe, aber NICHT blind")

# =====================================================================================
# CHECK G5: operationale Quanten-Gegenprobe zu G4
# =====================================================================================
print("\n== CHECK G5 (operational): 2 echte Qubit-Unitaries realisieren das G4-Argmax, K3 = Gram ==")
rng = np.random.default_rng(7)
_, (nA, n2, n3) = k3max_orbit(orbO, n1_fixed=n0)   # nA == n0
Q1 = sigvec(n0)
V1 = su2_rot_map(n2, n0)                            # R_{V1} bildet n2 -> n0  => V1^dag Q1 V1 hat Bloch n2
W = su2_rot_map(n3, n0)                             # Gesamt-Propagator V2 V1 realisiert n3
V2 = W @ V1.conj().T
Q2 = V1.conj().T @ Q1 @ V1
Q3 = W.conj().T @ Q1 @ W
# Konstruktion verifizieren: Heisenberg-Richtungen == Orbit-Ziele
assert np.linalg.norm(bloch(Q2) - n2) < 1e-10 and np.linalg.norm(bloch(Q3) - n3) < 1e-10
# V1, V2 unitaer
for U in (V1, V2, W):
    assert np.max(np.abs(U.conj().T @ U - I2)) < 1e-10


def rand_rho():
    G = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
    r = G @ G.conj().T
    return r / np.trace(r).real


def luders(rho, A, B):
    return float((0.5 * np.trace(rho @ (A @ B + B @ A))).real)


k3_vals = []
for _ in range(40):
    rho = rand_rho()
    C21 = luders(rho, Q2, Q1)
    C32 = luders(rho, Q3, Q2)
    C31 = luders(rho, Q3, Q1)
    k3_vals.append(C21 + C32 - C31)
for extra in (np.array([[1, 0], [0, 0]], dtype=complex), I2 / 2):
    C21 = luders(extra, Q2, Q1)
    C32 = luders(extra, Q3, Q2)
    C31 = luders(extra, Q3, Q1)
    k3_vals.append(C21 + C32 - C31)
k3_vals = np.array(k3_vals)
gram = float(np.dot(n0, n2) + np.dot(n2, n3) - np.dot(n0, n3))
spread = float(np.ptp(k3_vals))
assert spread < 1e-12, f"K3 rho-abhaengig: Spannweite {spread}"
assert abs(k3_vals[0] - gram) < 1e-12 and abs(gram - 1.5) < TOL
print(f"  V1,V2 in SU(2) (unitaer, konstruiert ueber n2,n3); Bloch(Q2)=n2, Bloch(Q3)=n3 verifiziert")
print(f"  operationaler Lueders-K3 = {k3_vals[0]:.15f} fuer alle {len(k3_vals)} rho "
      f"(Spannweite {spread:.1e}, rho-unabhaengig) == Gram-Formel {gram:.15f} == 3/2")

print("\nALLE GRAM-CHECKS BESTANDEN (deterministisch; Geometrie exakt, Referenz A-p07 reproduziert)")
