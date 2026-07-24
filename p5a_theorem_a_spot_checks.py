"""spot_checks.py — numerische Spot-Checks fuer die k=4-Blindheits-Note (Auftrag B2).

Deterministisch: np.random.seed(42). Jeder Check ist assert-gesichert; die
Druck-Ausgabe ist der Rohbeleg fuer THEOREM_NOTE.md (dort ans Ende kopiert).

Checks:
  1  Lemma 1  — alle 24 Elemente von 2T (binaer-tetraedrisch) konjugieren Z auf +-Pauli
  1b Korollar a — alle 48 Elemente von 2O (k=2-Bild) ebenso
  2  Lemma 2  — (1/2){A,B} = +-I oder 0 fuer alle 36 signierten Pauli-Paare;
                Lueders-Korrelator dadurch rho-unabhaengig in {-1,0,+1}
  3  Lemma 3  — Fall-Liste A-E symbolisch == Brute-Force ueber Zufalls-rho; max K3 = 1
  4  Theorem  — end-zu-end: alle 576 Propagator-Paare (V1,V2) aus 2T, Zufalls-rho: K3 <= 1
  5  Kontrast — Nicht-Clifford-Propagator (Praezession theta=pi/3) erreicht K3 = 3/2
"""
import numpy as np
from itertools import product, combinations

np.random.seed(42)
TOL = 1e-12

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PNAMES = ["X", "Y", "Z"]
PMATS = [X, Y, Z]


def su2(a, b, c, d):
    """Einheits-Quaternion q = a+bi+cj+dk -> SU(2) via q -> a*I - i(bX + cY + dZ)."""
    return a * I2 - 1j * (b * X + c * Y + d * Z)


def build_2T():
    """Binaer-tetraedrische Gruppe 2T: 8 Lipschitz- + 16 Hurwitz-Einheiten (24 Elemente)."""
    quats = []
    for idx in range(4):
        for s in (1.0, -1.0):
            v = [0.0] * 4
            v[idx] = s
            quats.append(tuple(v))
    for signs in product((0.5, -0.5), repeat=4):
        quats.append(signs)
    return [su2(*q) for q in quats]


def build_2O():
    """Binaer-oktaedrische Gruppe 2O = 2T + 24 Einheiten der Form (+-e_i +- e_j)/sqrt2 (48 Elemente)."""
    els = build_2T()
    r = 1.0 / np.sqrt(2.0)
    for (i, j) in combinations(range(4), 2):
        for si, sj in product((1.0, -1.0), repeat=2):
            v = [0.0] * 4
            v[i] = si * r
            v[j] = sj * r
            els.append(su2(*v))
    return els


def signed_pauli_of(M):
    """(sign, name) falls M == sign*P fuer P in {X,Y,Z}, sonst None."""
    for name, P in zip(PNAMES, PMATS):
        for s in (1, -1):
            if np.max(np.abs(M - s * P)) < 1e-9:
                return s, name
    return None


def in_set(M, els):
    return any(np.max(np.abs(M - E)) < 1e-9 for E in els)


def rand_rho(n):
    """n Ginibre-Zufalls-Dichtematrizen + 5 reine Zustaende + maximal gemischt."""
    out = []
    for _ in range(n):
        G = np.random.randn(2, 2) + 1j * np.random.randn(2, 2)
        r = G @ G.conj().T
        out.append(r / np.trace(r).real)
    for v in ([1, 0], [0, 1], [1, 1], [1, 1j], [1, 0.3 + 0.4j]):
        v = np.array(v, dtype=complex)
        v = v / np.linalg.norm(v)
        out.append(np.outer(v, v.conj()))
    out.append(I2 / 2)
    return out


def corr(rho, A, B):
    """Lueders-Korrelator fuer dichotome Observablen: C = (1/2) tr(rho {A,B})."""
    return (0.5 * np.trace(rho @ (A @ B + B @ A))).real


# ---------- CHECK 1: Lemma 1 auf 2T ----------
print("== CHECK 1 (Lemma 1): 2T konjugiert Z auf +-Pauli ==")
T24 = build_2T()
assert len(T24) == 24
for U in T24:
    assert np.max(np.abs(U.conj().T @ U - I2)) < TOL, "nicht unitaer"
    assert abs(np.linalg.det(U) - 1) < 1e-9, "det != 1"
closure = all(in_set(U @ V, T24) for U in T24 for V in T24)
assert closure
print(f"2T: 24 Elemente, alle unitaer mit det=1; abgeschlossen unter Produkt (24x24): {closure}")
dist = {}
for U in T24:
    r = signed_pauli_of(U.conj().T @ Z @ U)
    assert r is not None, "U^dag Z U ist kein +-Pauli"
    key = ("+" if r[0] > 0 else "-") + r[1]
    dist[key] = dist.get(key, 0) + 1
print(f"U^dag Z U in {{+-X,+-Y,+-Z}} fuer 24/24 Elemente; Verteilung: {sorted(dist.items())}")

# ---------- CHECK 1b: Korollar (a), 2O (k=2-Bild) ----------
print("\n== CHECK 1b (Korollar a): 2O (k=2-Bild) konjugiert Z auf +-Pauli ==")
O48 = build_2O()
assert len(O48) == 48
closure_O = all(in_set(U @ V, O48) for U in O48 for V in O48)
assert closure_O
n_ok = sum(1 for U in O48 if signed_pauli_of(U.conj().T @ Z @ U) is not None)
assert n_ok == 48
print(f"2O: 48 Elemente, abgeschlossen (48x48): {closure_O}; U^dag Z U +-Pauli fuer {n_ok}/48 Elemente")

# ---------- CHECK 2: Lemma 2 ----------
print("\n== CHECK 2 (Lemma 2): alle 36 signierten Pauli-Paare ==")
signed = [(s, n, s * P) for n, P in zip(PNAMES, PMATS) for s in (1, -1)]
n_id, n_zero = 0, 0
for sA, nA, A in signed:
    for sB, nB, B in signed:
        M = 0.5 * (A @ B + B @ A)
        if nA == nB:  # B = +-A
            assert np.max(np.abs(M - sA * sB * I2)) < TOL
            n_id += 1
        else:  # verschiedene Pauli-Achsen antikommutieren
            assert np.max(np.abs(M)) < TOL
            n_zero += 1
print(f"36 Paare: {n_id}x (1/2){{A,B}} = +-I (gleiche Achse), {n_zero}x = 0 (verschiedene Achse)")
rhos = rand_rho(200)
spread_max, allowed = 0.0, True
for _, _, A in signed:
    for _, _, B in signed:
        vals = np.array([corr(r, A, B) for r in rhos])
        spread_max = max(spread_max, float(np.ptp(vals)))
        allowed &= min(abs(vals[0] - t) for t in (-1.0, 0.0, 1.0)) < 1e-10
assert spread_max < 1e-10 and allowed
print(f"C = (1/2)tr(rho{{A,B}}) ueber {len(rhos)} Zufalls-rho: max. Spannweite = {spread_max:.2e} "
      f"(rho-unabhaengig); alle Werte in {{-1,0,+1}}: {allowed}")

# ---------- CHECK 3: Lemma 3, Fall-Liste vs. Brute-Force ----------
print("\n== CHECK 3 (Lemma 3): Fall-Liste A-E vs. Brute-Force ueber Zufalls-rho ==")
rhos = rand_rho(50)
cases, max_K3 = {}, -np.inf
for n2, P2 in zip(PNAMES, PMATS):
    for s2 in (1, -1):
        for n3, P3 in zip(PNAMES, PMATS):
            for s3 in (1, -1):
                Z2, Z3 = s2 * P2, s3 * P3
                C21f = s2 if n2 == "Z" else 0
                C32f = s2 * s3 if n3 == n2 else 0
                C31f = s3 if n3 == "Z" else 0
                K3f = C21f + C32f - C31f
                for r in rhos:
                    K3n = corr(r, Z2, Z) + corr(r, Z3, Z2) - corr(r, Z3, Z)
                    assert abs(K3n - K3f) < 1e-10
                if n2 == "Z" and n3 == "Z":
                    case = "A"
                elif n2 == "Z":
                    case = "B"
                elif n3 == n2:
                    case = "C"
                elif n3 == "Z":
                    case = "D"
                else:
                    case = "E"
                cases.setdefault(case, []).append(K3f)
                max_K3 = max(max_K3, K3f)
for c in "ABCDE":
    v = cases[c]
    print(f"  Fall {c}: {len(v)} Kombinationen, K3-Werte {sorted(set(v))}")
assert sum(len(v) for v in cases.values()) == 36 and max_K3 <= 1 + TOL
print(f"36/36 Kombinationen: Fall-Formel == Brute-Force ({len(rhos)} rho, Toleranz 1e-10); max K3 = {max_K3}")

# ---------- CHECK 4: Theorem end-zu-end ueber 2T x 2T ----------
print("\n== CHECK 4 (Theorem end-zu-end): alle 576 Propagator-Paare (V1,V2) aus 2T ==")
rhos = rand_rho(8)
max_K3, n_pairs, all_pauli = -np.inf, 0, True
for V1 in T24:
    Z2 = V1.conj().T @ Z @ V1
    for V2 in T24:
        W = V2 @ V1
        Z3 = W.conj().T @ Z @ W
        all_pauli &= (signed_pauli_of(Z2) is not None) and (signed_pauli_of(Z3) is not None)
        n_pairs += 1
        for r in rhos:
            K3 = corr(r, Z2, Z) + corr(r, Z3, Z2) - corr(r, Z3, Z)
            max_K3 = max(max_K3, K3)
assert all_pauli and n_pairs == 576 and max_K3 <= 1 + 1e-10
print(f"{n_pairs} Paare, je {len(rhos)} rho (Ginibre+rein+gemischt): Z2,Z3 stets +-Pauli: {all_pauli}; "
      f"max K3 = {max_K3:.12f} <= 1 (erreicht, z.B. V1=V2=I)")

# ---------- CHECK 5: Kontrast — Nicht-Clifford-Propagator ----------
print("\n== CHECK 5 (Korollar b, Kontrast): Nicht-Clifford-Praezession theta=pi/3 ==")
th = np.pi / 3
U = np.cos(th / 2) * I2 - 1j * np.sin(th / 2) * X  # exp(-i theta X / 2), kein Clifford-Element
r1 = signed_pauli_of(U.conj().T @ Z @ U)
Z2 = U.conj().T @ Z @ U
W = U @ U
Z3 = W.conj().T @ Z @ W
K3_vals = np.array([corr(r, Z2, Z) + corr(r, Z3, Z2) - corr(r, Z3, Z) for r in rand_rho(20)])
assert r1 is None
assert np.max(np.abs(K3_vals - 1.5)) < 1e-10
print(f"U = exp(-i pi/6 X): U^dag Z U ist +-Pauli? {r1}  -> Lemma 1 verletzt")
print(f"K3 = {K3_vals[0]:.12f} fuer alle {len(K3_vals)} rho (Spannweite {np.ptp(K3_vals):.2e}) "
      f"= 3/2 = Lueders-Schranke (Budroni-Emary, arXiv:1309.3678) > 1")

print("\nALLE CHECKS BESTANDEN (seed=42, deterministisch)")
