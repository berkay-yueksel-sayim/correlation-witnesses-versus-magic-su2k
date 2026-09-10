"""spot_checks.py — numerical spot checks for the k=4 blindness note (task B2).

Deterministic: np.random.seed(42). Every check is assert-guarded; the
printed output is the raw record for THEOREM_NOTE.md (copied to its end there).

Checks:
  1  Lemma 1  — all 24 elements of 2T (binary tetrahedral) conjugate Z to +-Pauli
  1b Corollary a — all 48 elements of 2O (k=2 image) likewise
  2  Lemma 2  — (1/2){A,B} = +-I or 0 for all 36 signed Pauli pairs;
                the Lueders correlator is thereby rho-independent in {-1,0,+1}
  3  Lemma 3  — case list A-E symbolic == brute force over random rho; max K3 = 1
  4  Theorem  — end to end: all 576 propagator pairs (V1,V2) from 2T, random rho: K3 <= 1
  5  Contrast — non-Clifford propagator (precession theta=pi/3) reaches K3 = 3/2
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
    """Binary tetrahedral group 2T: 8 Lipschitz + 16 Hurwitz units (24 elements)."""
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
    """Binary octahedral group 2O = 2T + 24 units of the form (+-e_i +- e_j)/sqrt2 (48 elements)."""
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
    """(sign, name) if M == sign*P for P in {X,Y,Z}, else None."""
    for name, P in zip(PNAMES, PMATS):
        for s in (1, -1):
            if np.max(np.abs(M - s * P)) < 1e-9:
                return s, name
    return None


def in_set(M, els):
    return any(np.max(np.abs(M - E)) < 1e-9 for E in els)


def rand_rho(n):
    """n Ginibre random density matrices + 5 pure states + maximally mixed."""
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
    """Lueders correlator for dichotomic observables: C = (1/2) tr(rho {A,B})."""
    return (0.5 * np.trace(rho @ (A @ B + B @ A))).real


# ---------- CHECK 1: Lemma 1 on 2T ----------
print("== CHECK 1 (Lemma 1): 2T conjugates Z to +-Pauli ==")
T24 = build_2T()
assert len(T24) == 24
for U in T24:
    assert np.max(np.abs(U.conj().T @ U - I2)) < TOL, "not unitary"
    assert abs(np.linalg.det(U) - 1) < 1e-9, "det != 1"
closure = all(in_set(U @ V, T24) for U in T24 for V in T24)
assert closure
print(f"2T: 24 elements, all unitary with det=1; closed under multiplication (24x24): {closure}")
dist = {}
for U in T24:
    r = signed_pauli_of(U.conj().T @ Z @ U)
    assert r is not None, "U^dag Z U is not a +-Pauli"
    key = ("+" if r[0] > 0 else "-") + r[1]
    dist[key] = dist.get(key, 0) + 1
print(f"U^dag Z U in {{+-X,+-Y,+-Z}} for 24/24 elements; distribution: {sorted(dist.items())}")

# ---------- CHECK 1b: corollary (a), 2O (k=2 image) ----------
print("\n== CHECK 1b (corollary a): 2O (k=2 image) conjugates Z to +-Pauli ==")
O48 = build_2O()
assert len(O48) == 48
closure_O = all(in_set(U @ V, O48) for U in O48 for V in O48)
assert closure_O
n_ok = sum(1 for U in O48 if signed_pauli_of(U.conj().T @ Z @ U) is not None)
assert n_ok == 48
print(f"2O: 48 elements, closed (48x48): {closure_O}; U^dag Z U +-Pauli for {n_ok}/48 elements")

# ---------- CHECK 2: Lemma 2 ----------
print("\n== CHECK 2 (Lemma 2): all 36 signed Pauli pairs ==")
signed = [(s, n, s * P) for n, P in zip(PNAMES, PMATS) for s in (1, -1)]
n_id, n_zero = 0, 0
for sA, nA, A in signed:
    for sB, nB, B in signed:
        M = 0.5 * (A @ B + B @ A)
        if nA == nB:  # B = +-A
            assert np.max(np.abs(M - sA * sB * I2)) < TOL
            n_id += 1
        else:  # different Pauli axes anticommute
            assert np.max(np.abs(M)) < TOL
            n_zero += 1
print(f"36 pairs: {n_id}x (1/2){{A,B}} = +-I (same axis), {n_zero}x = 0 (different axis)")
rhos = rand_rho(200)
spread_max, allowed = 0.0, True
for _, _, A in signed:
    for _, _, B in signed:
        vals = np.array([corr(r, A, B) for r in rhos])
        spread_max = max(spread_max, float(np.ptp(vals)))
        allowed &= min(abs(vals[0] - t) for t in (-1.0, 0.0, 1.0)) < 1e-10
assert spread_max < 1e-10 and allowed
print(f"C = (1/2)tr(rho{{A,B}}) over {len(rhos)} random rho: max spread = {spread_max:.2e} "
      f"(rho-independent); all values in {{-1,0,+1}}: {allowed}")

# ---------- CHECK 3: Lemma 3, case list vs. brute force ----------
print("\n== CHECK 3 (Lemma 3): case list A-E vs. brute force over random rho ==")
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
    print(f"  Case {c}: {len(v)} combinations, K3 values {sorted(set(v))}")
assert sum(len(v) for v in cases.values()) == 36 and max_K3 <= 1 + TOL
print(f"36/36 combinations: case formula == brute force ({len(rhos)} rho, tolerance 1e-10); max K3 = {max_K3}")

# ---------- CHECK 4: theorem end to end over 2T x 2T ----------
print("\n== CHECK 4 (theorem end to end): all 576 propagator pairs (V1,V2) from 2T ==")
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
print(f"{n_pairs} pairs, each {len(rhos)} rho (Ginibre+pure+mixed): Z2,Z3 always +-Pauli: {all_pauli}; "
      f"max K3 = {max_K3:.12f} <= 1 (attained, e.g. V1=V2=I)")

# ---------- CHECK 5: contrast — non-Clifford propagator ----------
print("\n== CHECK 5 (corollary b, contrast): non-Clifford precession theta=pi/3 ==")
th = np.pi / 3
U = np.cos(th / 2) * I2 - 1j * np.sin(th / 2) * X  # exp(-i theta X / 2), not a Clifford element
r1 = signed_pauli_of(U.conj().T @ Z @ U)
Z2 = U.conj().T @ Z @ U
W = U @ U
Z3 = W.conj().T @ Z @ W
K3_vals = np.array([corr(r, Z2, Z) + corr(r, Z3, Z2) - corr(r, Z3, Z) for r in rand_rho(20)])
assert r1 is None
assert np.max(np.abs(K3_vals - 1.5)) < 1e-10
print(f"U = exp(-i pi/6 X): U^dag Z U is +-Pauli? {r1}  -> Lemma 1 violated")
print(f"K3 = {K3_vals[0]:.12f} for all {len(K3_vals)} rho (spread {np.ptp(K3_vals):.2e}) "
      f"= 3/2 = Lueders bound (Budroni-Emary, arXiv:1309.3678) > 1")

print("\nALL CHECKS PASSED (seed=42, deterministic)")
