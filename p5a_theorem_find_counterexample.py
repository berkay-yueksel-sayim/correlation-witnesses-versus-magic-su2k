"""
B1 - Numerische Gegenbeispiel-Suche fuer das k=4-Blindheits-Theorem (2026-07-07).
Claim: SU(2)_k=4 Single-Qubit-Braiding (endliches Braid-Bild, erwartet 2T, projektiv T)
=> dichotomer Lueders-LGI-Zeuge K3 = C21 + C32 - C31 im festen Pauli-Z-Rahmen erreicht
NIE mehr als 1, fuer ALLE Braid-Woerter, ALLE rho, alle drei Zeitpunkte.
Methode: exakte Enumeration (Gruppe -> Orbit von Z -> alle Paare/Tripel), KEIN Sampling.
Konventionen fuer Braid-Gates identisch zur Referenz-Brick-Gate-Engine dieses Projekts
(unabhaengig re-implementiert, gleiche Formeln). Deterministisch, kein Zufall noetig.
C_ij = (1/2) tr(rho {Q_i, Q_j}), Q_i = Heisenberg-entwickeltes Q=Z (Budroni-Emary-Standard).
"""
import json
import numpy as np

TOL = 1e-10
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIG = [X, Y, Z]

# ---- SU(2)_k Braid-Gates, j=1/2 (Konvention wie engine.py r1) ----
def braid_data(k):
    d = 2 * np.cos(np.pi / (k + 2))
    h_half = 0.5 * 1.5 / (k + 2)
    h1 = 2.0 / (k + 2)
    R0 = -np.exp(1j * np.pi * (0.0 - 2 * h_half))
    R1 = +np.exp(1j * np.pi * (h1 - 2 * h_half))
    s1 = np.diag([R0, R1])
    root = np.sqrt(max(d * d - 1.0, 0.0)) / d
    F = np.array([[1.0 / d, root], [root, -1.0 / d]], dtype=complex)
    s2 = F @ s1 @ F
    return s1, s2, F

# ---- projektiv (mod Phase) via SO(3)-Bild ----
def to_SO3(U):
    R = np.zeros((3, 3))
    Ud = U.conj().T
    for i in range(3):
        for j in range(3):
            R[i, j] = 0.5 * np.real(np.trace(SIG[i] @ U @ SIG[j] @ Ud))
    return R

def canon(U, dec=6):
    return tuple(np.round(to_SO3(U).flatten(), dec))

def enumerate_group(gens, cap=5000):
    seen = {canon(I2): I2.copy()}
    frontier = [I2.copy()]
    allg = gens + [g.conj().T for g in gens]
    while frontier:
        U = frontier.pop()
        for g in allg:
            V = g @ U
            c = canon(V)
            if c not in seen:
                seen[c] = V
                frontier.append(V)
                if len(seen) > cap:
                    return list(seen.values()), False
    return list(seen.values()), True

def closure_check(elems, gens):
    """Abschluss mit harter Toleranz TOL: fuer jedes U, g liegt gU in der Menge (SO3-Frobenius)."""
    Rs = np.array([to_SO3(U) for U in elems])
    worst = 0.0
    for U in elems:
        for g in gens:
            R = to_SO3(g @ U)
            dmin = np.min(np.sqrt(np.sum((Rs - R) ** 2, axis=(1, 2))))
            worst = max(worst, dmin)
    # Trennschaerfe: minimaler Abstand zwischen verschiedenen Elementen
    n = len(elems)
    sep = np.inf
    for a in range(n):
        for b in range(a + 1, n):
            sep = min(sep, np.sqrt(np.sum((Rs[a] - Rs[b]) ** 2)))
    return worst, sep

def proj_order(U, maxn=200):
    P = U.copy()
    c0 = canon(I2)
    for n in range(1, maxn + 1):
        if canon(P) == c0:
            return n
        P = U @ P
    return None

# ---- Orbit von Z: Q = g^dagger Z g = n.sigma (spurfrei, hermitesch) ----
def bloch(Q):
    return np.array([0.5 * np.real(np.trace(SIG[i] @ Q)) for i in range(3)])

def orbit_of_Z(elems):
    vecs = []
    for g in elems:
        n = bloch(g.conj().T @ Z @ g)
        if not any(np.linalg.norm(n - v) < 1e-8 for v in vecs):
            vecs.append(n)
    return vecs

def is_pm_pauli_axes(vecs):
    """Lemma-1-Check: Orbit Teilmenge von +-{x,y,z}-Achsenvektoren (Toleranz TOL)."""
    axes = [np.array(e) * s for e in [(1., 0, 0), (0, 1., 0), (0, 0, 1.)] for s in (+1, -1)]
    dev = 0.0
    for v in vecs:
        dmin = min(np.linalg.norm(v - a) for a in axes)
        dev = max(dev, dmin)
    return dev < TOL, dev

# ---- K3 exakt ----
def K3_ops(Q1, Q2, Q3):
    """Lueders: M = (1/2)({Q2,Q1}+{Q3,Q2}-{Q3,Q1}); max_rho K3 = lambda_max(M).
    Rueckgabe: (lambda_max, Eigenwert-Spread) - Spread ~0 <=> rho-unabhaengig (M ~ c*I)."""
    A = lambda P, R: P @ R + R @ P
    M = 0.5 * (A(Q2, Q1) + A(Q3, Q2) - A(Q3, Q1))
    ev = np.linalg.eigvalsh(M)
    return ev[-1], ev[-1] - ev[0]

def run_k(k):
    s1, s2, F = braid_data(k)
    out = {"k": k}
    # Sanity: Unitaritaet + YBE
    out["sanity_YBE"] = float(np.max(np.abs(s1 @ s2 @ s1 - s2 @ s1 @ s2)))
    out["sanity_unitarity"] = float(max(np.max(np.abs(s1 @ s1.conj().T - I2)),
                                        np.max(np.abs(s2 @ s2.conj().T - I2))))
    # 1) Gruppe schliessen
    elems, closed = enumerate_group([s1, s2])
    worst, sep = closure_check(elems, [s1, s2, s1.conj().T, s2.conj().T])
    out["group_closed"] = bool(closed)
    out["group_order_projective"] = len(elems)
    out["closure_residual"] = float(worst)          # muss < TOL
    out["min_element_separation"] = float(sep)      # muss >> TOL
    out["sigma1_proj_order"] = proj_order(s1)
    # 2) Orbit von Z
    vecs = orbit_of_Z(elems)
    out["orbit_size"] = len(vecs)
    ok_axes, dev = is_pm_pauli_axes(vecs)
    out["orbit_is_pm_XYZ"] = bool(ok_axes)
    out["orbit_axes_deviation"] = float(dev)
    out["orbit_vectors"] = [[round(float(x), 6) for x in v] for v in vecs]
    # 3) K3-Maximum exakt: alle Paare (Q2,Q3) in OxO, Q1 = Z fest
    Qs = [v[0] * X + v[1] * Y + v[2] * Z for v in vecs]
    k3max, spread_max, arg = -np.inf, 0.0, None
    for i2, Q2 in enumerate(Qs):
        for i3, Q3 in enumerate(Qs):
            lam, spread = K3_ops(Z, Q2, Q3)
            spread_max = max(spread_max, spread)
            if lam > k3max:
                k3max, arg = lam, (i2, i3)
    out["K3_max_pairs_Q1eqZ"] = float(k3max)
    out["K3_argmax_bloch"] = {"n2": [round(float(x), 6) for x in vecs[arg[0]]],
                              "n3": [round(float(x), 6) for x in vecs[arg[1]]]}
    out["M_eigen_spread_max"] = float(spread_max)   # ~0 => K3-Operator ~ c*I => rho-unabhaengig
    # Zweite, unabhaengige Route: Bloch-Skalar K3 = n1.n2 + n2.n3 - n1.n3
    z = np.array([0., 0., 1.])
    k3_bloch = max(float(np.dot(z, n2) + np.dot(n2, n3) - np.dot(z, n3))
                   for n2 in vecs for n3 in vecs)
    out["K3_max_bloch_route"] = k3_bloch
    out["routes_agree"] = bool(abs(k3_bloch - k3max) < 1e-12)
    # Staerker: alle drei Zeitpunkte frei (Q1,Q2,Q3 in O^3)
    k3_tri = max(float(np.dot(n1, n2) + np.dot(n2, n3) - np.dot(n1, n3))
                 for n1 in vecs for n2 in vecs for n3 in vecs)
    out["K3_max_triples_free_Q1"] = k3_tri
    # 5) Staerkere Aussage: C in {0,+-1} fuer JEDES Paar (rho-unabhaengig ist es immer,
    #    da {n.sig, m.sig} = 2(n.m) I; Frage ist ob n.m in {0,+-1})
    allv = vecs + [z]
    cvals = sorted({round(float(np.dot(a, b)), 6) for a in allv for b in allv})
    out["correlator_values"] = cvals
    out["correlator_in_0pm1"] = bool(all(min(abs(c - t) for t in (-1.0, 0.0, 1.0)) < TOL
                                         for c in cvals))
    return out

def m2_seal_check():
    """Anker an die deponierte Referenz: max M2 ueber k=4-Gruppe x Stabilizer-Zustaende (soll 0.5585)."""
    s1, s2, _ = braid_data(4)
    elems, _ = enumerate_group([s1, s2])
    STAB = [np.array([1, 0], complex), np.array([0, 1], complex),
            np.array([1, 1], complex) / np.sqrt(2), np.array([1, -1], complex) / np.sqrt(2),
            np.array([1, 1j], complex) / np.sqrt(2), np.array([1, -1j], complex) / np.sqrt(2)]
    PAULIS = [I2, X, Y, Z]
    def M2(psi):
        psi = psi / np.linalg.norm(psi)
        return -np.log2(sum(np.real(np.vdot(psi, P @ psi)) ** 4 for P in PAULIS) / 2.0)
    return max(M2(U @ s) for U in elems for s in STAB)

if __name__ == "__main__":
    results = {}
    for k in [2, 4, 8]:
        r = run_k(k)
        results[f"k={k}"] = r
        print(f"--- k={k} ---")
        for key, val in r.items():
            if key != "k":
                print(f"  {key}: {val}")
    m2 = m2_seal_check()
    results["m2_seal_check_k4"] = float(m2)
    print(f"--- Anker: max M2 (k=4-Gruppe x Stab) = {m2:.4f} (B-Siegel: 0.5585)")

    # Befund-Logik
    k4 = results["k=4"]
    counterexample = k4["K3_max_triples_free_Q1"] > 1 + TOL or k4["K3_max_pairs_Q1eqZ"] > 1 + TOL
    results["befund"] = {
        "k4_counterexample_found": bool(counterexample),
        "k4_K3_max": k4["K3_max_triples_free_Q1"],
        "k2_K3_max": results["k=2"]["K3_max_triples_free_Q1"],
        "k8_K3_max": results["k=8"]["K3_max_triples_free_Q1"],
        "k8_exceeds_1": bool(results["k=8"]["K3_max_triples_free_Q1"] > 1 + TOL),
        "tol": TOL,
    }
    print("=== BEFUND ===")
    print(json.dumps(results["befund"], indent=2))
    with open("ergebnis_b1.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print("geschrieben: ergebnis_b1.json")
