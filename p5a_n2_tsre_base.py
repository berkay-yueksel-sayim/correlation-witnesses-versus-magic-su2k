"""
N2 — magic of topological order.  2026-06-23, seed 2026.

GOAL: compute the stabilizer Renyi entropy M_2 of the doubled-Fibonacci Levin-Wen GROUND STATE
on a concrete EDGE-QUBIT honeycomb torus (a real number), with the toric-code(Z2) GS on the
SAME geometry as a falsifiable SRE=0 null control.

GROUNDING (not guessed — real project sources read):
- Geometry (vertex incidence, plaquette hex cycles): adapted from the exact-diagonalization
  reference implementation of the doubled-Fibonacci Levin-Wen honeycomb torus
  (sim04e, GSD=4 pinned in advance + confirmed to machine precision).
- Levin-Wen B_p: B_p=(1/D^2)(1 + phi*B_p^tau), D^2=phi+2; B_p^tau = 6-F product (LW2005 Eq.13,
  NO inner sum). F^{ttt}_t=[[1/phi,1/sqrt(phi)],[1/sqrt(phi),-1/phi]] (PROJEKT.md).
  Ordering of the 6 F factors EXACTLY as in sim04e shotA build_cycle_operator (verified).
- M_2 enumerator + anchors (M2(|0>)=0, M2(|+>)=0, M2(T|+>)=log2(4/3), additivity, GHZ/Bell=0):
  project B_magic r1-r3 (verified). M_2(psi)=-log2( (1/2^E) sum_{P in 4^E} <psi|P|psi>^4 ).

METHOD: B_p (Eq.13) + Q_v (branching) on the full 2^E edge-qubit space; H=-sum_v Q_v -
sum_p B_p exactly diagonalized (eigh in the Q_v-allowed subspace), GS = joint +1 eigenspace,
ONE GS vector -> exact M_2 (4^E Paulis, fast iterated enumerator). Toric: same
geometry, A_v=Prod X / B_p=Prod Z stabilizers -> stabilizer GS -> M_2 EXACTLY 0 (kill switch).

PRIOR ART (never "first"): Wei&Liu arXiv:2503.04566 (doubled Fib = LRM phase, qualitatively via
GSD/Bravyi-Konig); Korbany-Ellison-Stephen-Piroli arXiv:2605.22424 (LRM via mutual information);
1D TSRE subtraction Nehra et al. arXiv:2512.16673. Contribution here: the FIRST explicitly COMPUTED
M_2 number of the real-space edge-qubit string-net GS + toric-SRE=0 control + a reusable lock.
"""
import numpy as np
import itertools, time, json, sys
from pathlib import Path

np.random.seed(2026)
PHI = (1 + np.sqrt(5)) / 2.0
INV_PHI = 1.0 / PHI
SQ_PHI = np.sqrt(PHI)
D2 = PHI + 2.0
OUT = Path(__file__).resolve().parent

# verifizierte Geometrie importieren (read-only; deposited alongside as
# p5a_n2_hex_cycle_generator.py)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p5a_n2_hex_cycle_generator import (A_id, B_id, edge_x, edge_y, edge_z,  # noqa
                                 get_plaquette_cycle, get_edge_endpoints,
                                 check_plaquette_cleanliness)

# =====================================================================
# 0. Pauli + EXAKTER M_2-Enumerator (iteriert, schnell) + Anker
# =====================================================================
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [I2, X, Y, Z]

def pauli_expectations_fast(psi, E):
    """Robust + fast: <psi|P|psi> for all 4^E P (real). Attach one Pauli axis (4) per qubit.
       Layout invariant: T ALWAYS has the form (2,)*E (out-ket axes, order of qubits
       0..E-1) followed by the Pauli axes attached so far, (4,)*done, at the END.
       For qubit q contract Pmat[in] with out axis q (position q), new out axis q + new
       Pauli axis at the end. Finally contract the out axes with the bra."""
    psi = (psi / np.linalg.norm(psi)).astype(complex)
    Pmat = np.stack(PAULI)                       # (4,2,2): [p,out,in]
    T = psi.reshape((2,) * E)                    # only out axes, no Pauli axis
    for q in range(E):
        # T axes: [out_0..out_{E-1}] (E of them) + [pauli_0..pauli_{q-1}] (q of them).
        # contract Pmat 'in' (axis 2) with out axis q (at position q).
        T = np.tensordot(T, Pmat, axes=([q], [2]))
        # tensordot removes out axis q and appends Pmat's remaining axes (p[4], out'[2]) AT THE END.
        # New form: [out_0..out_{q-1}, out_{q+1}..out_{E-1}, pauli_0..pauli_{q-1}, p_q, out'_q]
        # out'_q (last axis) must go back to position q (out block). p_q stays at the Pauli end,
        # BUT currently sits BEFORE out'_q -> after moving out'_q to q the order is fine.
        T = np.moveaxis(T, -1, q)     # out'_q ans out-Block-Ende-Slot q
        # now the last axis = p_q (at the end of the Pauli block) -> correct.
    # T form: (2,)*E [out axes 0..E-1] + (4,)*E [Pauli axes 0..E-1].
    bra = psi.conj().reshape(-1)
    T = T.reshape(2**E, 4**E)
    ev = bra @ T                                   # (4^E,) komplex
    return ev.real

def M2_exact(psi, E, max_axes=6):
    """Exact M_2. Direct for small E; chunked for large E so that the intermediate tensor
       (4^inner x 2^E) fits into RAM. Accumulates sum_P <P>^4 over outer Pauli chunks."""
    if E <= max_axes:
        ev = pauli_expectations_fast(psi, E)
        return -np.log2(np.sum(ev**4) / (2.0**E))
    # chunked: fixed Paulis on the first c qubits (outer 4^c), the rest enumerated exactly
    c = E - max_axes                      # # outer-Qubits
    s = _m2_sum_chunked(psi, E, c)
    return -np.log2(s / (2.0**E))

def _m2_sum_chunked(psi, E, c):
    """sum_{P in 4^E} <psi|P|psi>^4, split: outer = Paulis on qubits 0..c-1 (4^c),
       inner = exact 4^(E-c) enumeration via pauli_expectations_fast on the reduced
       operator. <P>=<psi| (P_out (x) P_in) |psi>. Apply P_out fixed -> phi_out (E-qubit
       vector), then inner enumeration of the Paulis on qubits c..E-1 of <psi| P_in |phi_out>."""
    psi = (psi/np.linalg.norm(psi)).astype(complex)
    total = 0.0
    rest = E - c
    Pmat = np.stack(PAULI)
    nout = 4**c; done = 0; t0 = time.perf_counter()
    for combo_out in itertools.product(range(4), repeat=c):
        # phi = (P_out (x) I) |psi>
        phi = psi.reshape((2,)*E)
        for q, p in enumerate(combo_out):
            if p == 0:
                continue
            phi = np.tensordot(PAULI[p], phi, axes=([1], [q]))
            phi = np.moveaxis(phi, 0, q)
        phi = phi.reshape(-1)
        # inner: all 4^rest expectation values <psi| (I_out (x) P_in) |phi>
        # = bra(psi) . (P_in on qubits c..E-1) . phi   ; outer qubits are the identity.
        ev = _inner_expectations(psi.conj(), phi, E, c)   # (4^rest,)
        total += np.sum(ev**4)
    return total

def _inner_expectations(bra, phi, E, c):
    """<bra| (I on 0..c-1) (x) (P on c..E-1) |phi> for all 4^rest Paulis. rest=E-c."""
    rest = E - c
    Pmat = np.stack(PAULI)
    T = phi.reshape((2,)*E)
    # build Pauli axes only for qubits c..E-1
    for q in range(c, E):
        T = np.tensordot(T, Pmat, axes=([q], [2]))
        T = np.moveaxis(T, -1, q)
    # T-Form: (2,)*E [out] + (4,)*rest [pauli c..E-1]
    T = T.reshape(2**E, 4**rest)
    ev = bra.reshape(-1) @ T
    return ev.real

def expval_one(psi, combo, E):
    """Single <psi|P_combo|psi> (for sampling)."""
    t = psi.reshape((2,) * E)
    out = t
    for q, p in enumerate(combo):
        if p == 0:
            continue
        out = np.tensordot(PAULI[p], out, axes=([1], [q]))
        out = np.moveaxis(out, 0, q)
    return np.vdot(psi.reshape(-1), out.reshape(-1)).real

def M2_sample(psi, E, nsamp, rng):
    psi = psi / np.linalg.norm(psi)
    vals = np.empty(nsamp)
    for i in range(nsamp):
        combo = tuple(rng.integers(0, 4, size=E))
        vals[i] = expval_one(psi, combo, E)**4
    mean = vals.mean(); sem = vals.std(ddof=1)/np.sqrt(nsamp)
    arg = (2.0**E)*mean
    m2 = -np.log2(arg)
    dm2 = abs((2.0**E)*sem/(arg*np.log(2)))
    return m2, dm2

def anchor_locks():
    res = {}
    z0 = np.array([1, 0], dtype=complex)
    plus = np.array([1, 1], dtype=complex)/np.sqrt(2)
    Tg = np.diag([1, np.exp(1j*np.pi/4)])
    Tp = Tg @ plus
    res['M2(|0>)'] = M2_exact(z0, 1)
    res['M2(|+>)'] = M2_exact(plus, 1)
    res['M2(T|+>)'] = M2_exact(Tp, 1)
    res['log2(4/3)'] = float(np.log2(4/3))
    res['M2(T+ (x) T+)'] = M2_exact(np.kron(Tp, Tp), 2)
    res['2*log2(4/3)'] = float(2*np.log2(4/3))
    ghz = np.zeros(8, dtype=complex); ghz[0] = ghz[7] = 1/np.sqrt(2)
    bell = np.zeros(4, dtype=complex); bell[0] = bell[3] = 1/np.sqrt(2)
    res['M2(GHZ3)'] = M2_exact(ghz, 3)
    res['M2(Bell)'] = M2_exact(bell, 2)
    return res

# =====================================================================
# 1. geometry (imported, verified), generic N x M
# =====================================================================
def build_geo(N, M):
    E = 3*N*M
    NV = 2*N*M
    # vertex incidence exactly as in hilbert_2x2.vertex_incident_edges
    def vinc(v):
        NM = N*M
        if v < NM:
            i, j = v//M, v % M
            return [edge_x(i, j, N, M), edge_y(i, j, N, M), edge_z(i, j, N, M)]
        b = v-NM; i, j = b//M, b % M
        return [edge_x(i, j, N, M), edge_y((i+1) % N, j, N, M), edge_z(i, (j+1) % M, N, M)]
    vert_edges = {v: vinc(v) for v in range(NV)}
    # Plaquetten + Sauberkeit
    plaqs = []
    clean_all = True
    for i in range(N):
        for j in range(M):
            ok, reason, _ = check_plaquette_cleanliness(i, j, N, M)
            edges, verts = get_plaquette_cycle(i, j, N, M)
            if not ok or len(set(edges)) != 6:
                clean_all = False
                continue
            ring_set = set(edges)
            legs = []
            legok = True
            for v in verts:
                non = [q for q in vert_edges[v] if q not in ring_set]
                if len(non) != 1:
                    legok = False; break
                legs.append(non[0])
            if not legok:
                clean_all = False; continue
            plaqs.append(dict(ring=edges, ring_v=verts, legs=legs))
    return dict(N=N, M=M, E=E, NV=NV, vert_edges=vert_edges, plaqs=plaqs, clean=clean_all)

# =====================================================================
# 2. F-Symbol + B_p^tau (LW Eq.13)
# =====================================================================
_F2 = np.array([[INV_PHI, 1/SQ_PHI], [1/SQ_PHI, -INV_PHI]])
def allowed_triple(a, b, c):
    return (a+b+c) != 1            # Fibonacci: forbidden = exactly one tau

def f_symbol(i, j, m, k, l, n):
    if not (allowed_triple(i, j, m) and allowed_triple(k, l, m)
            and allowed_triple(i, n, l) and allowed_triple(j, k, n)):
        return 0.0
    if i == 1 and j == 1 and k == 1 and l == 1:
        return _F2[m, n]
    return 1.0

def build_Bp_tau(geo, p):
    E = geo['E']
    ring = geo['plaqs'][p]['ring']; legs = geo['plaqs'][p]['legs']
    dim = 2**E
    Bt = np.zeros((dim, dim))
    for col in range(dim):
        cfg = [(col >> (E-1-q)) & 1 for q in range(E)]
        old = [cfg[ring[k]] for k in range(6)]
        legl = [cfg[legs[k]] for k in range(6)]
        for new in itertools.product((0, 1), repeat=6):
            amp = 1.0
            for k in range(6):
                amp *= f_symbol(legl[k], old[k-1], old[k], 1, new[k], new[(k-1) % 6])
                if amp == 0.0:
                    break
            if amp == 0.0:
                continue
            ncfg = cfg.copy()
            for k in range(6):
                ncfg[ring[k]] = new[k]
            row = 0
            for q in range(E):
                row |= ncfg[q] << (E-1-q)
            Bt[row, col] += amp
    return Bt

def Qv_mask(geo):
    E = geo['E']; dim = 2**E
    idxs = np.arange(dim)
    bits = ((idxs[:, None] >> (E-1-np.arange(E))[None, :]) & 1)
    allowed = np.ones(dim, dtype=bool)
    for v, qs in geo['vert_edges'].items():
        s = bits[:, qs].sum(axis=1)
        allowed &= (s != 1)
    return allowed

# =====================================================================
# 3. Fibonacci-GS via ED
# =====================================================================
def fib_groundstate(geo, verbose=True):
    E = geo['E']; dim = 2**E
    t0 = time.perf_counter()
    allowed = Qv_mask(geo)
    sub = np.where(allowed)[0]; subdim = len(sub)
    H = np.zeros((subdim, subdim)); Bps = []
    for p in range(len(geo['plaqs'])):
        Bt = build_Bp_tau(geo, p)[np.ix_(sub, sub)]
        Bp = (np.eye(subdim) + PHI*Bt)/D2
        Bps.append(Bp); H -= Bp
    Hs = (H+H.T)/2
    evals, evecs = np.linalg.eigh(Hs)
    E0 = evals[0]
    deg = int(np.sum(evals <= E0+1e-8))
    gap = float(evals[deg]-evals[deg-1]) if subdim > deg else float('nan')
    gs_sub = evecs[:, :deg]
    gs_full = np.zeros((dim, deg), dtype=complex); gs_full[sub, :] = gs_sub
    bp_dev = 0.0
    for Bp in Bps:
        for c in range(deg):
            bp_dev = max(bp_dev, np.max(np.abs(Bp @ gs_sub[:, c] - gs_sub[:, c])))
    info = dict(E=E, full_dim=dim, sub_dim=subdim, E0=float(E0), GSD=deg,
                gap=gap, bp_dev=float(bp_dev), build_s=time.perf_counter()-t0)
    if verbose:
        print(f"  [Fib]   E={E} subdim={subdim} | E0={E0:.6f} GSD={deg} "
              f"gap={gap:.4f} max|Bp*gs-gs|={bp_dev:.1e} ({info['build_s']:.1f}s)")
    return gs_full, info

# =====================================================================
# 4. toric-code GS (Z2) on the same geometry
# =====================================================================
def toric_groundstate(geo):
    E = geo['E']; dim = 2**E
    idxs = np.arange(dim)
    bits = ((idxs[:, None] >> (E-1-np.arange(E))[None, :]) & 1)
    H = np.zeros((dim, dim))
    for p in geo['plaqs']:
        zval = np.prod(1-2*bits[:, p['ring']], axis=1)
        H[idxs, idxs] -= zval
    for v, qs in geo['vert_edges'].items():
        mask = 0
        for q in qs:
            mask |= 1 << (E-1-q)
        perm = idxs ^ mask
        H[perm, idxs] -= 1.0
    Hs = (H+H.T)/2
    evals, evecs = np.linalg.eigh(Hs)
    E0 = evals[0]; deg = int(np.sum(evals <= E0+1e-8))
    gs = evecs[:, :deg].astype(complex)
    info = dict(E=E, E0=float(E0), GSD=deg,
                gap=float(evals[deg]-evals[deg-1]) if dim > deg else float('nan'))
    return gs, info

# =====================================================================
# 5. M_2 of a GS subspace (one deterministic vector)
# =====================================================================
def m2_of_gs(gs_full, E, exact_ceiling=12, nsamp=300000):
    psi = gs_full[:, 0].copy()
    k = np.argmax(np.abs(psi))
    psi = psi*np.exp(-1j*np.angle(psi[k]))
    imag = float(np.linalg.norm(psi.imag))
    if E <= exact_ceiling:
        return dict(method='exact', M2=float(M2_exact(psi, E)), E=E, imag_norm=imag)
    rng = np.random.default_rng(2026)
    m2, dm2 = M2_sample(psi, E, nsamp, rng)
    return dict(method='sample', M2=float(m2), dM2=float(dm2), nsamp=nsamp, E=E, imag_norm=imag)

# =====================================================================
# MAIN
# =====================================================================
def main():
    report = {'seed': 2026, 'phi': PHI, 'D2': D2}
    print("="*72)
    print("N2 — magic of topological order (doubled Fibonacci vs toric)")
    print("="*72)

    print("\n[LOCK] M_2 enumerator against anchors:")
    anc = anchor_locks()
    for kk in ['M2(|0>)', 'M2(|+>)', 'M2(T|+>)', 'log2(4/3)', 'M2(T+ (x) T+)',
               '2*log2(4/3)', 'M2(GHZ3)', 'M2(Bell)']:
        print(f"    {kk:>16} = {anc[kk]:.10f}")
    assert abs(anc['M2(|0>)']) < 1e-9
    assert abs(anc['M2(|+>)']) < 1e-9
    assert abs(anc['M2(T|+>)']-anc['log2(4/3)']) < 1e-9
    assert abs(anc['M2(T+ (x) T+)']-anc['2*log2(4/3)']) < 1e-9
    assert abs(anc['M2(GHZ3)']) < 1e-9 and abs(anc['M2(Bell)']) < 1e-9
    print("    -> Enumerator-Anker OK")
    report['anchors'] = {k: float(v) for k, v in anc.items()}
    report['lock_product_state_zero'] = float(anc['M2(|0>)'])
    report['lock_known_magic_state'] = {'M2(T|+>)': float(anc['M2(T|+>)']),
                                        'log2(4/3)': float(anc['log2(4/3)'])}

    report['sizes'] = []
    headline = None
    for (N, M) in [(2, 2)]:
        geo = build_geo(N, M)
        print(f"\n[Geometrie] {N}x{M}: E={geo['E']} Kanten, {len(geo['plaqs'])} Plaquetten, "
              f"{geo['NV']} Vertices, clean={geo['clean']}")
        if len(geo['plaqs']) == 0:
            print("  no well-formed plaquettes -> skipped")
            continue
        gs_fib, info_fib = fib_groundstate(geo)
        gs_tor, info_tor = toric_groundstate(geo)
        print(f"  [Toric] E0={info_tor['E0']:.4f} GSD={info_tor['GSD']} gap={info_tor['gap']:.4f}")
        t0 = time.perf_counter()
        m2_fib = m2_of_gs(gs_fib, geo['E'])
        m2_tor = m2_of_gs(gs_tor, geo['E'])
        dt = time.perf_counter()-t0
        print(f"  M_2(Fib-GS)   = {m2_fib['M2']:.6f}  ({m2_fib['method']}, imag={m2_fib['imag_norm']:.1e})")
        print(f"  M_2(Toric-GS) = {m2_tor['M2']:.3e}  ({m2_tor['method']})   [M2-Zeit {dt:.1f}s]")
        entry = dict(N=N, M=M, E=geo['E'], plaqs=len(geo['plaqs']),
                     fib=dict(info=info_fib, m2=m2_fib), tor=dict(info=info_tor, m2=m2_tor))
        report['sizes'].append(entry)
        if info_fib['GSD'] == 4 and headline is None:
            headline = entry

    if headline is None and report['sizes']:
        headline = report['sizes'][-1]
    if headline:
        report['headline_size'] = (headline['N'], headline['M'])
        report['lock_toric_stabilizer_zero'] = headline['tor']['m2']['M2']
        report['lock_groundstate_genuine'] = dict(
            GSD=headline['fib']['info']['GSD'], E0=headline['fib']['info']['E0'],
            bp_dev=headline['fib']['info']['bp_dev'])
        report['headline_M2_fib'] = headline['fib']['m2']['M2']

    print("\n" + "="*72)
    print("PFLICHT-LOCKS (Zahlen)")
    print("="*72)
    print(f"  product_state_zero : M2(|0>)      = {report['lock_product_state_zero']:.2e}")
    print(f"  known_magic_state  : M2(T|+>)     = {report['lock_known_magic_state']['M2(T|+>)']:.6f} "
          f"(log2(4/3)={report['lock_known_magic_state']['log2(4/3)']:.6f})")
    if headline:
        print(f"  toric_stabilizer_0 : M2(Toric-GS)= {report['lock_toric_stabilizer_zero']:.2e}  "
              f"(Geometrie {report['headline_size']})")
        g = report['lock_groundstate_genuine']
        print(f"  groundstate_genuine: GSD={g['GSD']} E0={g['E0']:.6f} max|Bp*gs-gs|={g['bp_dev']:.1e}")
        print(f"  >>> HEADLINE: M_2(doubled-Fib GS) = {report['headline_M2_fib']:.6f}  "
              f"vs  M_2(Toric GS) = {report['lock_toric_stabilizer_zero']:.2e}")

    with open(OUT/"result.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)
    print(f"\nWritten: {OUT/'result.json'}")
    return report

if __name__ == "__main__":
    main()
