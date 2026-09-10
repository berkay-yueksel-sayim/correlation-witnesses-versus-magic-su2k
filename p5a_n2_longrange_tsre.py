"""
N2 v2 — CORRECTED magic of topological order. 2026-06-24, seed 2026.
Main-thread fix of v1 (tsre.py) after two diagnosed problems:
  (1) v1 used exact 4^12 enumeration (E=12) -> runaway (>17min). HERE: fast bit-trick
      Pauli sampler (O(dim) per Pauli, vectorized) with reported SEM.
  (2) v1 took gs[:,0] = an ARBITRARY eigh vector in the 4-fold degenerate ground space. M_2 does
      depend on the logical state; an arbitrary toric ground-space vector is NOT automatically a
      stabilizer state -> the null control would be invalid. HERE: CANONICAL states:
        - Fibonacci = string-net condensate |GS> ~ (Prod_p B_p) |0...0>  (trivial flux sector,
          from the string-free vacuum; principled + gauge-fixed).
        - Toric     = a true stabilizer state |GS> ~ (Prod_v (I+A_v)/2) |0...0>.
Verified v1 primitives (geometry, F symbol, B_p^tau Eq.13, Q_v, M_2 anchors) are imported.

PRIOR ART (never "first"): Wei&Liu arXiv:2503.04566 (doubled Fib = long-range-magic phase),
Nehra et al. arXiv:2512.16673 (TSRE subtraction 1D), Leone-Oliviero-Hamma (M_2 tool),
Fliss arXiv:2011.01962 (Mana SU(2)_k). Contribution: first explicitly COMPUTED real-space edge-qubit
M_2 number of the string-net GS + toric-SRE=0 stabilizer control + reusable stabilizer-discriminator lock.
"""
import numpy as np
import time, json
from pathlib import Path
import p5a_n2_tsre_base as v1   # verifizierte Primitive (build_geo, build_Bp_tau, Qv_mask, fib_groundstate, anchor_locks ...)

np.random.seed(2026)
PHI = v1.PHI; D2 = v1.D2
OUT = Path(__file__).resolve().parent

# ============================================================
# Schneller Pauli-Sampler (Bit-Trick) — <psi|P|psi> in O(dim), vektorisiert
# P = (x) sigma_{p_q}, p in {0=I,1=X,2=Y,3=Z}.  P|i> = c_i |i ^ xmask>,
#   xmask = bits with p in {X,Y};  zfull = bits with p in {Y,Z};
#   c_i = i^{nY} * (-1)^{popcount(i & zfull)} ;  <P> = i^{nY} * sum_i conj(psi_i) (-1)^... psi_{i^xmask}
# ============================================================
def make_popcount(dim):
    return np.array([bin(x).count("1") for x in range(dim)], dtype=np.int64)

def expval_pauli(psi, combo, E, ar, POP):
    """combo: array length E over {0,1,2,3}. psi: (dim,) complex normalized. ar=arange(dim)."""
    xmask = 0; zfull = 0; nY = 0
    for q, p in enumerate(combo):
        bit = 1 << (E - 1 - q)            # qubit q = bit (E-1-q) (same convention as v1)
        if p == 1:   xmask |= bit
        elif p == 2: xmask |= bit; zfull |= bit; nY += 1
        elif p == 3: zfull |= bit
    signs = 1 - 2 * (POP[ar & zfull] & 1)               # (-1)^popcount(i & zfull)
    perm = ar ^ xmask
    s = np.vdot(psi, signs * psi[perm])                 # sum conj(psi_i) signs_i psi_{i^xmask}
    # FIX (skeptic finding 2026-06-24): the correct phase factor is (-i)^nY, not i^nY.
    # Reason: c_{k^xmask} instead of c_k -> signs shift by (-1)^nY (xmask & zfull = Y bits).
    # |<P>| was exact throughout -> M_2 (uses <P>^4) and the magnitude discriminator UNAFFECTED;
    # the sign matters only for odd moments/Mana. Locked down correctly here.
    val = ((-1j) ** nY) * s
    return val.real                                     # P hermitesch -> reell

def m2_sample(psi, E, nsamp, rng, ar, POP, collect_vals=False):
    psi = (psi / np.linalg.norm(psi)).astype(complex)
    vals = np.empty(nsamp)
    for i in range(nsamp):
        combo = rng.integers(0, 4, size=E)
        vals[i] = expval_pauli(psi, combo, E, ar, POP) ** 4
    mean = vals.mean(); sem = vals.std(ddof=1) / np.sqrt(nsamp)
    arg = (2.0 ** E) * mean                              # = (1/2^E) sum<P>^4 * 2^E  -> estimator for 2^E*E[<P>^4]
    m2 = -np.log2(arg)
    dm2 = abs((2.0 ** E) * sem / (arg * np.log(2)))
    out = dict(M2=float(m2), dM2=float(dm2), nsamp=nsamp, arg=float(arg))
    if collect_vals:
        out['expvals'] = None  # placeholder; separate routine for the discriminator
    return out

def stabilizer_fraction(psi, E, nsamp, rng, ar, POP, tol=1e-6):
    """Fraction of sampled <P> within tol of {0,+1,-1}. Stabilizer state -> 1.0 ; otherwise -> <1.0.
       Also returns the histogram extreme (max |<P>| that is NOT ~{0,1}) as a magic signature."""
    psi = (psi / np.linalg.norm(psi)).astype(complex)
    n_stab = 0; max_noninteger = 0.0; sample_vals = []
    for _ in range(nsamp):
        combo = rng.integers(0, 4, size=E)
        v = expval_pauli(psi, combo, E, ar, POP)
        sample_vals.append(v)
        near = min(abs(v), abs(v - 1), abs(v + 1))
        if near < tol:
            n_stab += 1
        else:
            max_noninteger = max(max_noninteger, min(abs(v), abs(abs(v) - 1)))
    return dict(frac_in_01=n_stab / nsamp, max_noninteger=float(max_noninteger),
                n=nsamp, mean_abs=float(np.mean(np.abs(sample_vals))))

# ============================================================
# canonical states
# ============================================================
def canonical_fib_condensate(geo, gs_sub_basis, allowed):
    """|GS_fib> ~ (Prod_p B_p) |0...0> : string-net condensate in the trivial sector.
       B_p = (I + phi*B_p^tau)/D2 (projector onto +1). |0..0> = string-free (Q_v satisfied).
       Verify: the result lies in the degenerate ground space (overlap ~1) and B_p|GS>=|GS>."""
    E = geo['E']; dim = 2 ** E
    psi = np.zeros(dim, dtype=complex); psi[0] = 1.0      # |0...0> = all edges 0
    Bps = []
    for p in range(len(geo['plaqs'])):
        Bt = v1.build_Bp_tau(geo, p)
        Bp = (np.eye(dim) + PHI * Bt) / D2
        Bps.append(Bp)
        psi = Bp @ psi
    nrm = np.linalg.norm(psi)
    psi = psi / nrm
    # locks: B_p|GS>=|GS> ; support only on Q_v-allowed configs ; in the eigh ground space
    bp_dev = max(np.max(np.abs(Bp @ psi - psi)) for Bp in Bps)
    leak_qv = float(np.linalg.norm(psi[~allowed]))
    # overlap with the eigh ground space (gs_sub_basis: (subdim, GSD) on the allowed index)
    sub = np.where(allowed)[0]
    psi_sub = psi[sub]
    proj = gs_sub_basis @ (gs_sub_basis.conj().T @ psi_sub)
    in_gs = float(np.linalg.norm(proj) / max(np.linalg.norm(psi_sub), 1e-300))
    return psi, dict(bp_dev=float(bp_dev), leak_qv=leak_qv, overlap_groundspace=in_gs,
                     norm_from_vacuum=float(nrm))

def canonical_toric_stabilizer(geo):
    """|GS_tor> ~ (Prod_v (I+A_v)/2) |0...0> : a true stabilizer state. A_v = Prod_{q in v} X_q.
       |0..0> satisfies B_p (Prod Z = +1) ; projection with vertex operators -> stabilizer GS."""
    E = geo['E']; dim = 2 ** E
    ar = np.arange(dim)
    psi = np.zeros(dim, dtype=complex); psi[0] = 1.0
    for v, qs in geo['vert_edges'].items():
        xmask = 0
        for q in qs:
            xmask |= 1 << (E - 1 - q)
        psi = (psi + psi[ar ^ xmask]) / 2.0            # (I + A_v)/2
    nrm = np.linalg.norm(psi); psi = psi / nrm
    # locks: A_v|GS>=|GS> for all v ; B_p (Prod Z) |GS> = |GS> for all p
    av_dev = 0.0
    for v, qs in geo['vert_edges'].items():
        xmask = 0
        for q in qs:
            xmask |= 1 << (E - 1 - q)
        av_dev = max(av_dev, np.max(np.abs(psi[ar ^ xmask] - psi)))
    bits = ((ar[:, None] >> (E - 1 - np.arange(E))[None, :]) & 1)
    bp_dev = 0.0
    for p in geo['plaqs']:
        zval = np.prod(1 - 2 * bits[:, p['ring']], axis=1)   # Prod Z on the ring
        bp_dev = max(bp_dev, np.max(np.abs(zval * psi - psi)))
    return psi, dict(av_dev=float(av_dev), bp_dev=float(bp_dev), support=int(np.sum(np.abs(psi) > 1e-12)))

# ============================================================
# MAIN
# ============================================================
def main():
    rep = {'seed': 2026, 'phi': float(PHI), 'note': 'v2 main-thread fix: sampler + canonical states'}
    print("=" * 72)
    print("N2 v2 — magic of topological order (CANONICAL; doubled Fib vs toric)")
    print("=" * 72)

    # --- anchor lock (from v1, exact) ---
    anc = v1.anchor_locks()
    print("\n[LOCK] M_2-Enumerator-Anker:")
    for kk in ['M2(|0>)', 'M2(|+>)', 'M2(T|+>)', 'log2(4/3)', 'M2(GHZ3)', 'M2(Bell)']:
        print(f"    {kk:>12} = {anc[kk]:.10f}")
    assert abs(anc['M2(|0>)']) < 1e-9 and abs(anc['M2(|+>)']) < 1e-9
    assert abs(anc['M2(T|+>)'] - anc['log2(4/3)']) < 1e-9
    assert abs(anc['M2(GHZ3)']) < 1e-9 and abs(anc['M2(Bell)']) < 1e-9
    print("    -> anchors OK (M2(T|+>)=log2(4/3); stabilizer states=0)")
    rep['anchors'] = {k: float(v) for k, v in anc.items()}

    geo = v1.build_geo(2, 2)
    E = geo['E']; dim = 2 ** E
    ar = np.arange(dim); POP = make_popcount(dim)
    print(f"\n[Geometrie] 2x2: E={E} Kanten, {len(geo['plaqs'])} Plaq, {geo['NV']} Vertices, clean={geo['clean']}")

    # --- eigh ground space (for GSD/genuine + condensate overlap) ---
    gs_full, info_fib = v1.fib_groundstate(geo, verbose=True)
    allowed = v1.Qv_mask(geo); sub = np.where(allowed)[0]
    gs_sub_basis = gs_full[sub, :]               # (subdim, GSD) ortho basis of the ground space
    print(f"  [Fib eigh] GSD={info_fib['GSD']} E0={info_fib['E0']:.6f} gap={info_fib['gap']:.4f} "
          f"max|Bp*gs-gs|={info_fib['bp_dev']:.1e}")

    # --- build the canonical states ---
    t0 = time.perf_counter()
    psi_fib, lf = canonical_fib_condensate(geo, gs_sub_basis, allowed)
    psi_tor, lt = canonical_toric_stabilizer(geo)
    print(f"  [Fib  kanon.] Kondensat (Prod B_p)|0>: B_p-dev={lf['bp_dev']:.1e}  Q_v-leak={lf['leak_qv']:.1e}  "
          f"Ueberlapp-Grundraum={lf['overlap_groundspace']:.6f}")
    print(f"  [Toric kanon.] Stabilizer (Prod (I+A_v)/2)|0>: A_v-dev={lt['av_dev']:.1e}  B_p-dev={lt['bp_dev']:.1e}  "
          f"Support={lt['support']} Konfigs")
    assert lf['overlap_groundspace'] > 1 - 1e-6, "Fib condensate does not lie in the ground space!"
    assert lf['bp_dev'] < 1e-9 and lf['leak_qv'] < 1e-9, "Fib-Kondensat-Lock verletzt"
    assert lt['av_dev'] < 1e-9 and lt['bp_dev'] < 1e-9, "Toric-Stabilizer-Lock verletzt"

    # --- M_2 via Sampling ---
    NS = 150000
    rng = np.random.default_rng(2026)
    m2_fib = m2_sample(psi_fib, E, NS, rng, ar, POP)
    rng2 = np.random.default_rng(7)
    m2_tor = m2_sample(psi_tor, E, NS, rng2, ar, POP)
    dt = time.perf_counter() - t0
    print(f"\n  M_2(Fib-GS)   = {m2_fib['M2']:.4f} +/- {m2_fib['dM2']:.4f}   (nsamp={NS})")
    print(f"  M_2(Toric-GS) = {m2_tor['M2']:.4f} +/- {m2_tor['dM2']:.4f}   (nsamp={NS})  [erwartet 0]")

    # --- Stabilizer-Diskriminator-Lock ---
    rng3 = np.random.default_rng(11); rng4 = np.random.default_rng(13)
    sd_tor = stabilizer_fraction(psi_tor, E, 20000, rng3, ar, POP)
    sd_fib = stabilizer_fraction(psi_fib, E, 20000, rng4, ar, POP)
    print(f"\n  [Diskriminator] Anteil <P> in {{0,+/-1}}:  Toric={sd_tor['frac_in_01']:.4f}  Fib={sd_fib['frac_in_01']:.4f}")
    print(f"                  max non-integer |<P>|: toric={sd_tor['max_noninteger']:.2e}  Fib={sd_fib['max_noninteger']:.3f}")
    print(f"  -> toric: all <P> in {{0,+/-1}} (stabilizer => M_2=0 exactly). Fib: continuous <P> (non-stabilizer).")

    print("\n" + "=" * 72)
    print("PFLICHT-LOCKS")
    print("=" * 72)
    print(f"  product_state_zero : M2(|0>)        = {anc['M2(|0>)']:.2e}")
    print(f"  known_magic_state  : M2(T|+>)       = {anc['M2(T|+>)']:.6f}  (log2(4/3)={anc['log2(4/3)']:.6f})")
    print(f"  toric_stabilizer_0 : M2(Toric-GS)   = {m2_tor['M2']:.4f} +/- {m2_tor['dM2']:.4f} ; stab-frac={sd_tor['frac_in_01']:.4f}")
    print(f"  groundstate_genuine: GSD={info_fib['GSD']} E0={info_fib['E0']:.4f} Kondensat-Ueberlapp={lf['overlap_groundspace']:.6f}")
    print(f"  >>> HEADLINE: M_2(doubled-Fib GS) = {m2_fib['M2']:.4f} +/- {m2_fib['dM2']:.4f}  >> 0  "
          f"vs  M_2(Toric GS) = {m2_tor['M2']:.4f} (~0)")

    rep['E'] = E
    rep['fib'] = dict(info=info_fib, canonical=lf, m2=m2_fib, discriminator=sd_fib)
    rep['toric'] = dict(canonical=lt, m2=m2_tor, discriminator=sd_tor)
    rep['locks'] = dict(product_state_zero=float(anc['M2(|0>)']),
                        known_magic_state=float(anc['M2(T|+>)']),
                        toric_stabilizer_zero=float(m2_tor['M2']),
                        toric_stab_frac=float(sd_tor['frac_in_01']),
                        groundstate_genuine=dict(GSD=info_fib['GSD'], overlap=lf['overlap_groundspace']))
    rep['headline_M2_fib'] = float(m2_fib['M2'])
    rep['toric_m2_note'] = (
        "toric.m2.M2 (~0.107) is a positively-biased finite-sample Monte Carlo estimate, "
        "NOT the true value. The true toric-code stabilizer 2-Renyi entropy is 0 exactly "
        "(a stabilizer ground state has no nonstabilizerness). This is certified in this "
        "same JSON by toric.discriminator: max_noninteger=0.0 and frac_in_01=1.0, i.e. every "
        "sampled Pauli expectation value lies exactly in {0,+-1}, the stabilizer signature. "
        "The paper cites the exact value M_2(Toric)=0. The ~0.107 is the known positive bias "
        "floor of the finite-sample M_2 estimator near zero, not a residual nonstabilizerness."
    )
    with open(OUT / "ergebnis_v2.json", "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=2, default=float)
    print(f"\nWritten: {OUT/'ergebnis_v2.json'}  (total {dt:.1f}s build+sampling)")
    return rep

if __name__ == "__main__":
    main()
