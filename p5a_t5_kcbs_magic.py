"""
T5 — contextuality ↔ magic on the SAME SU(2)_k anyon state, 2026-06-24, seed 2026.
KCBS needs d>=3 (a qubit is non-contextual) -> d=3 fusion space (n=4 strands, total=2, dim 3 = qutrit).
On the same qutrit: (i) KCBS pentagon witness (Klyachko), class. 2, quantum √5≈2.236 ; (ii) Mana(d=3, Gross-Wigner,
prime) as the magic measure. Question: do contextuality and magic couple over the braid orbit, k-resolved?

HEDGES (Volltext-Check 2026-06-24, PFLICHT):
- Chou arXiv:2506.14537 (peer-rev. MDPI) already has Fibonacci-braid KCBS contextuality (state-dependent) → NO priority claim
  for "anyon braiding shows KCBS". Our remaining contribution: a numerically reproducible number (Chou gives none) + k resolution +
  magic↔context coupling on the same state (Chou has no magic).
- Howard-Wallman-Veitch-Emerson 1401.4174 (context=magic) holds for qudits of ODD PRIME d → only analogy/motivation for our qutrit.
Pentagon oriented GENERICALLY (fixed random SO(3), seed 2026) — otherwise the axes align trivially. Own derivation.
"""
import sys

if __name__ == "__main__" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
np.random.seed(2026)

# ---------- SU(2)_k n-strand TLJ-braid (verified) ----------
def qint(m,k): return np.sin(m*np.pi/(k+2))/np.sin(np.pi/(k+2))
def braid_phases(k):
    d=2*np.cos(np.pi/(k+2)); h_half=0.5*1.5/(k+2); h1=2.0/(k+2)
    R0=-np.exp(1j*np.pi*(0.0-2*h_half)); R1=np.exp(1j*np.pi*(h1-2*h_half)); return d,R0,R1
def paths(n,k,total):
    res=[]
    def rec(p):
        if len(p)==n+1:
            if p[-1]==total: res.append(tuple(p))
            return
        for s in (+1,-1):
            q=p[-1]+s
            if 0<=q<=k: rec(p+[q])
    rec([0]); return res
def build_braid(n,k,total):
    d,R0,R1=braid_phases(k); B=sorted(paths(n,k,total)); idx={p:i for i,p in enumerate(B)}; D=len(B)
    gens=[]
    s1=np.zeros((D,D),complex)
    for p in B: s1[idx[p],idx[p]]=R0 if p[2]==0 else R1
    gens.append(s1)
    for m in range(2,n):
        em=np.zeros((D,D),complex)
        for p in B:
            if p[m-1]!=p[m+1]: continue
            base=p[m-1]; opts=[v for v in (base-1,base+1) if 0<=v<=k]
            for a in opts:
                for b in opts:
                    pa=list(p); pa[m]=a; pb=list(p); pb[m]=b
                    if tuple(pa) in idx and tuple(pb) in idx:
                        em[idx[tuple(pa)],idx[tuple(pb)]]=np.sqrt(qint(a+1,k)*qint(b+1,k))/qint(base+1,k)
        gens.append(R1*np.eye(D)+((R0-R1)/d)*em)
    return gens,B

# ---------- KCBS-Pentagon (Qutrit), generisch orientiert ----------
def rand_SO3(seed):
    rng=np.random.default_rng(seed); A=rng.normal(size=(3,3)); Q,R=np.linalg.qr(A); Q=Q@np.diag(np.sign(np.diag(R)))
    if np.linalg.det(Q)<0: Q[:,0]*=-1
    return Q
def kcbs_vectors():
    cosT=np.sqrt(np.cos(np.pi/5)/(1+np.cos(np.pi/5))); sinT=np.sqrt(1-cosT**2)
    V=[np.array([sinT*np.cos(4*np.pi*i/5), sinT*np.sin(4*np.pi*i/5), cosT]) for i in range(5)]
    R=rand_SO3(2026); return [ (R@v).astype(complex) for v in V ]
KV=kcbs_vectors()
def kcbs(psi):  # Sum_i |<v_i|psi>|^2 ; klass. <=2, quanten <=√5
    psi=psi/np.linalg.norm(psi); return sum(abs(np.vdot(v,psi))**2 for v in KV)

# ---------- Mana(d=3) Gross-Wigner ----------
def build_HW(d=3):
    w=np.exp(2j*np.pi/d); X=np.zeros((d,d),complex)
    for j in range(d): X[(j+1)%d,j]=1.0
    Z=np.diag([w**j for j in range(d)]); inv2=(d+1)//2
    def Dop(a,b): return (w**(-inv2*a*b))*np.linalg.matrix_power(Z,b)@np.linalg.matrix_power(X,a)
    A0=sum(Dop(a,b) for a in range(d) for b in range(d))/d
    A={(a,b):Dop(a,b)@A0@Dop(a,b).conj().T for a in range(d) for b in range(d)}
    return A
AW=build_HW(3)
def mana(psi):
    psi=psi/np.linalg.norm(psi); return np.log2(sum(abs(np.real(np.vdot(psi,Au@psi))/3) for Au in AW.values()))  # log2 base (series convention, RUN-3 follow-up 2 2026-07-05; was np.log/ln — normalization Sigma_u A(u)=d*I identical to r9, verified)

def canon(U):
    ph=U[np.argmax(np.abs(U[:,0])),0]; U=U*np.conj(ph)/abs(ph)
    return tuple(np.round(U.flatten(),4).view(float))
def group(gens,cap=600):
    seen={canon(np.eye(3,dtype=complex)):np.eye(3,dtype=complex)}; fr=[np.eye(3,dtype=complex)]
    gg=list(gens)+[np.linalg.inv(g) for g in gens]
    while fr:
        U=fr.pop()
        for g in gg:
            V=g@U; c=canon(V)
            if c not in seen: seen[c]=V; fr.append(V)
            if len(seen)>cap: return list(seen.values()),True
    return list(seen.values()),False

print("="*74); print("T5 — CONTEXTUALITY ↔ MAGIC on the same d=3 anyon qutrit (KCBS + Mana)")
print("="*74)
# Dim-Check + Locks
g4,B4=build_braid(4,8,2); print(f"Fusion space n=4,total=2: dim={len(B4)} (Qutrit) Pfade={B4}")
assert len(B4)==3, f"dim!={3}"
# LOCK: KCBS bounds — apex reaches √5, a generic comp.-basis state < √5
import itertools
apex=rand_SO3(2026)@np.array([0,0,1.0])  # the state that gives √5 at this orientation (KV=R@v_kanon -> apex=R@(0,0,1))
print(f"[LOCK] KCBS quanten-max (Apex) = {kcbs(apex.astype(complex)):.4f} (target √5={np.sqrt(5):.4f})")
print(f"[LOCK] KCBS comp. basis e0/e1/e2 = {[round(kcbs(np.eye(3,dtype=complex)[:,j]),3) for j in range(3)]} (generic, <√5)")
# LOCK: Mana stabilizer/comp. basis -> 0 ; random > 0
print(f"[LOCK] Mana(d=3) comp. basis = {[round(mana(np.eye(3,dtype=complex)[:,j]),3) for j in range(3)]} (=0) ; random={mana(np.random.randn(3)+1j*np.random.randn(3)):.3f}(>0)")

print("\n[K-RESOLVED] max KCBS + max Mana over the braid orbit + coupling:")
print(f"  {'k':>3} {'Type':>7} {'|grp|':>7} {'maxKCBS':>8} {'contextual?':>13} {'maxMana':>8} {'Mana@KCBSmax':>13} {'corr(K,Mana)':>13}")
rows=[]
for k in [2,4,8,3,5]:
    gens,B=build_braid(4,k,2)
    if len(B)!=3: print(f"  k={k}: dim={len(B)}!=3 skipped"); continue
    els,dense=group(gens)
    if dense:
        els=[]
        for _ in range(4000):
            U=np.eye(3,dtype=complex)
            for _ in range(22): U=gens[np.random.randint(len(gens))]@U
            els.append(U)
    Ks=[]; Ms=[]; states=[]
    for U in els:
        for j in range(3):
            psi=U[:,j]; states.append(psi); Ks.append(kcbs(psi)); Ms.append(mana(psi))
    Ks=np.array(Ks); Ms=np.array(Ms)
    maxK=Ks.max(); maxM=Ms.max(); mAtK=Ms[np.argmax(Ks)]
    corr=np.corrcoef(Ks,Ms)[0,1] if Ks.std()>1e-9 and Ms.std()>1e-9 else float('nan')
    ctx='YES (>2)' if maxK>2+1e-6 else 'no'
    kind='finite' if not dense else 'dense'
    rows.append((k,kind,len(els) if not dense else None,maxK,ctx,maxM,mAtK,corr))
    print(f"  {k:>3} {kind:>7} {str(len(els) if not dense else '∞'):>7} {maxK:>8.4f} {ctx:>13} {maxM:>8.4f} {mAtK:>13.4f} {corr:>13.3f}")

print("\n  READING: maxKCBS>2 = braiding reaches a contextual qutrit state; corr(K,Mana) = coupling of the axes.")
print("  HEDGES: NO priority claim for braid-KCBS (Chou 2506.14537 state-dep., peer-rev.); Howard 1401.4174 = qudit/odd-prime-d (analogy).")
import json
from pathlib import Path
if __name__ == "__main__":
    json.dump({'kcbs_classical':2,'kcbs_quantum':float(np.sqrt(5)),
               'rows':[{'k':k,'type':t,'maxKCBS':float(mk),'contextual':c,'maxMana':float(mm),'Mana_at_KCBSmax':float(ma),'corr':float(co) if co==co else None} for (k,t,_,mk,c,mm,ma,co) in rows]},
              open(Path(__file__).resolve().parent / 'p5a_ergebnis_t5.json','w'),indent=2)
    print("\nWritten: p5a_ergebnis_t5.json")
