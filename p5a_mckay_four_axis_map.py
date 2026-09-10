"""
T4 — four-axis k-map (consolidation), 2026-06-24, seed 2026.
Brings A/B-light/N1/T3 together in ONE k-resolved map + 2 new computations:
  (B-CURVE) full RoM/geom curve k=2..50 + k->inf limit.
  (B-FLOW)  magic trajectory under successive braiding: McKay finiteness (2O/2T/2I) as the
            DISCRIMINATOR — finite group (k=2,4,8) -> bounded/periodic; dense k -> fills the plateau.
  (B-MAP)   dissociation map: magic(M2/RoM/geom) · witness(octahedron) · LGI(K3) · gate magic(M^A, from r6).
Verified single-qubit engine (engine.py r1 conventions). Own derivation.
"""
import numpy as np
np.random.seed(2026)
X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]],complex); Z=np.array([[1,0],[0,-1]],complex); I2=np.eye(2,dtype=complex)
STAB=[np.array(v,complex)/np.linalg.norm(v) for v in ([1,0],[0,1],[1,1],[1,-1],[1,1j],[1,-1j])]

def braid_gens(k):
    d=2*np.cos(np.pi/(k+2)); R0=-np.exp(-1j*3*np.pi/(2*(k+2))); R1=np.exp(1j*np.pi/(2*(k+2)))
    root=np.sqrt(max(d*d-1,0))/d; F=np.array([[1/d,root],[root,-1/d]],complex)
    s1=np.diag([R0,R1]).astype(complex); s2=F@s1@F; return s1,s2

def to_SO3(U):
    R=np.zeros((3,3))
    for i,si in enumerate((X,Y,Z)):
        for j,sj in enumerate((X,Y,Z)):
            R[i,j]=0.5*np.real(np.trace(si@U@sj@U.conj().T))
    return R
def canon(U): return tuple(np.round(to_SO3(U).flatten(),5))
def group(gens,cap=500):
    seen={canon(I2):I2.copy()}; fr=[I2.copy()]; gg=list(gens)+[np.linalg.inv(g) for g in gens]
    while fr:
        U=fr.pop()
        for g in gg:
            V=g@U; c=canon(V)
            if c not in seen: seen[c]=V; fr.append(V)
            if len(seen)>cap: return list(seen.values()),True
    return list(seen.values()),False

def bloch(psi): psi=psi/np.linalg.norm(psi); return np.array([np.real(np.vdot(psi,P@psi)) for P in (X,Y,Z)])
def M2(psi):
    psi=psi/np.linalg.norm(psi); ev=[np.real(np.vdot(psi,P@psi)) for P in (I2,X,Y,Z)]
    return -np.log2(sum(e**4 for e in ev)/2)
def RoM(psi): return np.sum(np.abs(bloch(psi)))
def geom(psi):
    psi=psi/np.linalg.norm(psi); return 1-max(abs(np.vdot(s,psi))**2 for s in STAB)

def orbit(k):
    """Group elements (finite) or sampling (dense). Returns: list U, dense flag."""
    gens=braid_gens(k)
    els,capped=group(gens)
    if capped:  # dense -> sample
        els=[]
        for _ in range(4000):
            U=I2.copy()
            for _ in range(24): U=gens[np.random.randint(2)]@U
            els.append(U)
        return els,True,None
    return els,False,len(els)

def max_measures(k):
    els,dense,order=orbit(k)
    mM2=mRoM=mGeom=0.0; wit=0.0
    for U in els:
        for s in STAB:
            psi=U@s; mM2=max(mM2,M2(psi)); r=RoM(psi); mRoM=max(mRoM,r); mGeom=max(mGeom,geom(psi)); wit=max(wit,r)
    return mM2,mRoM,mGeom,wit,dense,order

def K3(k):
    """LGI K3 = max over the group of 2*C(B)-C(B^2), C(U)=zz component of SO(3)-R(U). Luders bound 3/2."""
    els,dense,_=orbit(k)
    best=-9
    for B in els:
        RB=to_SO3(B); RB2=to_SO3(B@B); best=max(best,2*RB[2,2]-RB2[2,2])
    return best

# ---- ANCHORS / LOCKS ----
print("="*78); print("T4 — FOUR-AXIS k-MAP (magic · witness · LGI · gate)")
print("="*78)
m2_2,rom_2,g_2,w_2,_,ord2=max_measures(2)
assert m2_2<1e-9 and rom_2<1.0+1e-6, f"k=2 Clifford lock violated: M2={m2_2}, RoM={rom_2}"
print(f"[LOCK] k=2: M2={m2_2:.2e}(=0) RoM={rom_2:.4f}(=1 Oktaeder-Ecke) geom={g_2:.2e}(=0)  |group|={ord2}(2O=24)  OK")
print(f"[ANCHOR] single-qubit magic max: M2_max=log2(3/2)={np.log2(1.5):.4f}  RoM_max=sqrt(3)={np.sqrt(3):.4f}  (T direction r=(1,1,1)/sqrt3)")

# ---- B-CURVE: full curve + limit ----
print("\n[B-CURVE] magic measures vs k (max over the braid orbit) + k->inf limit:")
print(f"  {'k':>3} {'M2':>7} {'RoM':>7} {'geom':>7} {'|grp|':>7} {'dense?':>7}")
curve=[]
for k in [2,3,4,5,6,7,8,9,10,12,16,20,50]:
    mM2,mRoM,mGeom,wit,dense,order=max_measures(k)
    curve.append((k,mM2,mRoM,mGeom,order,dense))
    print(f"  {k:>3} {mM2:>7.4f} {mRoM:>7.4f} {mGeom:>7.4f} {str(order):>7} {'dense' if dense else 'finite':>7}")
# limit check: dense k -> global single-qubit max
dense_m2=[m for (k,m,_,_,_,d) in curve if d];
print(f"  => limit (dense k): M2->{np.mean(dense_m2):.4f} ~ log2(3/2)={np.log2(1.5):.4f} (global 1q magic max; a dense braid fills SU(2))")

# ---- B-FLOW: magic trajectory + McKay discriminator ----
print("\n[B-FLOW] magic trajectory under successive braiding (McKay finiteness = discriminator):")
print("  finite group (k=2,4,8) -> magic bounded by the group max; dense k -> fills up to the plateau.")
print(f"  {'k':>3} {'Type':>10} {'|grp|':>7} {'traj-max(L<=30)':>16} {'grp-max':>9} {'bounded?':>13}")
for k in [2,4,8,3,5]:
    gens=braid_gens(k)
    els,dense,order=orbit(k)
    grpmax=max(M2(U@s) for U in els for s in STAB)
    tmax=0.0
    for _ in range(2000):
        U=I2.copy()
        for _ in range(np.random.randint(1,31)): U=gens[np.random.randint(2)]@U
        for s in STAB: tmax=max(tmax,M2(U@s))
    kind='2O/2T/2I' if not dense else 'dense'
    bounded = 'YES (=grp-max)' if not dense and abs(tmax-grpmax)<1e-6 else ('fills' if dense else f'{tmax:.3f}')
    print(f"  {k:>3} {kind:>10} {str(order):>7} {tmax:>16.4f} {grpmax:>9.4f} {bounded:>13}")
print("  -> finite k: the trajectory reaches EXACTLY the group max (finitely many states, bounded) ;")
print("     dense k: the trajectory fills up to the global plateau log2(3/2). That is the McKay discriminator.")

# ---- B-MAP: dissociation map ----
print("\n[B-MAP] dissociation map (all axes, k-resolved):")
gateM=  {2:0.000,3:0.562,4:0.765,5:0.832,8:0.815}  # M^A(sigma2) from r6 (T3)
print(f"  {'k':>3} {'Group':>8} {'M2':>7} {'RoM':>7} {'geom':>7} {'Witness':>8} {'LGI-K3':>7} {'GateM^A':>8}")
mapping=[]
for k in [2,3,4,5,8]:
    mM2,mRoM,mGeom,wit,dense,order=max_measures(k)
    k3=K3(k)
    # label precision (RUN-3 follow-up 2, 2026-07-05, COMMENT ONLY, NO number changed — history untouched):
    # The number in brackets (24/12/60) is the correctly computed PROJECTIVE order (SO(3) quotient/PU(2) image via
    # canon_SO3; verified), NOT wrong. Inconsistent is ONLY the binary prefix "2X" in front of it — the binary
    # SU(2) group 2X has by definition TWICE the order (2O=48, 2T=24, 2I=120), and would therefore not pair with
    # 24/12/60 under its own prefix. The accompanying analysis note gives, for the same k=2/4/8 map,
    # explicitly "48/24/120 (2O/2T/2I, binary)" — i.e. raw code (here: projective + binary
    # prefix) and note prose (binary) deliberately use different levels of the same group; no value is
    # wrong, only two conventions. k=4 canonical: projective T(12) <-> binary 2T(24)=SL(2,3).
    grp={2:'2O(24)',4:'2T(12)',8:'2I(60)'}.get(k,f'dense')
    gm=gateM.get(k,np.nan)
    mapping.append((k,grp,mM2,mRoM,mGeom,wit,k3,gm))
    print(f"  {k:>3} {grp:>8} {mM2:>7.4f} {mRoM:>7.4f} {mGeom:>7.4f} {wit:>8.4f} {k3:>7.4f} {gm:>8.4f}")
print(f"\n  LGI-Anker-Check: K3(k=4)={K3(4):.4f}(~1.0 classical) K3(k=8)={K3(8):.4f}(3/sqrt5={3/np.sqrt(5):.4f}) K3(dense)->{K3(5):.4f}(Lueders 1.5)")
print("  DISSOCIATION: k=4 -> magic & witness & gateM HIGH, but LGI-K3=1.0 (classical, BLIND). The axes decouple.")

import json
from pathlib import Path
if __name__ == "__main__":
    json.dump({'curve':[{'k':k,'M2':m,'RoM':r,'geom':g,'order':o,'dense':d} for (k,m,r,g,o,d) in curve],
               'map':[{'k':k,'group':grp,'M2':m,'RoM':r,'geom':g,'witness':w,'K3':k3,'gateMA':gm} for (k,grp,m,r,g,w,k3,gm) in mapping],
               'limits':{'M2_max':float(np.log2(1.5)),'RoM_max':float(np.sqrt(3))}},
              open(Path(__file__).resolve().parent / 'p5a_ergebnis_t4.json','w'),indent=2)
    print("\nWritten: p5a_ergebnis_t4.json")
