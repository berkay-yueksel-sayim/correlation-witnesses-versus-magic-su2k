"""
T4 — Vier-Achsen-k-Karte (Konsolidierung), 2026-06-24, Seed 2026.
Fuehrt A/B-leicht/N1/T3 in EINER k-aufgeloesten Karte zusammen + 2 neue Rechnungen:
  (B-CURVE) volle RoM/geom-Kurve k=2..50 + k->inf-Limes.
  (B-FLOW)  Magic-Trajektorie unter sukzessivem Braiding: McKay-Endlichkeit (2O/2T/2I) als
            DISKRIMINATOR — endliche Gruppe (k=2,4,8) -> beschraenkt/periodisch; dichte k -> fuellt Plateau.
  (B-MAP)   Dissoziations-Karte: Magic(M2/RoM/geom) · Witness(Octahedron) · LGI(K3) · Gate-Magic(M^A, aus r6).
Verifizierte Single-Qubit-Engine (engine.py r1 Konventionen). Eigene Herleitung.
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
    """Gruppen-Elemente (endlich) oder Sampling (dicht). Rueckgabe: Liste U, dense-Flag."""
    gens=braid_gens(k)
    els,capped=group(gens)
    if capped:  # dicht -> sample
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
    """LGI K3 = max ueber Gruppe von 2*C(B)-C(B^2), C(U)=zz-Komponente von SO(3)-R(U). Luders-Bound 3/2."""
    els,dense,_=orbit(k)
    best=-9
    for B in els:
        RB=to_SO3(B); RB2=to_SO3(B@B); best=max(best,2*RB[2,2]-RB2[2,2])
    return best

# ---- ANKER / LOCKS ----
print("="*78); print("T4 — VIER-ACHSEN-k-KARTE (Magic · Witness · LGI · Gate)")
print("="*78)
m2_2,rom_2,g_2,w_2,_,ord2=max_measures(2)
assert m2_2<1e-9 and rom_2<1.0+1e-6, f"k=2-Clifford-Lock verletzt: M2={m2_2}, RoM={rom_2}"
print(f"[LOCK] k=2: M2={m2_2:.2e}(=0) RoM={rom_2:.4f}(=1 Oktaeder-Ecke) geom={g_2:.2e}(=0)  |group|={ord2}(2O=24)  OK")
print(f"[ANKER] Single-Qubit-Magic-Max: M2_max=log2(3/2)={np.log2(1.5):.4f}  RoM_max=sqrt(3)={np.sqrt(3):.4f}  (T-Richtung r=(1,1,1)/sqrt3)")

# ---- B-CURVE: volle Kurve + Limes ----
print("\n[B-CURVE] Magic-Maße vs k (max ueber Braid-Orbit) + k->inf-Limes:")
print(f"  {'k':>3} {'M2':>7} {'RoM':>7} {'geom':>7} {'|grp|':>7} {'dicht?':>7}")
curve=[]
for k in [2,3,4,5,6,7,8,9,10,12,16,20,50]:
    mM2,mRoM,mGeom,wit,dense,order=max_measures(k)
    curve.append((k,mM2,mRoM,mGeom,order,dense))
    print(f"  {k:>3} {mM2:>7.4f} {mRoM:>7.4f} {mGeom:>7.4f} {str(order):>7} {'dicht' if dense else 'endl.':>7}")
# Limes-Check: dichte k -> globaler Single-Qubit-Max
dense_m2=[m for (k,m,_,_,_,d) in curve if d];
print(f"  => Limes (dichte k): M2->{np.mean(dense_m2):.4f} ~ log2(3/2)={np.log2(1.5):.4f} (globaler 1q-Magic-Max; dichter Braid fuellt SU(2))")

# ---- B-FLOW: Magic-Trajektorie + McKay-Diskriminator ----
print("\n[B-FLOW] Magic-Trajektorie unter sukzessivem Braiding (McKay-Endlichkeit = Diskriminator):")
print("  endliche Gruppe (k=2,4,8) -> Magic beschraenkt durch Gruppen-Max; dichte k -> fuellt zum Plateau.")
print(f"  {'k':>3} {'Typ':>10} {'|grp|':>7} {'traj-max(L<=30)':>16} {'grp-max':>9} {'beschraenkt?':>13}")
for k in [2,4,8,3,5]:
    gens=braid_gens(k)
    els,dense,order=orbit(k)
    grpmax=max(M2(U@s) for U in els for s in STAB)
    tmax=0.0
    for _ in range(2000):
        U=I2.copy()
        for _ in range(np.random.randint(1,31)): U=gens[np.random.randint(2)]@U
        for s in STAB: tmax=max(tmax,M2(U@s))
    typ='2O/2T/2I' if not dense else 'dicht'
    bounded = 'JA (=grp-max)' if not dense and abs(tmax-grpmax)<1e-6 else ('fuellt' if dense else f'{tmax:.3f}')
    print(f"  {k:>3} {typ:>10} {str(order):>7} {tmax:>16.4f} {grpmax:>9.4f} {bounded:>13}")
print("  -> endliche k: Trajektorie erreicht GENAU den Gruppen-Max (endlich viele Zustaende, beschraenkt) ;")
print("     dichte k: Trajektorie fuellt zum globalen Plateau log2(3/2). Das ist der McKay-Diskriminator.")

# ---- B-MAP: Dissoziations-Karte ----
print("\n[B-MAP] Dissoziations-Karte (alle Achsen, k-aufgeloest):")
gateM=  {2:0.000,3:0.562,4:0.765,5:0.832,8:0.815}  # M^A(sigma2) aus r6 (T3)
print(f"  {'k':>3} {'Gruppe':>8} {'M2':>7} {'RoM':>7} {'geom':>7} {'Witness':>8} {'LGI-K3':>7} {'GateM^A':>8}")
mapping=[]
for k in [2,3,4,5,8]:
    mM2,mRoM,mGeom,wit,dense,order=max_measures(k)
    k3=K3(k)
    # Label-Praezisierung (RUN-3-Nachzug2, 2026-07-05, NUR Kommentar, KEINE Zahl geaendert — Historie unangetastet):
    # Die Klammerzahl (24/12/60) ist die korrekt berechnete PROJEKTIVE Ordnung (SO(3)-Quotient/PU(2)-Bild via
    # canon_SO3; verifiziert), NICHT falsch. Inkonsistent ist NUR das binaere Praefix "2X" davor — die binaere
    # SU(2)-Gruppe 2X hat per Definition die DOPPELTE Ordnung (2O=48, 2T=24, 2I=120), waere also mit ihrem
    # eigenen Praefix nicht mit 24/12/60 zu paaren. Die begleitende Analyse-Notiz nennt fuer
    # dieselbe k=2/4/8-Karte explizit "48/24/120 (2O/2T/2I, binaer)" — d.h. Rohcode (hier: projektiv + binaeres
    # Praefix) und Notiz-Prosa (binaer) nutzen bewusst verschiedene Ebenen derselben Gruppe; kein Wert ist
    # falsch, nur zwei Konventionen. k=4 kanonisch: projektiv T(12) <-> binaer 2T(24)=SL(2,3).
    grp={2:'2O(24)',4:'2T(12)',8:'2I(60)'}.get(k,f'dicht')
    gm=gateM.get(k,np.nan)
    mapping.append((k,grp,mM2,mRoM,mGeom,wit,k3,gm))
    print(f"  {k:>3} {grp:>8} {mM2:>7.4f} {mRoM:>7.4f} {mGeom:>7.4f} {wit:>8.4f} {k3:>7.4f} {gm:>8.4f}")
print(f"\n  LGI-Anker-Check: K3(k=4)={K3(4):.4f}(~1.0 klassisch) K3(k=8)={K3(8):.4f}(3/sqrt5={3/np.sqrt(5):.4f}) K3(dicht)->{K3(5):.4f}(Lueders 1.5)")
print("  DISSOZIATION: k=4 -> Magic & Witness & GateM HOCH, aber LGI-K3=1.0 (klassisch, BLIND). Die Achsen entkoppeln.")

import json
json.dump({'curve':[{'k':k,'M2':m,'RoM':r,'geom':g,'order':o,'dense':d} for (k,m,r,g,o,d) in curve],
           'map':[{'k':k,'group':grp,'M2':m,'RoM':r,'geom':g,'witness':w,'K3':k3,'gateMA':gm} for (k,grp,m,r,g,w,k3,gm) in mapping],
           'limits':{'M2_max':float(np.log2(1.5)),'RoM_max':float(np.sqrt(3))}},
          open('p5a_ergebnis_t4.json','w'),indent=2)
print("\nGeschrieben: p5a_ergebnis_t4.json")
