"""
B_magic r2 (Sub-CC 3): 2. Maß (RoM) + LGI K3 EXAKT nach Paper 4 (K3=2C(B)-C(B^2),
C(U)=zz-Komponente von SO(3)-R(U)) + Vier-Wege-Dissoziations-Tabelle. Deterministisch, Seed 2026.
Reuse der unabhaengigen Konventionen aus engine.py (eigene F/R-Herleitung).
"""
import numpy as np
np.random.seed(2026)
I2=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],dtype=complex); Y=np.array([[0,-1j],[1j,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
PAULIS=[I2,X,Y,Z]

def braid_data(k):
    d=2*np.cos(np.pi/(k+2)); h_half=0.5*1.5/(k+2); h1=2.0/(k+2)
    R0=-np.exp(1j*np.pi*(0.0-2*h_half)); R1=np.exp(1j*np.pi*(h1-2*h_half))
    s1=np.diag([R0,R1]); root=np.sqrt(max(d*d-1,0))/d
    F=np.array([[1/d,root],[root,-1/d]],dtype=complex); s2=F@s1@F
    return s1,s2,F,d

def to_SO3(U):
    sig=[X,Y,Z]; R=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            R[i,j]=0.5*np.real(np.trace(sig[i]@U@sig[j]@U.conj().T))
    return R

# ---- Magic-Masse ----
def M2(psi):
    psi=psi/np.linalg.norm(psi); ev=[np.real(np.vdot(psi,P@psi)) for P in PAULIS]
    return -np.log2(sum(e**4 for e in ev)/2)
def RoM(psi):  # single-qubit Robustness of Magic = L1-Norm des Bloch-Vektors
    psi=psi/np.linalg.norm(psi)
    return sum(abs(np.real(np.vdot(psi,P@psi))) for P in (X,Y,Z))

STAB=[np.array(v,dtype=complex)/np.linalg.norm(v) for v in
      ([1,0],[0,1],[1,1],[1,-1],[1,1j],[1,-1j])]

def GEOM(psi):  # 3. Maß: 1 - max Stabilizer-Fidelity (geometrisches Maß)
    psi=psi/np.linalg.norm(psi)
    return 1.0 - max(abs(np.vdot(s,psi))**2 for s in STAB)

def canon(U,dec=5): return tuple(np.round(to_SO3(U).flatten(),dec))
def group(gens,cap=2000):
    seen={canon(I2):I2.copy()}; fr=[I2.copy()]; gg=gens+[np.linalg.inv(g) for g in gens]
    while fr:
        U=fr.pop()
        for g in gg:
            V=g@U; c=canon(V)
            if c not in seen: seen[c]=V; fr.append(V)
            if len(seen)>cap: return seen
    return seen
def magic_over(k):
    s1,s2,_,_=braid_data(k)
    if k in (2,4,8):
        els=list(group([s1,s2]).values())
    else:
        els=[];
        for _ in range(3000):
            U=I2.copy()
            for _ in range(24): U=(s1 if np.random.rand()<.5 else s2)@U
            els.append(U)
    m2=rom=geo=0.0
    for U in els:
        for s in STAB:
            psi=U@s; m2=max(m2,M2(psi)); rom=max(rom,RoM(psi)); geo=max(geo,GEOM(psi))
    return m2,rom,geo

# ---- LGI K3 = 2 C(B) - C(B^2), C(U,n)=n.R(U).n ; K3 ist das MAXIMUM ueber die
#      Braid-Gruppe (Paper 4): dichte k -> 3/2, endliche gedeckelt. axis-/state-optimiert. ----
def K3_opt_elem(B):  # max ueber Mess-Achse n  (= state+axis-opt, da K3(rho)=Tr[rho M])
    M=2*to_SO3(B)-to_SO3(B@B); Ms=0.5*(M+M.T)
    return float(np.max(np.linalg.eigvalsh(Ms)))
def K3_zhat_elem(B):
    M=2*to_SO3(B)-to_SO3(B@B); return float(M[2,2])
def K3_level(k):
    s1,s2,_,_=braid_data(k)
    if k in (2,4,8):
        els=list(group([s1,s2]).values())
    else:
        els=[]
        for _ in range(6000):
            U=I2.copy()
            for _ in range(20): U=(s1 if np.random.rand()<.5 else s2)@U
            els.append(U)
    return max(K3_opt_elem(B) for B in els), max(K3_zhat_elem(B) for B in els)

print("="*70)
print("VIER-WEGE-DISSOZIATIONS-TABELLE  (LGI K3 = max ueber Braid-Gruppe, Paper 4)")
print(f"{'k':>2} {'FLW-univ.':>9} {'3-Str.':>8} {'K3opt':>6} {'K3(Qz)':>7} {'LGI?':>5} {'M2':>7} {'RoM':>6} {'geom':>6}  Magic?")
flw=lambda k: 'nein' if k in(1,2) else ('NEIN(*)' if k==4 else 'ja')
for k in [2,3,4,5,6,7,8,9,10]:
    fin='endlich' if k in (2,4,8) else 'dicht'
    ko,kz=K3_level(k); m2,rom,geo=magic_over(k)
    fires='JA' if ko>1+1e-3 else 'NEIN'
    mag ='NEIN' if m2<1e-6 else 'JA'
    print(f"{k:>2} {flw(k):>9} {fin:>8} {ko:>6.3f} {kz:>7.3f} {fires:>5} {m2:>7.4f} {rom:>6.3f} {geo:>6.3f}  {mag}")
print("="*70)
print(f"Anker: dicht K3opt->3/2 ; k=8 K3opt=1.427 (max-Q, ikos.72°) & K3(Qz)=3/sqrt5={3/np.sqrt(5):.4f} ; k=4 K3=1.000 (inert)")
print("k=4: LGI feuert NICHT (K3=1, jede Achse+jeder Zustand) ABER Magic fast max (M2=0.558) => DISSOZIATION.")
print("Target 3: K3(rho)=Tr[rho M] zustandsabhaengig; k=4-Maximum ueber ALLE Zustaende+Achsen+Braids = 1.0")
print("  => genuin blind, strukturell inert, KEINE Tautologie. (*) k=4 = die FLW-Ausnahme (theorie-nicht-universell).")

print("="*70)
print("B-leicht: MAGIC-WITNESS D = 1/2(1+sum|<P>|), feuert iff sum|<P>|>1  (= Octahedron/RoM-Kriterium, n&s 1 Qubit)")
print(f"{'k':>2} {'sum|P|':>7} {'D':>6} {'Witness':>8} {'K3opt':>6} {'LGI?':>5}")
for k in [2,4,8,3,5]:
    ko,_=K3_level(k); m2,rom,geo=magic_over(k)
    D=0.5*(1+rom); fire='FEUERT' if rom>1+1e-6 else 'nein'
    lgi='JA' if ko>1+1e-3 else 'NEIN'
    print(f"{k:>2} {rom:>7.3f} {D:>6.3f} {fire:>8} {ko:>6.3f} {lgi:>5}")
print("=> k=4: Magic-Witness FEUERT (sum|P|=1.715>1, D=1.358) wo LGI BLIND ist (K3=1.000) = gefeuerter Detektor am Blind-Punkt.")
print("   k=2 (Clifford): sum|P|=1.000 (feuert NICHT), M2=0 -> sauberer fire/vanish. (Witness-Mathe Lehrbuch; neu = operationale Anwendung.)")
