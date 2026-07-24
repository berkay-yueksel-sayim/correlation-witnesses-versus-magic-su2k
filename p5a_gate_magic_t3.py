"""
T3 — GATE-/Channel-Magic des SU(2)_k-Single-Qubit-Braids, 2026-06-24, Seed 2026.
Bisher: nur ZUSTANDS-Magic. NEU: die Magie des GATTERS U(sigma) selbst, via amortisierte
Stabilizer-Renyi-Entropie (Zhu, Chen, Shen, Liu, Yu, Wang, arXiv:2409.06659 — Volltext-VERIFIZIERT).

DEFINITION (Zhu et al.):
  M_a^A(U) := sup_m  max_{|phi> in H_{n+m}}  [ M_a((U (x) I_2^m)|phi>) - M_a(|phi>) ]
  (Output-SRE minus Input-SRE, ueber Ancilla-Zahl m + reine Eingaenge; "amortisiert").
  Monoton fuer alle a>=0; faithful: M^A(U)=0 <=> U Clifford; T-count-Schranke t(U) >= M^A(U)/M^A(T).
  M_2^A(T) = 2 - log2(3) ~ 0.4150 (verifizierter Anker).

EHRLICHE LABELS (Pflicht, aus dem Volltext-Check):
- Wir rechnen eine UNTERE SCHRANKE (m in {0,1}; sup_m ist unbeschraenkt -> exakter Wert nur >=).
- Nur der SINGLE-QUBIT-Braid (sigma = 2x2, fuer alle k wohldefiniertes 1-Qubit-Gatter). Der
  Multi-Qubit-Fusionsraum (k!=2 nicht 2er-Potenz) braucht eine Encoding-Bridge -> hier NICHT.
- Anwendung des Maszes von Zhu et al. auf U(sigma) (kein im Paper belegtes Braid-Resultat).
Tracer/Locks: k=2 (Ising-Braid = Clifford) -> M^A=0 ; T-Gate -> 0.4150 ; H/S (Clifford) -> 0.
"""
import numpy as np
from scipy.optimize import minimize
np.random.seed(2026)

X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]],complex); Z=np.array([[1,0],[0,-1]],complex); I2=np.eye(2,dtype=complex)
P1=[I2,X,Y,Z]
def paulis(n):
    mats=[np.array([[1]],complex)]
    for _ in range(n):
        mats=[np.kron(M,P) for M in mats for P in P1]
    return mats
PA={n:paulis(n) for n in (1,2)}

def M2(psi):
    n=int(round(np.log2(len(psi)))); d=2**n
    s=sum((np.vdot(psi,M@psi).real)**4 for M in PA[n])
    return -np.log2(s/d)

def amortized_lb(U, m=0, restarts=40):
    """untere Schranke fuer M_2^A(U): max ueber reine (1+m)-Qubit-Eingaenge von M2(out)-M2(in)."""
    dim=2**(1+m); Uf=np.kron(U,np.eye(2**m,dtype=complex))
    def neg(x):
        v=x[:dim]+1j*x[dim:]; nrm=np.linalg.norm(v)
        if nrm<1e-12: return 0.0
        v=v/nrm
        return -(M2(Uf@v)-M2(v))
    best=0.0
    for _ in range(restarts):
        x0=np.random.randn(2*dim)
        r=minimize(neg,x0,method='Nelder-Mead',options={'maxiter':4000,'xatol':1e-8,'fatol':1e-10})
        best=max(best,-r.fun)
    return best

def amortized(U, restarts=40):
    return max(amortized_lb(U,0,restarts), amortized_lb(U,1,restarts))

# ---- SU(2)_k Single-Qubit-Braid-Generatoren (verifizierte Konventionen, engine.py r1) ----
def braid_gens(k):
    d=2*np.cos(np.pi/(k+2)); R0=-np.exp(-1j*3*np.pi/(2*(k+2))); R1=np.exp(1j*np.pi/(2*(k+2)))
    root=np.sqrt(max(d*d-1,0))/d; F=np.array([[1/d,root],[root,-1/d]],complex)
    s1=np.diag([R0,R1]).astype(complex); s2=F@s1@F
    return s1,s2

# ---- ANKER / LOCKS ----
print("="*64); print("T3 — GATE-MAGIC (amortized SRE, untere Schranke m<=1)")
print("="*64)
T=np.diag([1,np.exp(1j*np.pi/4)]).astype(complex)
H=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
S=np.diag([1,1j]).astype(complex)
mT=amortized(T); mH=amortized(H); mS=amortized(S)
print(f"[ANKER] M_2^A(T)   = {mT:.4f}   (Soll 2-log2(3)={2-np.log2(3):.4f})  {'OK' if abs(mT-(2-np.log2(3)))<2e-2 else 'CHECK'}")
print(f"[LOCK ] M_2^A(H)   = {mH:.4f}   (Clifford -> 0)  {'OK' if mH<2e-2 else 'CHECK'}")
print(f"[LOCK ] M_2^A(S)   = {mS:.4f}   (Clifford -> 0)  {'OK' if mS<2e-2 else 'CHECK'}")
s1_2,s2_2=braid_gens(2)
mk2=max(amortized(s1_2),amortized(s2_2))
print(f"[LOCK ] M_2^A(sigma,k=2) = {mk2:.4f}  (Ising-Braid = Clifford -> 0)  {'OK' if mk2<3e-2 else 'CHECK'}")

# ---- GATE-MAGIC-KURVE ueber k ----
print("\n[KURVE] Gate-Magic des Single-Qubit-Braids vs k (untere Schranke):")
print(f"  {'k':>2} {'M^A(s1)':>9} {'M^A(s2)':>9} {'M^A(s1s2)':>10} {'T-count>=':>10}")
rows=[]
for k in [2,3,4,5,8]:
    s1,s2=braid_gens(k)
    m1=amortized(s1); m2=amortized(s2); m12=amortized(s2@s1)
    tc=m12/(2-np.log2(3))
    rows.append((k,m1,m2,m12,tc))
    print(f"  {k:>2} {m1:>9.4f} {m2:>9.4f} {m12:>10.4f} {tc:>10.3f}")
print("\n  Labels: untere Schranke (m<=1); Single-Qubit-Braid; Anwendung von Zhu et al. 2409.06659.")
print("  Lock k=2->0 = Ising-Clifford-Konsistenz. T-count>= = M^A(word)/M^A(T) (untere Schranke).")
import json
json.dump({'anchors':{'T':float(mT),'H':float(mH),'S':float(mS),'sigma_k2':float(mk2),'M2A_T_exact':float(2-np.log2(3))},
           'curve':[{'k':k,'M_A_s1':float(a),'M_A_s2':float(b),'M_A_s1s2':float(c),'Tcount_lb':float(t)} for (k,a,b,c,t) in rows]},
          open('p5a_ergebnis_t3.json','w'),indent=2)
print("\nGeschrieben: p5a_ergebnis_t3.json")
