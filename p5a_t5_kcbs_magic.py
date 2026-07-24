"""
T5 — Kontextualitaet ↔ Magic am SELBEN SU(2)_k-Anyon-State, 2026-06-24, Seed 2026.
KCBS braucht d>=3 (Qubit ist nicht-kontextuell) -> d=3-Fusionsraum (n=4 strands, total=2, dim 3 = Qutrit).
Auf demselben Qutrit: (i) KCBS-Pentagon-Witness (Klyachko), klass. 2, quanten √5≈2.236 ; (ii) Mana(d=3, Gross-Wigner,
prime) als Magic-Mass. Frage: koppeln Kontextualitaet und Magic ueber den Braid-Orbit, k-aufgeloest?

HEDGES (Volltext-Check 2026-06-24, PFLICHT):
- Chou arXiv:2506.14537 (peer-rev. MDPI) hat Fibonacci-Braid-KCBS-Kontextualitaet (state-dependent) bereits → KEINE Erstheit
  fuer "Anyon-Braiding zeigt KCBS". Unser Restbeitrag: numer. reproduzierbare Zahl (Chou zeigt keine) + k-Aufloesung +
  Magic↔Kontext-Kopplung am selben State (Chou hat null Magic).
- Howard-Wallman-Veitch-Emerson 1401.4174 (Kontext=Magic) gilt fuer qudits ODD PRIME d → nur Analogie/Motivation fuer unser Qutrit.
Pentagon GENERISCH orientiert (fixe Zufalls-SO(3), Seed 2026) — sonst triviale Achsen-Ausrichtung. Eigene Herleitung.
"""
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
    psi=psi/np.linalg.norm(psi); return np.log2(sum(abs(np.real(np.vdot(psi,Au@psi))/3) for Au in AW.values()))  # log2-Basis (Serien-Konvention, RUN-3-Nachzug2 2026-07-05; war np.log/ln — Normierung Sigma_u A(u)=d*I identisch zu r9, verifiziert)

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

print("="*74); print("T5 — KONTEXTUALITAET ↔ MAGIC am selben d=3-Anyon-Qutrit (KCBS + Mana)")
print("="*74)
# Dim-Check + Locks
g4,B4=build_braid(4,8,2); print(f"Fusionsraum n=4,total=2: dim={len(B4)} (Qutrit) Pfade={B4}")
assert len(B4)==3, f"dim!={3}"
# LOCK: KCBS-Bounds — Apex erreicht √5, generischer Komp.-Basiszustand < √5
import itertools
apex=rand_SO3(2026)@np.array([0,0,1.0])  # der State, der bei dieser Orientierung √5 gibt (KV=R@v_kanon -> Apex=R@(0,0,1))
print(f"[LOCK] KCBS quanten-max (Apex) = {kcbs(apex.astype(complex)):.4f} (Soll √5={np.sqrt(5):.4f})")
print(f"[LOCK] KCBS Komp.-Basis e0/e1/e2 = {[round(kcbs(np.eye(3,dtype=complex)[:,j]),3) for j in range(3)]} (generisch, <√5)")
# LOCK: Mana stabilizer/Komp.-Basis -> 0 ; Zufall > 0
print(f"[LOCK] Mana(d=3) Komp.-Basis = {[round(mana(np.eye(3,dtype=complex)[:,j]),3) for j in range(3)]} (=0) ; Zufall={mana(np.random.randn(3)+1j*np.random.randn(3)):.3f}(>0)")

print("\n[K-RESOLVED] max KCBS + max Mana ueber Braid-Orbit + Kopplung:")
print(f"  {'k':>3} {'Typ':>7} {'|grp|':>7} {'maxKCBS':>8} {'kontextuell?':>13} {'maxMana':>8} {'Mana@KCBSmax':>13} {'corr(K,Mana)':>13}")
rows=[]
for k in [2,4,8,3,5]:
    gens,B=build_braid(4,k,2)
    if len(B)!=3: print(f"  k={k}: dim={len(B)}!=3 uebersprungen"); continue
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
    ctx='JA (>2)' if maxK>2+1e-6 else 'nein'
    typ='endl.' if not dense else 'dicht'
    rows.append((k,typ,len(els) if not dense else None,maxK,ctx,maxM,mAtK,corr))
    print(f"  {k:>3} {typ:>7} {str(len(els) if not dense else '∞'):>7} {maxK:>8.4f} {ctx:>13} {maxM:>8.4f} {mAtK:>13.4f} {corr:>13.3f}")

print("\n  DEUTUNG: maxKCBS>2 = Braiding erreicht kontextuellen Qutrit-State; corr(K,Mana) = Kopplung der Achsen.")
print("  HEDGES: KEINE Erstheit fuer Braid-KCBS (Chou 2506.14537 state-dep., peer-rev.); Howard 1401.4174 = qudit/odd-prime-d (Analogie).")
import json
json.dump({'kcbs_classical':2,'kcbs_quantum':float(np.sqrt(5)),
           'rows':[{'k':k,'typ':t,'maxKCBS':float(mk),'contextual':c,'maxMana':float(mm),'Mana_at_KCBSmax':float(ma),'corr':float(co) if co==co else None} for (k,t,_,mk,c,mm,ma,co) in rows]},
          open('p5a_ergebnis_t5.json','w'),indent=2)
print("\nGeschrieben: p5a_ergebnis_t5.json")
