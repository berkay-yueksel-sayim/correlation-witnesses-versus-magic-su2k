"""
B_magic neighbour 1 — N1 HARDENING, 2026-06-23.
Two complementary magic witnesses that either HARDEN the soft truncation finding ("k=4 above
additive") or honestly expose it as soft. Built on the verified multi-strand engine
(multi_strand.py): exactly the same conventions (q=e^{i pi/(k+2)}, R symbols, sigma generators,
path basis B=sorted(paths(...))). NOTHING changed in the old engine.

(a) LEAKAGE-FREE 2-qubit encoding (exactly dim-4) — no renormalization artefact any more.
(b) TRUNCATION-FREE Mana(d=5) on the native fusion space (n=6, t=0, dim 5).

PRIOR ART (referenced in the text/AUSWERTUNG, never "first"):
- Universal braiding -> non-Clifford magic is a corollary of Freedman-Larsen-Wang (universality)
  + Gottesman-Knill. M_2 has already been tabulated explicitly for braiding magic in the D(S_3)
  quantum double (Byles et al., arXiv:2408.03377: M_2=log(16/13)). Our contribution: M_2/RoM of
  leakage-free 2-qubit SU(2)_k braid GATES as a function of k (worked example, no priority claim).
- Mana on SU(2)_k fusion/conformal-block spaces has been established since Fliss, JHEP 04 (2021) 090
  (arXiv:2011.01962) (odd prime d). Our contribution: Mana/Wigner negativity as a
  function of a concretely applied braid word on the prime-dim (d=5) fusion space.

Deterministic, seed 2026. Own derivation (no foreign code copied).
"""
import numpy as np
np.random.seed(2026)

# ======================================================================
# ENGINE (taken exactly from multi_strand.py — NO change of conventions)
# ======================================================================
def qint(m,k):  # quantum integer [m]_q, q=e^{i pi/(k+2)}
    return np.sin(m*np.pi/(k+2))/np.sin(np.pi/(k+2))

def braid_phases(k):
    d=2*np.cos(np.pi/(k+2)); h_half=0.5*1.5/(k+2); h1=2.0/(k+2)
    R0=-np.exp(1j*np.pi*(0.0-2*h_half))
    R1=np.exp(1j*np.pi*(h1-2*h_half))
    return d,R0,R1

def paths(n,k,total):
    res=[]
    def rec(p):
        i=len(p)
        if i==n+1:
            if p[-1]==total: res.append(tuple(p))
            return
        for step in (+1,-1):
            q=p[-1]+step
            if 0<=q<=k: rec(p+[q])
    rec([0])
    return res

def build_braid(n,k,total):
    d,R0,R1=braid_phases(k)
    B=sorted(paths(n,k,total)); idx={p:i for i,p in enumerate(B)}; D=len(B)
    gens=[]
    s1=np.zeros((D,D),dtype=complex)
    for p in B:
        s1[idx[p],idx[p]] = R0 if p[2]==0 else R1
    gens.append(s1)
    for m in range(2,n):
        em=np.zeros((D,D),dtype=complex)
        for p in B:
            if p[m-1]!=p[m+1]: continue
            base=p[m-1]
            opts=[v for v in (base-1,base+1) if 0<=v<=k]
            for a in opts:
                for b in opts:
                    pa=list(p); pa[m]=a; pb=list(p); pb[m]=b
                    if tuple(pa) in idx and tuple(pb) in idx:
                        em[idx[tuple(pa)],idx[tuple(pb)]] = np.sqrt(qint(a+1,k)*qint(b+1,k))/qint(base+1,k)
        sm = R1*np.eye(D) + ((R0-R1)/d)*em
        gens.append(sm)
    return gens,B

# sanity: all sigma_m unitary (a unitary action is mandatory for M2/Mana validity)
def assert_unitary(gens, tol=1e-9, tag=""):
    for i,g in enumerate(gens):
        err=np.max(np.abs(g.conj().T@g - np.eye(g.shape[0])))
        assert err<tol, f"sigma_{i+1} not unitary ({tag}): {err:.1e}"

# ======================================================================
# (a) LEAKAGE-FREE 2-QUBIT ENCODING (exactly dim-4)
# ======================================================================
# n=8, total=0 braid rep. Subspace S = {paths with p[4]==0}.
# Structure (verified): the 4 paths factorize into
#   qubit A := p[2] in {0,2}  (anyons 1-4 fuse to 0)
#   qubit B := p[6] in {0,2}  (anyons 5-8 fuse to 0)
# Komputational-Reihenfolge EXPLIZIT deklariert:
#   |00> p[2]=0,p[6]=0 ; |01> p[2]=0,p[6]=2 ; |10> p[2]=2,p[6]=0 ; |11> p[2]=2,p[6]=2
# (a=0 -> Bit 0 ; a=2 -> Bit 1)

I2=np.eye(2,dtype=complex); X=np.array([[0,1],[1,0]],dtype=complex)
Yp=np.array([[0,-1j],[1j,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
P1=[I2,X,Yp,Z]; P2q=[np.kron(a,b) for a in P1 for b in P1]   # 16 Zwei-Qubit-Paulis

def M2_2q(psi4):
    psi4=psi4/np.linalg.norm(psi4)
    ev=[np.real(np.vdot(psi4,P@psi4)) for P in P2q]
    return -np.log2(sum(e**4 for e in ev)/4.0)

# additive (local-block) single-qubit values for the baseline
STAB1=[np.array(v,dtype=complex)/np.linalg.norm(v) for v in ([1,0],[0,1],[1,1],[1,-1],[1,1j],[1,-1j])]
def M2_1q(psi):
    psi=psi/np.linalg.norm(psi); ev=[np.real(np.vdot(psi,P@psi)) for P in (I2,X,Yp,Z)]
    return -np.log2(sum(e**4 for e in ev)/2.0)

def setup_2q(k):
    gens,B=build_braid(8,k,0)            # gens[i] = sigma_{i+1}
    assert_unitary(gens, tag=f"n=8,k={k}")
    idx={p:i for i,p in enumerate(B)}; D=len(B)
    sub=[p for p in B if p[4]==0]
    assert len(sub)==4, f"k={k}: subspace p[4]==0 has dim {len(sub)} != 4"
    # Komputational-Basis in expliziter Reihenfolge
    def comp_path(a,b):  # a,b in {0,2}
        for p in sub:
            if p[2]==a and p[6]==b: return p
        raise RuntimeError("comp path fehlt")
    order=[(0,0),(0,2),(2,0),(2,2)]      # |00>,|01>,|10>,|11>
    comp_idx=[idx[comp_path(a,b)] for (a,b) in order]
    # projector P_S (onto the 4 subspace paths)
    PS=np.zeros((D,D),dtype=complex)
    for p in sub: PS[idx[p],idx[p]]=1.0
    return gens,B,idx,D,comp_idx,PS,sub,order

def leakage_norm(U,PS):
    # ||(I-P_S) U P_S|| : exact preservation of S  <=> 0
    D=PS.shape[0]; I=np.eye(D)
    return np.linalg.norm((I-PS)@U@PS, 2)

def U_on_comp(U,comp_idx):
    # 4x4 action on the computational subspace
    return U[np.ix_(comp_idx,comp_idx)]

def rand_word(gens,length):
    U=np.eye(gens[0].shape[0],dtype=complex)
    for _ in range(length): U=gens[np.random.randint(len(gens))]@U
    return U

# --- additive (local-block) baseline: words ONLY from {s1,s2,s3} x {s5,s6,s7} ---
# These preserve S exactly (block-diagonal in A x B). Max additive 2-qubit M2 = 2 x max single-qubit M2.
def additive_baseline(k, comp_idx, PS, gens, nwords=3000, length=24):
    # local generators: sigma_1,2,3 (group A) and sigma_5,6,7 (group B). gens[i]=sigma_{i+1}.
    locA=[gens[0],gens[1],gens[2]]            # s1,s2,s3
    locB=[gens[4],gens[5],gens[6]]            # s5,s6,s7
    loc = locA+locB
    # leakage check for the local generators: must be EXACTLY 0
    leakmax=max(leakage_norm(g,PS) for g in loc)
    # additive baseline = best possible 2-qubit M2 as a product of two single-qubit M2 (separable)
    # -> we measure the maximal single-qubit M2 over block-A words and over block-B, additively.
    # Block A acts on qubit A (p[2]), block B on qubit B (p[6]); on comp 4x4 they are kron factors.
    def block_single_m2(blockgens):
        best=0.0
        for _ in range(nwords):
            U=np.eye(gens[0].shape[0],dtype=complex)
            for _ in range(length): U=blockgens[np.random.randint(len(blockgens))]@U
            U4=U_on_comp(U,comp_idx)
            for s in STAB1:
                # start |s>_A (x) |0>_B  -> after block A, qubit B stays = |0>
                psiA=np.kron(s,np.array([1,0],dtype=complex))
                out=U4@psiA
                # reduce to qubit A (qubit B untouched) -> single-qubit M2 of qubit A
                # since separable: M2_1q of the qubit-A state
                a=out.reshape(2,2)  # index (A,B)
                # B should stay |0>; take column 0
                psA=a[:,0]
                if np.linalg.norm(psA)>1e-6: best=max(best,M2_1q(psA))
        return best
    mA=block_single_m2(locA); mB=block_single_m2(locB)
    return mA+mB, leakmax, mA, mB

# --- leakage-free ENTANGLERS: words from ALL generators (incl. s4) that preserve S exactly ---
def order_of_s4(k, gens, comp_idx, PS, maxord=64):
    # sigma_4 = gens[3]. Finite order at k=2,4,8. Check the powers.
    s4=gens[3]; D=s4.shape[0]; I=np.eye(D)
    U=I.copy(); res=[]
    for m in range(1,maxord+1):
        U=s4@U
        leak=leakage_norm(U,PS)
        res.append((m,leak,U.copy()))
        if np.max(np.abs(U-I))<1e-9:  # sigma_4^m = I
            res.append(('ord',m)); break
    return res

def op_schmidt_rank(U4, tol=1e-7):
    """Operator Schmidt rank of a 4x4 on 2 qubits. ==1 <=> a pure product gate kron(A,B)
       (=> NOT entangling). >1 <=> genuinely entangling."""
    T=U4.reshape(2,2,2,2).transpose(0,2,1,3).reshape(4,4)
    sv=np.linalg.svd(T,compute_uv=False)
    return int(np.sum(sv>tol)), sv

def product_error(U4):
    """||U4 - e^{i phi} kron(A,B)|| : 0 <=> exaktes Produktgatter."""
    T=U4.reshape(2,2,2,2).transpose(0,2,1,3).reshape(4,4)
    u,s,vh=np.linalg.svd(T)
    A=np.sqrt(s[0])*u[:,0].reshape(2,2); Bm=np.sqrt(s[0])*vh[0,:].reshape(2,2)
    prod=np.kron(A,Bm); ph=np.vdot(prod.flatten(),U4.flatten()); ph/=abs(ph) if abs(ph)>0 else 1
    return np.max(np.abs(U4-ph*prod)), A, Bm

def leakfree_entanglers(k, gens, comp_idx, PS, sub, order,
                        n_rand=20000, len_rand=30, structured=True):
    """Collect leakage-free words U (||(I-PS)U PS||<1e-9). KEY QUESTION: is a leakage-free
       word ever ENTANGLING (op Schmidt rank>1), or always a product gate kron(A,B)?
       Strategies: (i) long random + filter ; (ii) structured s4^m, commutators, conjugates.
       Returns: best_m2 (max 2-qubit M2 over comp. stabilizer inputs, leakage-free),
       best_leak, n_lf, n_entangling, max_prod_err, max_single (max single-qubit M2 per block)."""
    D=gens[0].shape[0]
    LEAK_TOL=1e-9
    best_m2=0.0; best_leak=None; n_found=0; n_ent=0; best_word=None
    max_prod_err=0.0; max_single=0.0
    comp_stab=[]  # 2-qubit stabilizer product inputs as comp. vectors
    for sa in STAB1:
        for sb in STAB1:
            comp_stab.append(np.kron(sa,sb))

    def eval_U(U,tag):
        nonlocal best_m2,best_leak,n_found,n_ent,best_word,max_prod_err,max_single
        leak=leakage_norm(U,PS)
        if leak>=LEAK_TOL: return False
        n_found+=1
        U4=U_on_comp(U,comp_idx)
        rank,_=op_schmidt_rank(U4)
        if rank>1: n_ent+=1
        else:
            perr,A,Bm=product_error(U4)
            max_prod_err=max(max_prod_err,perr)
            for sst in STAB1:
                if np.linalg.norm(A@sst)>1e-6:  max_single=max(max_single,M2_1q(A@sst))
                if np.linalg.norm(Bm@sst)>1e-6: max_single=max(max_single,M2_1q(Bm@sst))
        # max 2-qubit M2 over comp. stabilizer inputs (leakage-free: stays exactly in S)
        for psi4 in comp_stab:
            out=U4@psi4
            m2=M2_2q(out)
            if m2>best_m2:
                best_m2=m2; best_leak=leak; best_word=tag+(" [ENT]" if rank>1 else " [prod]")
        return True

    if structured:
        s4=gens[3]; I=np.eye(D,dtype=complex)
        U=I.copy()
        for m in range(1,65):
            U=s4@U; eval_U(U,f"s4^{m}")
            if np.max(np.abs(U-I))<1e-9: break
        gg=gens+[np.linalg.inv(g) for g in gens]
        for wlen in (1,2,3):
            for _ in range(2000):
                w=I.copy()
                for _ in range(wlen): w=gg[np.random.randint(len(gg))]@w
                wi=np.linalg.inv(w); m=np.random.randint(1,6)
                eval_U(w@np.linalg.matrix_power(s4,m)@wi,f"conj w(len{wlen}) s4^{m}")
        for _ in range(4000):
            a=gg[np.random.randint(len(gg))]; b=gg[np.random.randint(len(gg))]
            eval_U(a@b@np.linalg.inv(a)@np.linalg.inv(b),"commutator")

    for _ in range(n_rand):
        eval_U(rand_word(gens,len_rand),f"rand(len{len_rand})")

    if best_leak is None: best_leak=0.0
    return best_m2, best_leak, n_found, n_ent, max_prod_err, max_single, best_word

def run_part_a():
    print("="*72)
    print("(a) LEAKAGE-FREE 2-QUBIT ENCODING (exactly dim-4, n=8, S={p[4]==0})")
    print("    Komp.-Reihenfolge: |00>=(pA0,pB0) |01>=(pA0,pB2) |10>=(pA2,pB0) |11>=(pA2,pB2)")
    print("    add.base = 2 x max Single-Qubit-M2 (lokale Block-Woerter s1,2,3 x s5,6,7)")
    print("    n_ent = #leakage-free words with op Schmidt rank>1 (= genuinely entangling)")
    print("    lf_max = max 2-qubit M2 over leakage-free words (stabilizer inputs)")
    print(f"{'k':>2} {'add.base':>9} {'leak_loc':>9} {'lf_max':>8} {'n_lf':>6} {'n_ent':>6} {'prodErr':>9} {'2xsingle':>9}")
    rows=[]
    for k in [2,4,8,3,5]:
        gens,B,idx,D,comp_idx,PS,sub,order = setup_2q(k)
        addbase, leak_loc, mA, mB = additive_baseline(k, comp_idx, PS, gens)
        rnd = 9000 if k in (3,5) else 5000
        lf_max, lf_leak, n_lf, n_ent, prod_err, max_single, bw = leakfree_entanglers(
            k, gens, comp_idx, PS, sub, order, n_rand=rnd, len_rand=30, structured=True)
        two_single=2*max_single
        rows.append((k,addbase,leak_loc,lf_max,lf_leak,n_lf,n_ent,prod_err,two_single,max_single,bw))
        print(f"{k:>2} {addbase:>9.4f} {leak_loc:>9.1e} {lf_max:>8.4f} {n_lf:>6} {n_ent:>6} {prod_err:>9.1e} {two_single:>9.4f}")
        # LOCK k=2: 2-qubit M2 = 0 (Clifford), leakage exactly 0
        if k==2:
            assert addbase<1e-9, f"k=2 additive baseline != 0 (Clifford-Lock verletzt): {addbase:.2e}"
            assert lf_max<1e-9,  f"k=2 leakfree M2 != 0 (Clifford-Lock verletzt): {lf_max:.2e}"
            assert leak_loc<1e-9, f"k=2 lokale Leakage != 0: {leak_loc:.2e}"
            print("    LOCK k=2: additive=0, leakfree=0, lokale Leakage=0  -> Clifford-Konsistenz OK")
        else:
            # HARD LOCK: every leakage-free word is an EXACT product gate kron(A,B)
            assert n_ent==0, f"k={k}: {n_ent} verschraenkende leakage-freie Woerter gefunden (Befund waere HART)!"
            assert prod_err<1e-9, f"k={k}: leakage-free words not exactly a product (prodErr {prod_err:.1e})"
            # lf_max must not exceed 2xsingle (product bound)
            assert lf_max <= two_single+1e-6, f"k={k}: lf_max {lf_max:.4f} > 2xsingle {two_single:.4f}"
    print("    FINDING (a): n_ent==0 for all k -> EVERY leakage-free word is a product gate kron(A,B).")
    print("    => NO super-additive 2-qubit magic. lf_max <= 2 x single-qubit (= additive). (see AUSWERTUNG)")
    return rows

# ======================================================================
# (b) TRUNCATION-FREIE Mana(d=5) — Gross-2006 diskrete Wigner-Funktion, prime d=5
# ======================================================================
# The native fusion space n=6,t=0 has dim 5 (for k>=3). 5 is prime -> Wigner well defined.
# Heisenberg-Weyl d=5: X|j>=|j+1 mod 5>, Z|j>=w^j|j>, w=exp(2pi i/5).
# Displacement D(a,b)=w^{ -inv2 * a*b } Z^b X^a  (inv2 = (d+1)/2 = 3 for d=5)  [Gross 2006, prime d]
# Phasenpunkt-Operator A(0)=(1/d) sum_u D(u) ; A(u)=D(u) A(0) D(u)^dagger.
# W_psi(u)=<psi|A(u)|psi>/d (reell). Mana=log( sum_u |W_psi(u)| ).

d5=5
w5=np.exp(2j*np.pi/d5)
def build_HW(d=5):
    w=np.exp(2j*np.pi/d)
    X=np.zeros((d,d),dtype=complex)
    for j in range(d): X[(j+1)%d, j]=1.0
    Z=np.diag([w**j for j in range(d)])
    inv2=(d+1)//2  # multiplicative inverse of 2 mod d (d odd prime)
    # Displacement D(a,b)=w^{-inv2*a*b} Z^b X^a
    def D(a,b):
        return (w**(-inv2*a*b)) * np.linalg.matrix_power(Z,b) @ np.linalg.matrix_power(X,a)
    # A(0)=(1/d) sum_{a,b} D(a,b)
    A0=np.zeros((d,d),dtype=complex)
    for a in range(d):
        for b in range(d):
            A0+=D(a,b)
    A0/=d
    A={}
    for a in range(d):
        for b in range(d):
            Da=D(a,b); A[(a,b)]=Da@A0@Da.conj().T
    return X,Z,D,A0,A

Xd,Zd,Dd,A0d,Ad = build_HW(d5)

def wigner(psi):
    psi=psi/np.linalg.norm(psi)
    W={}
    for u,Au in Ad.items():
        W[u]=np.real(np.vdot(psi,Au@psi))/d5
    return W
def mana(psi):
    W=wigner(psi)
    return np.log2(sum(abs(v) for v in W.values()))  # log2 base (series convention, RUN-3 follow-up 2 2026-07-05; was np.log/ln — normalization Sigma_u A(u)=d*I identical to r9, verified)

def part_b_locks():
    print("="*72)
    print("(b) TRUNCATION-FREIE Mana(d=5) — Gross-2006 Wigner, prime d=5")
    # LOCK 0: A(u) hermitesch, sum_u A(u)=d*I, tr A(u)=1, A0^2 ~ projektor-artig (Wigner-Konsistenz)
    herm=max(np.max(np.abs(A-A.conj().T)) for A in Ad.values())
    sumA=np.max(np.abs(sum(Ad.values()) - d5*np.eye(d5)))
    trA=max(abs(np.trace(A)-1.0) for A in Ad.values())
    print(f"    HW-Locks: max|A-A^dag|={herm:.1e}  |sum A - dI|={sumA:.1e}  max|trA-1|={trA:.1e}")
    assert herm<1e-9 and sumA<1e-9 and trA<1e-9, "HW/Wigner-Konstruktion fehlerhaft"
    # LOCK (mandatory): every basis state |j> -> Wigner >=0 -> Mana=0
    maxbase=0.0
    for j in range(d5):
        e=np.zeros(d5,dtype=complex); e[j]=1.0
        m=mana(e); maxbase=max(maxbase,abs(m))
        W=wigner(e); negmin=min(W.values())
        assert m<1e-9, f"basis state |{j}> has Mana {m:.2e} != 0 (lock violated)"
        assert negmin>-1e-12, f"basis state |{j}> has negative Wigner {negmin:.2e}"
    print(f"    STABILIZER-INPUT-LOCK: max|Mana(|j>)|={maxbase:.1e} (all 5 basis states, Wigner>=0) -> 0 OK")
    # generischer Zufallszustand -> Mana > 0
    np.random.seed(2026)
    v=np.random.randn(d5)+1j*np.random.randn(d5); mg=mana(v)
    print(f"    Kontroll: generischer Zufallszustand Mana={mg:.4f} (>0 erwartet)  {'OK' if mg>1e-3 else 'FAIL'}")
    assert mg>1e-3, "generic state not magical -> Mana implementation suspicious"
    return maxbase, mg

def run_part_b():
    maxbase, mg = part_b_locks()
    print("    APPLICATION: max Mana over the braid orbit of the dim-5 states (encoding-free, HW-structure-dep.)")
    print("    NOTE: k=2 has dim(n=6,t=0)=4 != 5 -> d=5 Wigner NOT natively applicable (honestly skipped).")
    print(f"    NOTE: the k=2->Mana-0 lock does NOT apply here (it was a d=2 statement). The k=2 value is NOT enforced.")
    print(f"{'k':>2} {'dim':>4} {'maxMana':>9} {'#orbit':>7}  {'note':>10}")
    rows=[]
    for k in [2,3,4,5,8]:
        gens,B=build_braid(6,k,0); D=len(B)
        assert_unitary(gens, tag=f"n=6,k={k}")
        if D!=d5:
            print(f"{k:>2} {D:>4} {'   n/a':>9} {0:>7}  {'dim!=5':>10}")
            rows.append((k,D,None,0,"dim!=5"))
            continue
        # Orbit: starte in jedem Basiszustand, wende Braid-Woerter an
        best=0.0; cnt=0
        if k in (2,4,8):
            # finite/structured: enumerate short words via BFS up to cap
            from itertools import product
            gg=gens+[np.linalg.inv(g) for g in gens]
            seen=set(); frontier=[np.eye(D,dtype=complex)]
            # BFS over the group (canonicalize via the rounded matrix)
            def key(U): return tuple(np.round(U.flatten(),5).view(float))
            seenU={key(np.eye(D,dtype=complex)):np.eye(D,dtype=complex)}
            fr=[np.eye(D,dtype=complex)]; CAP=2000
            while fr and len(seenU)<CAP:
                U=fr.pop()
                for g in gg:
                    V=g@U; ky=key(V)
                    if ky not in seenU:
                        seenU[ky]=V; fr.append(V)
            els=list(seenU.values()); cnt=len(els)
            for U in els:
                for j in range(D):
                    e=np.zeros(D,dtype=complex); e[j]=1.0
                    best=max(best,mana(U@e))
        else:
            NW=6000; L=24
            for _ in range(NW):
                U=np.eye(D,dtype=complex)
                for _ in range(L): U=gens[np.random.randint(len(gens))]@U
                for j in range(D):
                    e=np.zeros(D,dtype=complex); e[j]=1.0
                    best=max(best,mana(U@e))
                cnt+=1
        rows.append((k,D,best,cnt,"orbit"))
        print(f"{k:>2} {D:>4} {best:>9.4f} {cnt:>7}  {'orbit':>10}")
    return rows, maxbase, mg

# ======================================================================
if __name__=="__main__":
    print("#"*72)
    print("# N1 HARDENING — seed 2026 — two magic witnesses")
    print("#"*72)
    rows_a = run_part_a()
    rows_b, mana_lock, mana_rand = run_part_b()

    print("="*72)
    print("ZUSAMMENFASSUNG (echte Zahlen)")
    print("-- (a) leakage-frei 2-Qubit M2 (n_ent = #verschraenkende leakage-freie Woerter) --")
    print(f"{'k':>2} {'2xsingle':>9} {'lf_max':>8} {'n_ent':>6} {'prodErr':>9}")
    for (k,addbase,leak_loc,lf_max,lf_leak,n_lf,n_ent,prod_err,two_single,max_single,bw) in rows_a:
        print(f"{k:>2} {two_single:>9.4f} {lf_max:>8.4f} {n_ent:>6} {prod_err:>9.1e}")
    print("-- (b) Mana(d=5) --")
    print(f"{'k':>2} {'dim':>4} {'maxMana':>9}")
    for (k,D,best,cnt,note) in rows_b:
        bs = f"{best:.4f}" if best is not None else "n/a"
        print(f"{k:>2} {D:>4} {bs:>9}  ({note})")
    print(f"Mana stabilizer-input-Lock: max|Mana(|j>)|={mana_lock:.1e} (=0) ; random={mana_rand:.4f}")

    # provenance deposit for (b) Mana(d=5), log2 base (RUN-3 follow-up 2 2026-07-05).
    # Part (a) has no JSON (not part of this release, untouched).
    import json
    from pathlib import Path
    json.dump({'basis':'log2','stabilizer_input_lock_maxabs':float(mana_lock),'random_control':float(mana_rand),
               'rows_b':[{'k':k,'dim':D,'maxMana':(float(best) if best is not None else None),'orbit':cnt,'note':note}
                         for (k,D,best,cnt,note) in rows_b]},
              open(Path(__file__).resolve().parent / 'p5a_ergebnis_hardening_b.json','w'),indent=2)
    print("Written: p5a_ergebnis_hardening_b.json")
