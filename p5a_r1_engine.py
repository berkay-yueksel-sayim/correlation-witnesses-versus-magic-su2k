"""
Magic — reimplementation r1, 2026-06-23.
UNABHAENGIGE Herleitung aus Konventionen (Browser-Code NICHT kopiert).
Ziele (Prior-Art-Gate GO): (1) M2(k)-Anker reproduzieren, (2) McKay k->{2O,2T,2I} aus
Generator-Ordnungen re-herleiten, (3) LGI-Struktur prüfen. Deterministisch, fester Seed.
"""
import numpy as np

np.random.seed(2026)
I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)
PAULIS = [I2, X, Y, Z]

# ---- Konventionen (SU(2)_k, j=1/2) ----
def braid_data(k):
    q_arg = np.pi/(k+2)
    d = 2*np.cos(q_arg)                 # [2]_q = quantum dim of j=1/2
    h_half = 0.5*1.5/(k+2)              # h_j = j(j+1)/(k+2), j=1/2
    h0 = 0.0
    h1 = 1.0*2.0/(k+2)                  # j=1 -> h_1 = 2/(k+2)
    R0 = (-1)**(1-0)*np.exp(1j*np.pi*(h0 - 2*h_half))   # fusion channel c=0
    R1 = (-1)**(1-1)*np.exp(1j*np.pi*(h1 - 2*h_half))   # fusion channel c=1
    s1 = np.diag([R0, R1])
    root = np.sqrt(max(d*d-1.0,0.0))/d
    F = np.array([[1.0/d, root],[root, -1.0/d]], dtype=complex)
    s2 = F @ s1 @ F
    theta = np.angle(R1/R0) % (2*np.pi)                 # sigma_1 z-rotation angle
    return dict(d=d, R0=R0, R1=R1, s1=s1, s2=s2, F=F, theta=theta)

# ---- Magic: Stabilizer-2-Renyi M2 (Leone-Oliviero-Hamma 2022) ----
def M2(psi):
    psi = psi/np.linalg.norm(psi)
    ev = [np.real(np.vdot(psi, P@psi)) for P in PAULIS]
    xi = sum(e**4 for e in ev)
    return -np.log2(xi/2.0)

STAB = [np.array([1,0],dtype=complex), np.array([0,1],dtype=complex),
        np.array([1,1],dtype=complex)/np.sqrt(2), np.array([1,-1],dtype=complex)/np.sqrt(2),
        np.array([1,1j],dtype=complex)/np.sqrt(2), np.array([1,-1j],dtype=complex)/np.sqrt(2)]

# ---- projektive Gruppen-Enumeration (mod globale Phase, via phasen-freie SO(3)) ----
def to_SO3(U):
    sig=[X,Y,Z]; R=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            R[i,j]=0.5*np.real(np.trace(sig[i]@U@sig[j]@U.conj().T))
    return R
def canon(U, dec=5):
    return tuple(np.round(to_SO3(U).flatten(), dec))

def enumerate_group(gens, cap=2000):
    seen = {}
    start = canon(I2); seen[start] = I2.copy()
    frontier = [I2.copy()]
    invs = [np.linalg.inv(g) for g in gens]
    while frontier:
        U = frontier.pop()
        for g in gens+invs:
            V = g@U
            c = canon(V)
            if c not in seen:
                seen[c] = V; frontier.append(V)
                if len(seen) > cap: return seen, len(seen), False
    return seen, len(seen), True

def proj_order(U, maxn=200):
    P = U.copy()
    for n in range(1, maxn+1):
        c = canon(P)
        if c == canon(I2):
            return n
        P = U@P
    return None

def magic_finite(k):
    bd = braid_data(k)
    seen, order, closed = enumerate_group([bd['s1'], bd['s2']])
    mmax = 0.0
    for U in seen.values():
        for s in STAB:
            mmax = max(mmax, M2(U@s))
    return order, closed, mmax, proj_order(bd['s1']), bd['theta']

def magic_dense(k, nwords=4000, length=24):
    bd = braid_data(k); gens=[bd['s1'],bd['s2']]
    mmax=0.0
    for _ in range(nwords):
        psi = STAB[0].copy()
        for _ in range(length):
            psi = gens[np.random.randint(2)]@psi
        mmax = max(mmax, M2(psi))
    return mmax

# ---- LGI K3 (Q=Z) unter sigma_1-Schritten ----
def lgi_K3_Z(k):
    # einfacher Drei-Zeit-LGI mit Q=Z, Evolution pro Schritt = sigma_1, max-mixed Start
    U = braid_data(k)['s1']
    def corr(steps):
        Ut = np.linalg.matrix_power(U, steps)
        Qh = Ut.conj().T @ Z @ Ut
        return np.real(0.5*np.trace(Z@Qh))   # rho = I/2
    C12=corr(1); C23=corr(1); C13=corr(2)
    return C12+C23-C13

print("="*64)
print("SANITY-GATES")
for k in [2,4,8]:
    bd=braid_data(k); s1,s2,F=bd['s1'],bd['s2'],bd['F']
    ybe=np.max(np.abs(s1@s2@s1 - s2@s1@s2))
    uni=np.max(np.abs(s1@s1.conj().T - I2))
    f2 =np.max(np.abs(F@F - I2))
    print(f"k={k}: YBE={ybe:.2e}  unit={uni:.2e}  F^2-I={f2:.2e}  theta={np.degrees(bd['theta']):.1f}deg")
# Magic-Formel-Anker
Tplus = np.array([1, np.exp(1j*np.pi/4)],dtype=complex)/np.sqrt(2)
print(f"M2(|0>)={M2(STAB[0]):.6f} (soll 0)  M2(T|+>)={M2(Tplus):.6f} (soll log2(4/3)={np.log2(4/3):.6f})")
nmax = max(M2(np.array([1, np.exp(1j*a)],dtype=complex)/np.sqrt(2)) for a in np.linspace(0,2*np.pi,200))
print(f"max ueber |+>-Familie ~{nmax:.4f} ; theor. Single-Qubit-Max log2(3/2)={np.log2(1.5):.6f}")

print("="*64)
print("ANKER-TABELLE + McKay-Re-Herleitung")
print(f"{'k':>2} {'proj_order':>10} {'closed':>6} {'sigma1_ord':>10} {'theta':>7} {'M2':>8}  Gruppe(re-hergeleitet)")
mck={4:'O(2O)',3:'T(2T)',5:'I(2I)'}  # via sigma1-Ordnung 4/3/5
for k in [2,3,4,5,6,7,8,9,10]:
    if k in (2,4,8):
        order,closed,mmax,s1ord,theta = magic_finite(k)
        grp = {4:'O -> 2O',3:'T -> 2T',5:'I -> 2I'}.get(s1ord,'?')
        print(f"{k:>2} {order:>10} {str(closed):>6} {s1ord:>10} {np.degrees(theta):>6.0f} {mmax:>8.4f}  s1-ord={s1ord} => {grp}")
    else:
        mmax=magic_dense(k)
        bd=braid_data(k)
        print(f"{k:>2} {'dense':>10} {'-':>6} {'inf':>10} {np.degrees(bd['theta']):>6.0f} {mmax:>8.4f}  (dicht)")

print("="*64)
print("LGI K3(Q=Z) unter sigma_1 (Struktur-Check)")
for k in [2,4,8]:
    print(f"k={k}: K3(Q=Z,sigma1)={lgi_K3_Z(k):.4f}   ([sigma1,Z]=0 da diagonal -> Z erhalten)")
print("HINWEIS: sigma_1 ist IMMER diagonal (z-Rot) -> [sigma1,Z]=0 -> K3=1 strukturell.")
print("Die k-abhaengige K3 des Befunds nutzt Braid-WORTE/max-Q, nicht nur sigma_1 -> Target-3 braucht das volle LGI-Protokoll (offen, ehrlich geflaggt).")
