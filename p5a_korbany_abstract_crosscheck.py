# -*- coding: utf-8 -*-
# Manager cross-check: Korbany MI witness for doubled Fibonacci (independent, own numpy)
import numpy as np

def H2(p):
    p = np.asarray(p, float); p = p[p > 1e-15]
    return float(-(p*np.log2(p)).sum())

phi = (1+np.sqrt(5))/2
# doubled Fibonacci: sectors (a,b) in {1,tau}x{1,tau}, d_(a,b)=d_a*d_b
d = np.array([1.0, phi, phi, phi**2])         # (1,1),(1,tau),(tau,1),(tau,tau)
D2 = (d**2).sum()
p = d**2 / D2
H_dfib = H2(p)
# bonus closed form: H = log2(D^2) - sum p*log2(d^2)
H_closed = np.log2(D2) - (p*np.log2(d**2)).sum()

# Anchor 1: Z2 toric (4 sectors d=1) -> stabilizer -> integer -> no certificate (Thm 2)
H_toric = H2([1,1,1,1]/np.array(4.0))
# Anchor 2: toric T-state, binary entropy {cos^2(pi/8), sin^2(pi/8)}
c = np.cos(np.pi/8)**2; s = np.sin(np.pi/8)**2
H_Tstate = H2([c, s])

print("=== doubled-Fibonacci ===")
print("p (sector weights d^2/D^2):", np.round(p,4), " (expected 0.0764,0.2,0.2,0.5236)")
print(f"D^2 = {D2:.4f}  (expected 13.09 = (1+phi^2)^2)")
print(f"H_dFib (direct)      = {H_dfib:.6f}   (expert: 1.700979)")
print(f"H_dFib (closed-form) = {H_closed:.6f}   |diff| = {abs(H_dfib-H_closed):.2e}")
print(f"noninteger? -> LRN zertifiziert: {abs(H_dfib-round(H_dfib))>1e-6}")
print("\n=== Anchors ===")
print(f"Z2-Toric  H = {H_toric:.6f}   (expected 2.0 integer -> NO certificate)")
print(f"Toric-T   H = {H_Tstate:.6f}   (expert says 0.6009; slip value 0.5436 -> ?)")
print(f"  -> 0.601 correct? {abs(H_Tstate-0.601)<2e-3} | 0.5436 wrong? {abs(H_Tstate-0.5436)>1e-2}")
