# -*- coding: utf-8 -*-
# Manager-Gegenrechnung: Korbany-MI-Witness fuer doubled-Fibonacci (unabhaengig, eigener numpy)
import numpy as np

def H2(p):
    p = np.asarray(p, float); p = p[p > 1e-15]
    return float(-(p*np.log2(p)).sum())

phi = (1+np.sqrt(5))/2
# doubled-Fibonacci: Sektoren (a,b) in {1,tau}x{1,tau}, d_(a,b)=d_a*d_b
d = np.array([1.0, phi, phi, phi**2])         # (1,1),(1,tau),(tau,1),(tau,tau)
D2 = (d**2).sum()
p = d**2 / D2
H_dfib = H2(p)
# Bonus-Closed-Form: H = log2(D^2) - sum p*log2(d^2)
H_closed = np.log2(D2) - (p*np.log2(d**2)).sum()

# Anker 1: Z2-Toric (4 Sektoren d=1) -> Stabilizer -> integer -> kein Zertifikat (Thm 2)
H_toric = H2([1,1,1,1]/np.array(4.0))
# Anker 2: Toric-T-State, binaere Entropie {cos^2(pi/8), sin^2(pi/8)}
c = np.cos(np.pi/8)**2; s = np.sin(np.pi/8)**2
H_Tstate = H2([c, s])

print("=== doubled-Fibonacci ===")
print("p (Sektor-Gewichte d^2/D^2):", np.round(p,4), " (erwartet 0.0764,0.2,0.2,0.5236)")
print(f"D^2 = {D2:.4f}  (erwartet 13.09 = (1+phi^2)^2)")
print(f"H_dFib (direkt)      = {H_dfib:.6f}   (Experte: 1.700979)")
print(f"H_dFib (closed-form) = {H_closed:.6f}   |diff| = {abs(H_dfib-H_closed):.2e}")
print(f"noninteger? -> LRN zertifiziert: {abs(H_dfib-round(H_dfib))>1e-6}")
print("\n=== Anker ===")
print(f"Z2-Toric  H = {H_toric:.6f}   (erwartet 2.0 integer -> KEIN Zertifikat)")
print(f"Toric-T   H = {H_Tstate:.6f}   (Experte sagt 0.6009; Slip-Wert 0.5436 -> ?)")
print(f"  -> 0.601 korrekt? {abs(H_Tstate-0.601)<2e-3} | 0.5436 falsch? {abs(H_Tstate-0.5436)>1e-2}")
