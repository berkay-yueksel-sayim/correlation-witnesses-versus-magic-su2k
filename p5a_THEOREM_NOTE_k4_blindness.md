# Beweis-Note: Blindheit des dichotomen Lüders-LGI-Zeugen unter endlichen Braid-Bildern — Alignment-Geometrie (Orbit-Gram), mit Pauli-aligned-Clifford-Spezialfall

**Lauf:** 2026-07-07_B_k4blind_r1 (Auftrag B2), Reparatur 2026-07-08 · **Status:** zwei Theoreme, jedes Lemma numerisch spot-gecheckt.
**Theorem A** (Pauli-aligned-Clifford-Spezialfall): `spot_checks.py` (seed=42, deterministisch). **Theorem B** (Orbit-Gram-Geometrie, deckt das reale k=4): `spot_checks_gram.py` (deterministisch, Referenz A-p07-Reimpl reproduziert). Rohausgaben beider am Ende dieser Note.
**Wording:** Sagbar ist „dieser Zeuge ist strukturell blind, weil Q-Achse und Orbit ausgerichtet sind" — nicht mehr, keine Erstheit. Fremde Resultate sind zugeschrieben (Jones; Tuba–Wenzl; Freedman–Larsen–Wang; Budroni–Emary).

---

## Setup und Definitionen

- Qubit, festes dichotomes Q = n₀·σ (Standardrahmen: Q = Pauli-Z, n₀ = e_z). Drei Messzeiten t1 < t2 < t3; Propagation zwischen den Zeiten durch Elemente einer Gruppe G ⊂ SU(2): V1 (t1→t2), V2 (t2→t3). Heisenberg-Observablen:
  Q1 = Q, Q2 = V1† Q V1, Q3 = (V2V1)† Q (V2V1). Deren Bloch-Richtungen n1 = n₀, n2, n3 liegen im Orbit von n₀ unter dem SO(3)-Bild von G; **n2 und n3 sind unabhängig frei wählbar, weil V1 und V2 unabhängig sind.**
- **Lüders-Korrelator** für dichotome projektive Messungen (Standard-Konvention, vgl. Fritz 2010; Budroni–Emary arXiv:1309.3678):
  C_ij = Σ_{a,b=±1} a·b · P(a bei t_j, b bei t_i) = (1/2) tr(ρ {Q_i, Q_j}).
- **Zeuge:** K3 = C21 + C32 − C31; makrorealistische Schranke K3 ≤ 1.

---

# THEOREM A (Pauli-aligned-Clifford-Spezialfall)

Sei **G eine Untergruppe der Single-Qubit-Clifford-Gruppe im Rahmen von Q** — d. h. G normalisiert die von Q aufgespannte Pauli-Gruppe, äquivalent: Q liegt im signierten Pauli-Achsen-Raster ±{X, Y, Z} und das SO(3)-Bild von G permutiert dieses Raster. Dann gilt für den dichotomen Lüders-LGI-Zeugen mit festem Q = Pauli-Z und Propagatoren aus G:

**K3 ≤ 1 für ALLE ρ.** Der Zeuge ist strukturell blind, zustandsunabhängig.

> **⚠ WARNHINWEIS (Kern der Reparatur 2026-07-08):** Das **reale** k=4-Braid-Bild von SU(2)_k=4 erfüllt diese Voraussetzung **im Fusions-Z-Rahmen NICHT.** Das SO(3)-Bild von σ₁ ist dort eine projektive **3-fach-Drehung um z** (cos θ = −1/2), also ist z eine 3-fach-Achse der tetraedrischen Gruppe T, keine Pauli-2-fach-Achse; Orbit(Z) besteht aus **4 Tetraeder-Richtungen** (z + 3 Richtungen mit z-Komponente −1/3, paarweise Bloch-Dots exakt −1/3), NICHT aus ±{X, Y, Z}. Insbesondere σ₂† Z σ₂ = Bloch(−√2/3, −√6/3, −1/3) ist **kein** ±Pauli. Die k=4-Gruppe ist zwar zu einer Clifford-Untergruppe konjugiert, aber die Konjugation fixiert Z nicht — Lemma 1 bricht. **Die k=4-Blindheit gilt trotzdem, aber über die Orbit-Gram-Geometrie (THEOREM B), nicht über THEOREM A.** (Erstfassung dieser Note behauptete fälschlich k=4 ⊂ Clifford im Z-Rahmen; siehe Ehrlichkeits-Block.)

## Lemma A1 (Clifford-Konjugation hält Z im signierten Pauli-Orbit)

**Behauptung:** U Clifford (im Q-Rahmen) ⇒ U†ZU = ±P für ein P ∈ {X, Y, Z}.

**Beweis:** Die Clifford-Gruppe ist per Definition der Normalisator der Pauli-Gruppe: U†ZU liegt in der Pauli-Gruppe, also U†ZU = ω·P' mit ω Phase, P' ∈ {I, X, Y, Z}. Hermitezität von U†ZU erzwingt ω = ±1; Spurtreue der Konjugation (tr(U†ZU) = tr Z = 0) schließt P' = I aus. Also U†ZU = ±P, P ∈ {X, Y, Z}. Da G Untergruppe ist, sind V1 und V2V1 beide in G ⊂ Clifford ⇒ **alle drei Heisenberg-Observablen Z1, Z2, Z3 sind signierte Paulis.** ∎

**Spot-Check 1:** Alle 24 Elemente von 2T (8 Lipschitz- + 16 Hurwitz-Einheitsquaternionen → SU(2)) explizit aufgebaut, Unitarität + det=1 + Gruppenabschluss (24×24 Produkte) verifiziert; U†ZU ∈ {±X, ±Y, ±Z} für **24/24** Elemente, Verteilung gleichmäßig 4 pro Ziel. **Hinweis:** 2T tritt hier als **abstrakte** Pauli-aligned-Clifford-Untergruppe auf (Achsen im Pauli-Raster) — das ist eine ANDERE Einbettung als das reale k=4-Fusions-Bild (siehe Warnhinweis / THEOREM B). Siehe Rohausgabe CHECK 1.

## Lemma A2 (Lüders-Korrelatoren sind diskret und ρ-unabhängig)

**Behauptung:** Für signierte Paulis A, B gilt (1/2){A, B} = ±I falls B = ±A, sonst 0. Folglich ist C_ij = (1/2) tr(ρ {Z_i, Z_j}) entweder ±1 (**ρ-unabhängig**) oder 0.

**Beweis:** In d=2 gibt es für zwei signierte Paulis nur zwei Fälle: (i) gleiche Pauli-Achse, B = ±A ⇒ (1/2){A, B} = ±A² = ±I; (ii) verschiedene Achsen ⇒ die Paulis antikommutieren, {A, B} = 0. Einsetzen: C_ij = ±tr(ρ) = ±1 bzw. C_ij = 0 — ρ fällt vollständig heraus. ∎

**Spot-Check 2:** Alle 36 signierten Pauli-Paare: 12× (1/2){A,B} = ±I (gleiche Achse), 24× = 0 (verschiedene Achse), Toleranz 1e-12. Korrelator über 206 Zufallszustände: maximale Spannweite 4.44e-16 (ρ-unabhängig), alle Werte in {−1, 0, +1}. Siehe Rohausgabe CHECK 2.

## Lemma A3 (Vorzeichen-Konsistenz: alle Fälle enden bei K3 ≤ 1)

Schreibe Z2 = s2·P2, Z3 = s3·P3 (Lemma A1). Aus Lemma A2 folgt exakt:
C21 = s2·δ(P2=Z), C32 = s2·s3·δ(P3=P2), C31 = s3·δ(P3=Z). Vollständige Fall-Liste (36 Kombinationen):

| Fall | Bedingung | C21 | C32 | C31 | K3 = C21+C32−C31 | # |
|---|---|---|---|---|---|---|
| A | P2 = Z, P3 = Z | e1 := s2 | e1e2 (e2 := s3s2) | e1e2·e1 = s3 | e1 + e2 − e1e2 ∈ {1, −3} | 4 |
| B | P2 = Z, P3 ≠ Z | ±1 | 0 | 0 | ±1 | 8 |
| C | P2 ≠ Z, P3 = P2 | 0 | ±1 | 0 (da P3 = P2 ≠ Z) | ±1 | 8 |
| D | P2 ≠ Z, P3 = Z | 0 | 0 | ±1 | ∓1 | 8 |
| E | P2 ≠ Z, P3 ∉ {Z, P2} | 0 | 0 | 0 | 0 | 8 |

**Fall A ausgeschrieben** (der einzige mit drei nichttrivialen Korrelatoren): Z2 = e1·Z, Z3 = e2·Z2 ⇒ Z3 = e1e2·Z ⇒ C31 = C21·C32, also K3 = e1 + e2 − e1e2. Werte: (+,+) → 1, (+,−) → 1, (−,+) → 1, (−,−) → −3. Maximum 1.
**Fälle B–E:** Sobald mindestens ein Korrelator 0 ist, propagieren die (Anti-)Kommutations-Relationen und höchstens ein Term von K3 ist ±1. **Keine der 36 Kombinationen erreicht K3 > 1.** ∎

**Spot-Check 3:** Alle 36 Kombinationen: Fall-Formel gegen Brute-Force über 56 Zufalls-ρ, Abweichung < 1e-10; K3-Wertemengen pro Fall exakt wie in der Tabelle; max K3 = 1. Siehe Rohausgabe CHECK 3.

## Beweis von Theorem A

Lemma A1 (mit Untergruppen-Abschluss: V1, V2V1 ∈ G) ⇒ Z1, Z2, Z3 signierte Paulis. Lemma A2 ⇒ jeder der drei Lüders-Korrelatoren liegt ρ-unabhängig in {−1, 0, +1}. Lemma A3 ⇒ jede der 36 Konfigurationen ergibt K3 ≤ 1. Da die Korrelatoren ρ-unabhängig sind, gilt die Schranke für **alle** ρ. ∎

**Spot-Check 4 (end-zu-end):** Alle 576 Propagator-Paare (V1, V2) ∈ 2T × 2T, je 14 Zufalls-ρ: Z2, Z3 stets signierte Paulis; **max K3 = 1.000000000000** — Schranke exakt erreicht (z. B. V1 = V2 = I), nie überschritten. Siehe Rohausgabe CHECK 4.

**Korollar A-a (k=2, abstrakt Pauli-aligned):** Das Braid-Bild bei k=2 ist binär-oktaedrisch 2O; **als abstrakte Gruppe** ist 2O eine Pauli-aligned-Clifford-Untergruppe (SO(3)-Bild = volle Oktaedergruppe, permutiert ±{X,Y,Z}) — Theorem A gilt wörtlich. **Spot-Check 1b:** alle 48 Elemente von 2O konjugieren Z auf ±Pauli (48/48), Gruppenabschluss 48×48 verifiziert. (Ob das reale k=2-Fusions-Bild diese Achsen-Ausrichtung trägt, ist eine Konventionsfrage; THEOREM B G2 rechnet den k=2-Orbit direkt und findet ebenfalls K3max = 1.)

**Spot-Check 5 (Kontrast, Nicht-Clifford):** Nicht-Clifford-Präzession U = exp(−iπX/6) verletzt Lemma A1 (U†ZU kein ±Pauli) und liefert K3 = 1.500000000000 für alle getesteten ρ = Lüders-Schranke 3/2 (Budroni–Emary). Siehe Rohausgabe CHECK 5.

---

# THEOREM B (Orbit-Gram-Geometrie — deckt das reale k=4)

Sei G ⊂ SU(2) beliebig und Q = n₀·σ fest. Dann sind die Lüders-Korrelatoren **Bloch-Skalarprodukte der Heisenberg-Richtungen und ρ-unabhängig:**

**C_ij = n_i · n_j für JEDES ρ**, wobei n_i ∈ Orbit(n₀) unter dem SO(3)-Bild von G, mit n1 = n₀ fest und n2, n3 unabhängig frei (V1, V2 unabhängig).

**Beweis:** Für Single-Qubit-Dichotome ist {n·σ, m·σ} = 2(n·m)·I (Standard-Pauli-Algebra), also C_ij = (1/2) tr(ρ · 2(n_i·n_j) I) = (n_i·n_j)·tr(ρ) = n_i·n_j — ρ fällt exakt heraus (kein Clifford, kein Pauli-Raster nötig). ∎

**Folge:** K3max = max über Orbit-Paare (n2, n3) von (n₀·n2 + n2·n3 − n₀·n3). Blindheit (K3max = 1) ⇔ die Q-Achse n₀ ist mit einer Symmetrieachse des Orbits **ausgerichtet**; Endlichkeit oder Clifford-Eigenschaft allein genügen NICHT (siehe Anwendung iv).

**Anwendungen (exakte Werte, `spot_checks_gram.py`; Referenz A-p07-Reimpl reproduziert):**

| # | System | Orbit | Bloch-Dots (offdiag) | K3max | Deutung |
|---|---|---|---|---|---|
| (i) | **k=4 / Fusions-Z** | 4 Tetraeder-Richtungen | {−1/3} | **1 exakt** | **blind** (Alignment: z = 3-fach-Achse) |
| (ii) | k=2 (O, z = 4-fach-Achse) | 6 (±{x,y,z}) | {0, ±1} | 1 exakt | blind |
| (iii) | k=8 (I) | 12 Ikosaeder-Ecken | {±1, ±1/√5} | 3/√5 ≈ 1.3416 | **sieht** (> 1) |
| (iv) | **Alignment-Korollar:** Oktaeder-Gruppe, gedrehtes Q=(1,1,0)/√2 | 12 Kantenmitten | {0, ±½, ±1} | **3/2 exakt** | **nicht blind** |

**Fallliste zu (i)** (Tetraeder-Gram, Dot = 1 bei gleichem Index, sonst −1/3): über alle 4³ = 64 Orbit-Tripel ist max(d12 + d23 − d13) = 1 exakt (erreicht z. B. n1 = n2, n3 ≠ n1: 1 − 1/3 + 1/3 = 1). **Spot-Check G1/G2** (unten): reales k=4-Bild → 4 Tetraeder-Richtungen, alle Paar-Dots −1/3, σ₂†Zσ₂ kein ±Pauli; K3max = 1.

**Anwendung (iv) — Kern der zweiten Reparatur:** Dieselbe abstrakte Gruppe (chirale Oktaedergruppe O, 24 Elemente; als 2O sogar **volles Clifford-Bild**), aber mit um 45° **gedrehter** Q-Achse (1,1,0)/√2 liefert K3max = **3/2 exakt** — die Lüders-Schranke, weit über 1. **„Endliche Gruppe ⇒ blind" ist als Slogan FALSCH.** Blindheit ist Eigenschaft des **ausgerichteten Paars (Gruppe, Q-Achse)**, nicht der Endlichkeit und nicht der Clifford-Eigenschaft per se. **Spot-Check G4** (Gram) + **G5** (operationale Quanten-Gegenprobe mit zwei echten SU(2)-Qubit-Unitaries): K3 = 3/2, ρ-unabhängig, == Gram-Formel.

**Scope-Hinweis zu (iii):** Im **festen Z-Rahmen** sieht der Zeuge bei k=8 K3max = 3/√5 ≈ 1.3416. Der Referenzwert K3opt = 1.427 = 2cos72° − cos144° stammt aus dem **max-Q-Protokoll** (Q-Richtung ⊥ zur 5-fach-Achse optimiert); im festen Z-Rahmen liegt keine 5-fach-Achse ⊥ z, daher dort unerreichbar. Kein Widerspruch, sondern Protokoll-Scope-Unterschied — konsistent mit Anwendung (iv): auch bei k=8 hängt K3max an der Q-Achsen-Ausrichtung.

---

## Korollare (übergreifend)

**(a) Dichte Braid-Bilder (Kontrast).** Für k=3 und k≥5 ist das Braid-Bild dicht in der projektiven Unitärgruppe (Freedman–Larsen–Wang-Universalität; Endlichkeits-/Dichtheits-Klassifikation Jones 1986, Tuba–Wenzl math/9912013 — fremde Resultate). Dann ist der Orbit von Q dicht auf der Bloch-Sphäre, C_ij wird kontinuierlich, und K3 kann bis zur Lüders-Schranke **3/2** steigen (Budroni–Emary, arXiv:1309.3678). **Spot-Check 5:** U = exp(−iπX/6) liefert K3 = 3/2.

**(b) Treiber ist das ALIGNMENT, nicht Clifford/Endlichkeit (Schärfung).** Der Zeuge ist genau dann blind, wenn die feste Q-Achse mit einer hinreichend hohen Symmetrieachse des G-Orbits ausgerichtet ist, sodass die Orbit-Gram-Struktur max(d12+d23−d13) = 1 erzwingt (k=4: 3-fach-Achse → Tetraeder-Gram {1,−1/3}; k=2: 4-fach-Achse → {0,±1}). Endlichkeit genügt NICHT (k=8 endlich, aber K3max = 3/√5 > 1); Clifford-Eigenschaft genügt NICHT (Oktaeder-Gruppe mit gedrehtem Q: K3max = 3/2, Anwendung iv). Der Pauli-aligned-Clifford-Fall (Theorem A) ist der **Spezialfall**, in dem die Ausrichtung automatisch das ±Pauli-Raster trifft.

**(c) Verbindung zum Magic-Befund (korrigiert).** Die Referenz misst bei k=4 Magic M₂ = 0.5585 (Lauf r1, deponiert via `p5a_r1_engine.py`) bei gleichzeitig K3 = 1.000 (inert). Das Theorem erklärt die Dissoziation strukturell: der Zeuge sieht ausschließlich die **Orbit-Gram-Struktur** der Heisenberg-Richtungen — bei k=4 den Tetraeder-Gram mit Dots {1, −1/3} (NICHT ein „diskretes ±Pauli-Raster {0,±1}"; das gilt nur bei k=2), was K3 auf 1 deckelt — nicht die Magic-Ressource. Blindheit trotz M₂ = 0.5585 ist also zustandsunabhängige Alignment-Struktur — **dieser Zeuge ist blind für die k=4-Magic**; über andere Zeugen oder andere Q-Achsen sagt das nichts.

## Abgrenzung (ehrlich): Gottesman–Knill ist NICHT das Argument

Gottesman–Knill liefert klassische Simulierbarkeit von Clifford-Schaltkreisen auf Stabilisator-Zuständen — eine Aussage über Rechenkomplexität, keine über temporale Korrelationsschranken. Unser Argument benutzt an keiner Stelle Stabilisator-Simulierbarkeit; es benutzt die **Orbit-Gram-Geometrie** der Heisenberg-Richtungen (Theorem B) bzw., im Spezialfall, die Pauli-Orbit-Diskretheit (Theorem A). Eine Brücke „simulierbar ⇒ LGI-inert" müsste separat gebaut werden und wird hier weder gebraucht noch behauptet.

---

## Ehrlichkeits-Block (Historie der Reparatur, transparent)

- **Erstfassung (2026-07-07):** behauptete das Theorem allein über „G ⊂ Clifford" und nannte als Beispiel „das endliche k=4-Braid-Bild 2T ⊂ Clifford im Z-Rahmen". Korollar formulierte Blindheit über ein „diskretes ±Pauli-Raster ⇒ Korrelatoren {0,±1}".
- **Zwei unabhängige Refuter-Pässe (2026-07-08)** fingen zwei Fehler (Fehlerklasse **R6.2** — selbstbewusst-plausibler, aber faktisch falscher Zwischenschritt):
  1. Das **reale** k=4-Fusions-Z-Bild ist KEINE Clifford-Untergruppe im Z-Rahmen: σ₁ wirkt projektiv als 3-fach-Drehung um z, Orbit(Z) = 4 Tetraeder-Richtungen mit Dots −1/3, σ₂†Zσ₂ kein ±Pauli. Korrelatoren sind {1, −1/3}, nicht {0, ±1}. (Refuter B1, `find_counterexample.py`/`ergebnis_b1.json`.)
  2. **Misalignment:** dieselbe abstrakte 2T/2O-Gruppe mit gedrehtem Q = (X+Y)/√2 liefert K3max = 3/2 — „endliche Gruppe ⇒ blind" bzw. „Treiber = Clifford" ist falsch nuanciert; Treiber ist das **Alignment**. (Refuter A-p07-Reimpl, `reimpl_ergebnis.json['gram_checks']`.)
- **Reparatur = dieses Update (2026-07-08):** Theorem A mit expliziter Pauli-aligned-Voraussetzung + Warnhinweis; neues Theorem B (Orbit-Gram) deckt das reale k=4 korrekt; Korollare (b)/(c) auf Alignment umgestellt; Spot-Checks `spot_checks_gram.py` (G1–G5), reproduziert die unabhängige A-p07-Reimpl (tet=1, octa=1, ico=3/√5, octa_misaligned=3/2) bis auf Maschinengenauigkeit. `spot_checks.py` (Theorem A) blieb unverändert und grün.

---

## Verifikations-Status

| Baustein | Status |
|---|---|
| Theorem A · Lemma A1 (24/24 Elemente 2T abstrakt; Gruppenabschluss) | verifiziert (CHECK 1) |
| Theorem A · Korollar A-a (48/48 Elemente 2O abstrakt) | verifiziert (CHECK 1b) |
| Theorem A · Lemma A2 (36/36 Pauli-Paare; ρ-Spannweite 4.44e-16) | verifiziert (CHECK 2) |
| Theorem A · Lemma A3 (36/36 Fälle, Formel == Brute-Force, 56 ρ) | verifiziert (CHECK 3) |
| Theorem A end-zu-end (576 Paare × 14 ρ, max K3 = 1.0 exakt) | verifiziert (CHECK 4) |
| Theorem A · Kontrast Nicht-Clifford (K3 = 3/2) | verifiziert (CHECK 5) |
| **Theorem B · C_ij = n_i·n_j ρ-unabhängig (Eigenwert-Spread ≤ 6.7e-16)** | **verifiziert (CHECK G1/G5)** |
| **Theorem B (i) k=4 Tetraeder-Orbit, Dots −1/3, K3max = 1 exakt** | **verifiziert (CHECK G1/G2), reales Braid-Bild** |
| **Theorem B (ii) k=2 Oktaeder aligned, Dots {0,±1}, K3max = 1** | **verifiziert (CHECK G2)** |
| **Theorem B (iii) k=8 Ikosaeder, Dots {±1,±1/√5}, K3max = 3/√5** | **verifiziert (CHECK G3), reales Braid-Bild** |
| **Theorem B (iv) Alignment: Oktaeder-Gruppe, gedrehtes Q → K3max = 3/2** | **verifiziert (CHECK G4 Gram + G5 operational)** |
| **Referenz A-p07-Reimpl (tet=1, octa=1, ico=3/√5, octa_mis=3/2) reproduziert** | **verifiziert (|diff| = 0 bzw. ≤ 1e-15)** |
| Endlichkeit 2T/2O/2I als k=4/2/8-Braid-Bild | Literatur (Jones 1986; Tuba–Wenzl math/9912013), zugeschrieben |
| Lüders-Schranke 3/2 | Literatur (Budroni–Emary 1309.3678), zugeschrieben |
| K3opt(k=8, max-Q) = 1.427 · M₂(k=4) = 0.5585 | übernommen aus der Referenz (r1/r2), hier nicht neu gerechnet |

**Scope-Grenzen (ehrlich):** Aussage gilt für den dichotomen Lüders-LGI-Zeugen K3 mit **festem** Q und Propagatoren aus einer Gruppe G ⊂ SU(2). Blindheit ist an das **Alignment** von Q mit einer Orbit-Symmetrieachse gebunden (Theorem B) — nicht an Endlichkeit oder Clifford per se; bei gedrehtem Q kann derselbe endliche/Clifford-G K3 = 3/2 erreichen. Keine Aussage über andere Messschemata (schwache/POVM-Auslese), andere Zeugen oder Multi-Qubit-Erweiterungen.

---

## Rohausgaben Theorem A (spot_checks.py, seed=42, PYTHONIOENCODING=utf-8)

```
== CHECK 1 (Lemma 1): 2T konjugiert Z auf +-Pauli ==
2T: 24 Elemente, alle unitaer mit det=1; abgeschlossen unter Produkt (24x24): True
U^dag Z U in {+-X,+-Y,+-Z} fuer 24/24 Elemente; Verteilung: [('+X', 4), ('+Y', 4), ('+Z', 4), ('-X', 4), ('-Y', 4), ('-Z', 4)]

== CHECK 1b (Korollar a): 2O (k=2-Bild) konjugiert Z auf +-Pauli ==
2O: 48 Elemente, abgeschlossen (48x48): True; U^dag Z U +-Pauli fuer 48/48 Elemente

== CHECK 2 (Lemma 2): alle 36 signierten Pauli-Paare ==
36 Paare: 12x (1/2){A,B} = +-I (gleiche Achse), 24x = 0 (verschiedene Achse)
C = (1/2)tr(rho{A,B}) ueber 206 Zufalls-rho: max. Spannweite = 4.44e-16 (rho-unabhaengig); alle Werte in {-1,0,+1}: True

== CHECK 3 (Lemma 3): Fall-Liste A-E vs. Brute-Force ueber Zufalls-rho ==
  Fall A: 4 Kombinationen, K3-Werte [-3, 1]
  Fall B: 8 Kombinationen, K3-Werte [-1, 1]
  Fall C: 8 Kombinationen, K3-Werte [-1, 1]
  Fall D: 8 Kombinationen, K3-Werte [-1, 1]
  Fall E: 8 Kombinationen, K3-Werte [0]
36/36 Kombinationen: Fall-Formel == Brute-Force (56 rho, Toleranz 1e-10); max K3 = 1

== CHECK 4 (Theorem end-zu-end): alle 576 Propagator-Paare (V1,V2) aus 2T ==
576 Paare, je 14 rho (Ginibre+rein+gemischt): Z2,Z3 stets +-Pauli: True; max K3 = 1.000000000000 <= 1 (erreicht, z.B. V1=V2=I)

== CHECK 5 (Korollar b, Kontrast): Nicht-Clifford-Praezession theta=pi/3 ==
U = exp(-i pi/6 X): U^dag Z U ist +-Pauli? None  -> Lemma 1 verletzt
K3 = 1.500000000000 fuer alle 26 rho (Spannweite 8.88e-16) = 3/2 = Lueders-Schranke (Budroni-Emary, arXiv:1309.3678) > 1

ALLE CHECKS BESTANDEN (seed=42, deterministisch)
```

## Rohausgaben Theorem B (spot_checks_gram.py, deterministisch, PYTHONIOENCODING=utf-8, Reparatur 2026-07-08)

```
== CHECK G1 (THEOREM B): Tetraeder-Orbit -- Dots -1/3, K3max = 1 (Blindheit) ==
  4 Einheitsvektoren (1,1,1)-Typ/sqrt3; alle 6 Paar-Dots = -1/3 (max|Abw.| = 1.1e-16); K3max ueber 64 Tripel = 1.000000000000000 = 1 EXAKT (Ref reimpl tet_aligned=1.0000000000000002, |diff|=0.0e+00)
  reales k=4-Bild: proj. Ordnung 12 (T); sigma1 = proj. 3-fach-Drehung um z (Ordnung 3, cos(theta)=-0.500000=-1/2); Orbit(Z) = 4 Tetraeder-Richtungen, Paar-Dots = -1/3 (max|Abw.|=2.2e-16)
  sigma2^dag Z sigma2 -> Bloch(-0.471405,-0.816497,-0.333333) = (-sqrt2/3,-sqrt6/3,-1/3), KEIN +-Pauli; K3max (Q1=Z fest) = 1.000000000000000 = 1

== CHECK G2 (THEOREM B, k=2): Oktaeder aligned z -- Dots {0,+-1}, K3max = 1 ==
  6 Achsen +-{x,y,z}; Paar-Dots (offdiag) = [-1.0, 0.0] in {0,+-1}; K3max (frei) = 1.000000000000000, (Q1=z fest) = 1.000000000000000 = 1 (Ref reimpl octa_aligned=1.0, |diff|=0.0e+00)

== CHECK G3 (THEOREM B, k=8): Ikosaeder -- Dots {+-1,+-1/sqrt5}, K3max = 3/sqrt5 ==
  12 Ecken (0,+-1,+-phi)&zykl.; Paar-Dots (offdiag) = [-1.0, -0.447214, 0.447214] in {+-1,+-1/sqrt5}; K3max = 1.341640786499874 = 3/sqrt5 = 1.341640786499874 > 1 (Ref reimpl ico=1.341640786499874, |diff|=0.0e+00)
  reales k=8-Bild: proj. Ordnung 60 (I); Orbit(Z) = 12 Ikosaeder-Richtungen; K3max (Q1=Z fest) = 1.341640786499873 = 3/sqrt5 -> Zeuge SIEHT k=8

== CHECK G4 (Alignment-Korollar): Oktaeder-Rotationsgruppe, Q=(1,1,0)/sqrt2 -> K3max = 3/2 ==
  24 Vorzeichen-Permutationsmatrizen (det=1) = chirale Oktaedergruppe O; Orbit(Q) = 12 Kantenmitten-Richtungen; Dots (offdiag) = [-1.0, -0.5, -0.0, 0.5]
  K3max (Q1=Q=n0 fest) = 1.500000000000000 = 3/2 EXAKT (Ref reimpl octa_misaligned_110=1.4999999999999996, |diff|=0.0e+00) -> ENDLICHE Gruppe, aber NICHT blind

== CHECK G5 (operational): 2 echte Qubit-Unitaries realisieren das G4-Argmax, K3 = Gram ==
  V1,V2 in SU(2) (unitaer, konstruiert ueber n2,n3); Bloch(Q2)=n2, Bloch(Q3)=n3 verifiziert
  operationaler Lueders-K3 = 1.500000000000000 fuer alle 42 rho (Spannweite 6.7e-16, rho-unabhaengig) == Gram-Formel 1.500000000000000 == 3/2

ALLE GRAM-CHECKS BESTANDEN (deterministisch; Geometrie exakt, Referenz A-p07 reproduziert)
```
