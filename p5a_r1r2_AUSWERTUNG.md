# Auswertung — B_magic r2 (LGI + 2. Maß + Dissoziations-Tabelle), 2026-06-23

**Reimplementation r1.** Engine: `magic_r2.py` (LGI exakt nach Paper 4, RoM, Vier-Wege-Tabelle).
Aufbauend auf r1 (`engine.py` + `AUSWERTUNG.md`: M₂-Tabelle + McKay verifiziert). Seed 2026, deterministisch.

## Target 1+2 (r1+r2, bestätigt): Magic-Tabelle + McKay — DREI Maße unabhängig einig
| Maß | k=2 | k=4 | k=8 | dicht |
|---|---|---|---|---|
| **M₂** (Stabilizer-2-Rényi) | 0 | 0.5585 | 0.5530 | 0.585 |
| **RoM** (Robustness of Magic) | 1.000 | 1.715 | 1.714 | √3=1.732 |
| **geom** (1−max-Stab-Fidelity) | 0.000 | 0.167 | 0.156 | 0.207 |
**Alle drei unabhängig in der CC-Engine reimplementiert** (kein Browser-Wert übernommen) → drei Stufen
(Clifford 0 / endlich-non-Clifford / dicht-Maximum), dreifach bestätigt. McKay (proj. Ord. 24/12/60, σ₁-Ord 4/3/5) → 2O/2T/2I.

## Target 3 (r2): LGI exakt nach Paper 4 — gelöst
**K₃ = 2C(B) − C(B²)**, C(U,n)=nᵀR(U)n, **Maximum über die Braid-Gruppe** (axis-/state-optimiert; K₃(ρ)=Tr[ρM]).
Reproduziert die **Paper-4-Anker exakt:** dicht→3/2, **k=8 K₃opt=1.427 (max-Q, ikos. 72°), K₃(Q=ẑ)=3/√5=1.3416**, k=2 & k=4 = 1.000.
**Zustandsabhängigkeit:** K₃(ρ)=Tr[ρM] ist zustandsabhängig; das k=4-Maximum über **alle** Zustände + Achsen + Braid-Elemente = **1.000** → genuin **blind**, strukturell inert, **KEINE Tautologie.**

## ★ Vier-Wege-Dissoziations-Tabelle (alles unabhängig reproduziert)
| k | FLW-Theorie-univ. | 3-strand image | LGI K₃opt (feuert?) | K₃(Q=ẑ) | Magic M₂ | RoM |
|---|---|---|---|---|---|---|
| 2 | nein | endlich (2O) | 1.000 (nein) | 1.000 | **0.000** | 1.000 |
| 3 | ja | dicht | 1.500 (JA) | 1.498 | 0.585 | 1.732 |
| **4** | **NEIN** (FLW-Ausnahme) | endlich (2T) | **1.000 (NEIN)** | 1.000 | **0.558** | 1.715 |
| 5 | ja | dicht | 1.500 (JA) | 1.479 | 0.585 | 1.732 |
| 6 | ja | dicht | 1.500 (JA) | 1.495 | 0.584 | 1.732 |
| 7 | ja | dicht | 1.498 (JA) | 1.494 | 0.585 | 1.732 |
| 8 | ja | endlich (2I) | **1.427 (JA)** | **1.342=3/√5** | 0.553 | 1.714 |
| 9 | ja | dicht | 1.500 (JA) | 1.497 | 0.585 | 1.732 |
| 10 | ja | dicht | 1.500 (JA) | 1.491 | 0.585 | 1.732 |

**Der Showcase k=4:** LGI feuert **NICHT** (K₃=1, jede Achse, jeder Zustand), aber Magic ist **fast maximal** (M₂=0.558, RoM=1.715). **Der Zeuge ist blind für die Ressource.** Vier dissoziable Achsen: LGI=3-strand density · Magic=Non-Cliffordness · FLW-Universalität · endlich/dicht.

## Gate-Status (magic.md)
✅ Anker reproduziert (unabhängige Engine) · ✅ 3 Maße einig (M₂/RoM/geom.) · ✅ Prior-Art-Pass grün (datiert) · ✅ Sanity (YBE≤1e-15, k=2→0, Ord. 24/12/60, K₃(k=8)=3/√5) · ✅ LGI-Zustandsabh. geklärt. → **Magic-Material steht (Sektions-Qualität).** Sammeln, nicht vorab publizieren.

## Dateien
`engine.py` (r1: Konventionen, M₂, McKay) · `magic_r2.py` (LGI, RoM, Tabelle) · `AUSWERTUNG.md` (r1) · diese (r2).
