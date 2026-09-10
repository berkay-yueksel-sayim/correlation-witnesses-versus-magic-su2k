# Evaluation — B_magic r2 (LGI + 2nd measure + dissociation table), 2026-06-23

**Reimplementation r1.** Engine: `p5a_r2_dissociation_table.py` (LGI exactly as in Paper 4, RoM, four-way table).
Building on r1 (`p5a_r1_engine.py` + `AUSWERTUNG.md`: M₂ table + McKay verified). Seed 2026, deterministic.

## Target 1+2 (r1+r2, confirmed): magic table + McKay — THREE measures independently in agreement
| Measure | k=2 | k=4 | k=8 | dense |
|---|---|---|---|---|
| **M₂** (stabilizer 2-Rényi) | 0 | 0.5585 | 0.5530 | 0.585 |
| **RoM** (Robustness of Magic) | 1.000 | 1.715 | 1.714 | √3=1.732 |
| **geom** (1−max stab. fidelity) | 0.000 | 0.167 | 0.156 | 0.207 |
**All three independently reimplemented in the CC engine** (no browser value adopted) → three levels
(Clifford 0 / finite non-Clifford / dense maximum), confirmed threefold. McKay (proj. order 24/12/60, σ₁ order 4/3/5) → 2O/2T/2I.

## Target 3 (r2): LGI exactly as in Paper 4 — solved
**K₃ = 2C(B) − C(B²)**, C(U,n)=nᵀR(U)n, **maximum over the braid group** (axis-/state-optimized; K₃(ρ)=Tr[ρM]).
Reproduces the **Paper 4 anchors exactly:** dense→3/2, **k=8 K₃opt=1.427 (max-Q, icos. 72°), K₃(Q=ẑ)=3/√5=1.3416**, k=2 & k=4 = 1.000.
**State dependence:** K₃(ρ)=Tr[ρM] is state-dependent; the k=4 maximum over **all** states + axes + braid elements = **1.000** → genuinely **blind**, structurally inert, **NOT a tautology.**

## ★ Four-way dissociation table (everything independently reproduced)
| k | FLW theory-univ. | 3-strand image | LGI K₃opt (fires?) | K₃(Q=ẑ) | magic M₂ | RoM |
|---|---|---|---|---|---|---|
| 2 | no | finite (2O) | 1.000 (no) | 1.000 | **0.000** | 1.000 |
| 3 | yes | dense | 1.500 (YES) | 1.498 | 0.585 | 1.732 |
| **4** | **NO** (FLW exception) | finite (2T) | **1.000 (NO)** | 1.000 | **0.558** | 1.715 |
| 5 | yes | dense | 1.500 (YES) | 1.479 | 0.585 | 1.732 |
| 6 | yes | dense | 1.500 (YES) | 1.495 | 0.584 | 1.732 |
| 7 | yes | dense | 1.498 (YES) | 1.494 | 0.585 | 1.732 |
| 8 | yes | finite (2I) | **1.427 (YES)** | **1.342=3/√5** | 0.553 | 1.714 |
| 9 | yes | dense | 1.500 (YES) | 1.497 | 0.585 | 1.732 |
| 10 | yes | dense | 1.500 (YES) | 1.491 | 0.585 | 1.732 |

**The showcase k=4:** LGI does **NOT** fire (K₃=1, every axis, every state), but magic is **nearly maximal** (M₂=0.558, RoM=1.715). **The witness is blind to the resource.** Four dissociable axes: LGI=3-strand density · magic=non-Cliffordness · FLW universality · finite/dense.

## Gate status (magic.md — internal working document, not deposited)
✅ anchors reproduced (independent engine) · ✅ 3 measures in agreement (M₂/RoM/geom.) · ✅ prior-art pass green (dated) · ✅ sanity (YBE≤1e-15, k=2→0, order 24/12/60, K₃(k=8)=3/√5) · ✅ LGI state dependence clarified. → **The magic material stands (section quality).** Collect, do not pre-publish.

## Files
`p5a_r1_engine.py` (r1: conventions, M₂, McKay) · `p5a_r2_dissociation_table.py` (LGI, RoM, table) · `AUSWERTUNG.md` (r1 predecessor version, not deposited) · this one (r2).
