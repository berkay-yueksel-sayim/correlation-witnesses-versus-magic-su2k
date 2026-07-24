"""
sim04c TEIL A — Hex-Zyklus-Generator (isoliertes Pre-Requisite).

Geometrische Konstruktion einer Hexagon-Plakette auf einem N×M Honeycomb-Torus.
NICHT eine ungeordnete Kantenliste — sondern ein geordneter Zyklus
A → B → A → B → A → B → (zurück zu A), 6 distinkte Vertices.

Konvention:
    Vertex-IDs:
        A(i,j) für i ∈ [0,N), j ∈ [0,M):  ID = i*M + j         (range: 0..N*M-1)
        B(i,j) für i ∈ [0,N), j ∈ [0,M):  ID = N*M + i*M + j   (range: N*M..2*N*M-1)

    Kanten-IDs (3 pro A-Site):
        x-Kante  (A(i,j) — B(i,j)):                ID = 3*(i*M + j) + 0
        y-Kante  (A(i,j) — B((i-1)%N, j)):         ID = 3*(i*M + j) + 1
        z-Kante  (A(i,j) — B(i, (j-1)%M)):         ID = 3*(i*M + j) + 2

    Damit: jede Kante ist eindeutig durch (Typ, A-Site-Index) bestimmt.
    B-Sites haben jeweils:
        x-Kante: zu A(i,j)
        y-Kante: zu A((i+1)%N, j)
        z-Kante: zu A(i, (j+1)%M)

Hexagon-Plakette (i,j) — geometrische Definition:
    Sechs-Schritt-Zyklus startend bei A(i,j), abwechselnd A→B (x/y/z) und B→A
    (entgegengesetzte Richtung), so dass nach 6 Schritten wieder bei A(i,j):

        Schritt 1: A(i,j) -- x-Kante --> B(i,j)
        Schritt 2: B(i,j) -- y-rev   --> A((i+1)%N, j)
        Schritt 3: A((i+1)%N, j) -- z-Kante --> B((i+1)%N, (j-1)%M)
        Schritt 4: B((i+1)%N, (j-1)%M) -- x-rev --> A((i+1)%N, (j-1)%M)
        Schritt 5: A((i+1)%N, (j-1)%M) -- y-Kante --> B(i, (j-1)%M)
        Schritt 6: B(i, (j-1)%M) -- z-rev --> A(i, j)   [zurück]

Diese Konstruktion garantiert:
    * 6 Schritte
    * Vertices alternieren A-B-A-B-A-B
    * Auf nicht-zu-kleinen Tori (N≥2, M≥2): 6 distinkte Vertices
    * Geschlossener Ring
"""
import json
from pathlib import Path
from collections import namedtuple


CheckResult = namedtuple("CheckResult",
                         ["passed", "n_clean", "n_total", "details"])


# === Vertex und Edge ID Funktionen ===

def A_id(i, j, N, M):
    return (i % N) * M + (j % M)

def B_id(i, j, N, M):
    return N * M + (i % N) * M + (j % M)

def edge_x(i, j, N, M):
    """A(i,j) -- B(i,j), x-Richtung"""
    return 3 * ((i % N) * M + (j % M)) + 0

def edge_y(i, j, N, M):
    """A(i,j) -- B((i-1)%N, j), y-Richtung"""
    return 3 * ((i % N) * M + (j % M)) + 1

def edge_z(i, j, N, M):
    """A(i,j) -- B(i, (j-1)%M), z-Richtung"""
    return 3 * ((i % N) * M + (j % M)) + 2


# === Hex-Zyklus-Generator (geometrisch) ===

def get_plaquette_cycle(i, j, N, M):
    """
    Geometrische Konstruktion der Plakette (i,j) als geordneter A-B-A-B-A-B Zyklus.

    Returns:
        edges: list of 6 edge-IDs in cyclic order
        vertices: list of 6 vertex-IDs in cyclic order (alternating A-B-A-B-A-B)
    """
    ip = (i + 1) % N
    jm = (j - 1) % M

    # 6 Vertices in cyclic order: A B A B A B
    A0 = A_id(i,  j,  N, M)
    B0 = B_id(i,  j,  N, M)
    A1 = A_id(ip, j,  N, M)
    B1 = B_id(ip, jm, N, M)
    A2 = A_id(ip, jm, N, M)
    B2 = B_id(i,  jm, N, M)

    vertices = [A0, B0, A1, B1, A2, B2]

    # 6 Edges in cyclic order:
    # A0 - B0:  x-Kante von A(i,j)
    # B0 - A1:  y-Kante von A(ip, j)        [B(i,j) hat y-Kante zu A(ip, j)]
    # A1 - B1:  z-Kante von A(ip, j)        [A(ip, j) z-Kante zu B(ip, (j-1)%M) = B(ip, jm) = B1]
    # B1 - A2:  x-Kante von A(ip, jm)       [B(ip, jm) hat x-Kante zu A(ip, jm)]
    # A2 - B2:  y-Kante von A(ip, jm)       [A(ip, jm) y-Kante zu B((ip-1)%N, jm) = B(i, jm) = B2]
    # B2 - A0:  z-Kante von A(i, j)         [B(i, jm) hat z-Kante zu A(i, (jm+1)%M) = A(i, j) = A0]
    edges = [
        edge_x(i,  j,  N, M),
        edge_y(ip, j,  N, M),
        edge_z(ip, j,  N, M),
        edge_x(ip, jm, N, M),
        edge_y(ip, jm, N, M),
        edge_z(i,  j,  N, M),
    ]

    return edges, vertices


# === Sauberkeits-Check ===

def get_edge_endpoints(edge_id, N, M):
    """Gibt (vertex_a, vertex_b) für eine Kante zurück."""
    e_kind = edge_id % 3
    a_idx = edge_id // 3
    i = a_idx // M
    j = a_idx % M

    A = A_id(i, j, N, M)
    if e_kind == 0:    # x: A(i,j) — B(i,j)
        B = B_id(i, j, N, M)
    elif e_kind == 1:  # y: A(i,j) — B((i-1)%N, j)
        B = B_id((i - 1) % N, j, N, M)
    else:              # z: A(i,j) — B(i, (j-1)%M)
        B = B_id(i, (j - 1) % M, N, M)
    return A, B


def check_plaquette_cleanliness(i, j, N, M):
    """
    Verifiziert die fünf Sauberkeits-Bedingungen für Plakette (i,j).
    Returns: (passed: bool, reason: str, info: dict)
    """
    edges, vertices = get_plaquette_cycle(i, j, N, M)

    # 1. Genau 6 Kanten
    if len(edges) != 6:
        return False, f"Anzahl Kanten = {len(edges)}, erwartet 6", {}

    # 2. Genau 6 Kanten distinkt
    if len(set(edges)) != 6:
        return False, f"Kanten nicht distinkt: {edges}", {}

    # 3. Genau 6 Vertices und alle distinkt
    if len(vertices) != 6:
        return False, f"Anzahl Vertices = {len(vertices)}", {}
    if len(set(vertices)) != 6:
        return False, f"Vertices nicht distinkt: {vertices}", {}

    # 4. Jede aufeinanderfolgende Kante teilt genau einen Vertex mit der vorigen
    #    UND die Endpunkte der Kanten matchen die Vertex-Sequenz
    for k in range(6):
        v_curr = vertices[k]
        v_next = vertices[(k + 1) % 6]
        e = edges[k]
        endpoints = set(get_edge_endpoints(e, N, M))
        if endpoints != {v_curr, v_next}:
            return False, (
                f"Kante {e} hat Endpunkte {endpoints}, "
                f"erwartet {{{v_curr}, {v_next}}} an Position {k}"
            ), {}

    # 5. Vertices alternieren A-B-A-B-A-B
    NM = N * M
    is_A = [v < NM for v in vertices]
    if is_A != [True, False, True, False, True, False]:
        return False, f"Vertices alternieren nicht A-B-A-B-A-B: A={is_A}", {}

    # 6. Plus: kein Vertex berührt mehr als 2 Plaketten-Kanten
    incidence = {v: 0 for v in vertices}
    for e in edges:
        for v in get_edge_endpoints(e, N, M):
            if v in incidence:
                incidence[v] += 1
    max_inc = max(incidence.values())
    if max_inc != 2:
        return False, (
            f"Vertex-Inzidenz nicht uniform 2: max={max_inc}, "
            f"verteilung={dict(incidence)}"
        ), {}

    return True, "OK", {"edges": edges, "vertices": vertices,
                        "incidence": incidence}


def check_all_plaquettes(N, M):
    """
    Prüft alle N*M Plaketten auf Sauberkeit.
    Returns: CheckResult
    """
    n_total = N * M
    n_clean = 0
    details = []
    for i in range(N):
        for j in range(M):
            ok, reason, info = check_plaquette_cleanliness(i, j, N, M)
            details.append({
                "plaquette": [i, j],
                "passed": ok,
                "reason": reason,
            })
            if ok:
                n_clean += 1
    all_clean = (n_clean == n_total)
    return CheckResult(passed=all_clean, n_clean=n_clean,
                       n_total=n_total, details=details)


# === Pflicht-Test-Matrix ===

def main():
    print("=" * 70)
    print("sim04c TEIL A — Hex-Zyklus-Generator")
    print("=" * 70)
    print()

    test_grids = [(2, 2), (2, 3), (3, 2), (3, 3), (4, 4), (5, 5)]
    summary = {}
    smallest_clean = None

    for N, M in test_grids:
        result = check_all_plaquettes(N, M)
        status = "PASS" if result.passed else "FAIL"
        key = f"{N}x{M}"
        summary[key] = {
            "all_clean": result.passed,
            "n_clean": result.n_clean,
            "n_total": result.n_total,
        }
        if result.passed and smallest_clean is None:
            smallest_clean = key
        print(f"NxM={N}x{M}: {status}  ({result.n_clean}/{result.n_total} Plaketten sauber)")
        if not result.passed:
            for d in result.details:
                if not d["passed"]:
                    print(f"    p{d['plaquette']}: {d['reason']}")
                    break  # nur eine Beispiel-Diagnose

    summary["smallest_clean_NxM"] = smallest_clean

    print()
    if smallest_clean is None:
        print("STOP-KRITERIUM: Kein NxM ≤ 5x5 sauber.")
        print("Diagnose erforderlich, kein weiterer Code.")
    else:
        print(f"Kleinstes sauberes Gitter: {smallest_clean}")
        # Visualisiere eine Plakette des kleinsten sauberen Gitters
        N_s, M_s = map(int, smallest_clean.split("x"))
        print(f"\nBeispiel-Plakette (0,0) auf {smallest_clean}:")
        edges, vertices = get_plaquette_cycle(0, 0, N_s, M_s)
        NM = N_s * M_s
        for k in range(6):
            v = vertices[k]
            kind = "A" if v < NM else "B"
            idx = v if v < NM else v - NM
            i, j = idx // M_s, idx % M_s
            e = edges[k]
            print(f"  Step {k+1}: {kind}({i},{j}) [vertex {v}] -- edge {e}")
        v0 = vertices[0]
        kind0 = "A" if v0 < NM else "B"
        idx0 = v0 if v0 < NM else v0 - NM
        i0, j0 = idx0 // M_s, idx0 % M_s
        print(f"  -- back to {kind0}({i0},{j0}) [vertex {v0}]")

    # JSON output
    out_path = Path(__file__).parent / "hex_cycle_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"\nWritten: {out_path.name}")
    return summary


if __name__ == "__main__":
    main()
