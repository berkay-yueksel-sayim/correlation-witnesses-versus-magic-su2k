"""
sim04c PART A — hex cycle generator (isolated pre-requisite).

Geometric construction of a hexagon plaquette on an N×M honeycomb torus.
NOT an unordered edge list — but an ordered cycle
A → B → A → B → A → B → (back to A), 6 distinct vertices.

Convention:
    Vertex IDs:
        A(i,j) for i ∈ [0,N), j ∈ [0,M):  ID = i*M + j         (range: 0..N*M-1)
        B(i,j) for i ∈ [0,N), j ∈ [0,M):  ID = N*M + i*M + j   (range: N*M..2*N*M-1)

    Edge IDs (3 per A site):
        x edge  (A(i,j) — B(i,j)):                 ID = 3*(i*M + j) + 0
        y edge  (A(i,j) — B((i-1)%N, j)):          ID = 3*(i*M + j) + 1
        z edge  (A(i,j) — B(i, (j-1)%M)):          ID = 3*(i*M + j) + 2

    Thus: every edge is uniquely determined by (type, A site index).
    B sites each have:
        x edge: to A(i,j)
        y edge: to A((i+1)%N, j)
        z edge: to A(i, (j+1)%M)

Hexagon plaquette (i,j) — geometric definition:
    Six-step cycle starting at A(i,j), alternating A→B (x/y/z) and B→A
    (opposite direction), such that after 6 steps we are back at A(i,j):

        Step 1: A(i,j) -- x edge --> B(i,j)
        Step 2: B(i,j) -- y-rev  --> A((i+1)%N, j)
        Step 3: A((i+1)%N, j) -- z edge --> B((i+1)%N, (j-1)%M)
        Step 4: B((i+1)%N, (j-1)%M) -- x-rev --> A((i+1)%N, (j-1)%M)
        Step 5: A((i+1)%N, (j-1)%M) -- y edge --> B(i, (j-1)%M)
        Step 6: B(i, (j-1)%M) -- z-rev --> A(i, j)   [back]

This construction guarantees:
    * 6 steps
    * vertices alternate A-B-A-B-A-B
    * on not-too-small tori (N≥2, M≥2): 6 distinct vertices
    * a closed ring
"""
import json
from pathlib import Path
from collections import namedtuple


CheckResult = namedtuple("CheckResult",
                         ["passed", "n_clean", "n_total", "details"])


# === vertex and edge ID functions ===

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
    Geometric construction of plaquette (i,j) as an ordered A-B-A-B-A-B cycle.

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
    # A0 - B0:  x edge of A(i,j)
    # B0 - A1:  y edge of A(ip, j)          [B(i,j) has a y edge to A(ip, j)]
    # A1 - B1:  z edge of A(ip, j)          [A(ip, j) z edge to B(ip, (j-1)%M) = B(ip, jm) = B1]
    # B1 - A2:  x edge of A(ip, jm)         [B(ip, jm) has an x edge to A(ip, jm)]
    # A2 - B2:  y edge of A(ip, jm)         [A(ip, jm) y edge to B((ip-1)%N, jm) = B(i, jm) = B2]
    # B2 - A0:  z edge of A(i, j)           [B(i, jm) has a z edge to A(i, (jm+1)%M) = A(i, j) = A0]
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
    """Returns (vertex_a, vertex_b) for an edge."""
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
    Verifies the five cleanliness conditions for plaquette (i,j).
    Returns: (passed: bool, reason: str, info: dict)
    """
    edges, vertices = get_plaquette_cycle(i, j, N, M)

    # 1. Genau 6 Kanten
    if len(edges) != 6:
        return False, f"number of edges = {len(edges)}, expected 6", {}

    # 2. Genau 6 Kanten distinkt
    if len(set(edges)) != 6:
        return False, f"edges not distinct: {edges}", {}

    # 3. exactly 6 vertices and all distinct
    if len(vertices) != 6:
        return False, f"number of vertices = {len(vertices)}", {}
    if len(set(vertices)) != 6:
        return False, f"vertices not distinct: {vertices}", {}

    # 4. each consecutive edge shares exactly one vertex with the previous one
    #    AND the edge endpoints match the vertex sequence
    for k in range(6):
        v_curr = vertices[k]
        v_next = vertices[(k + 1) % 6]
        e = edges[k]
        endpoints = set(get_edge_endpoints(e, N, M))
        if endpoints != {v_curr, v_next}:
            return False, (
                f"edge {e} has endpoints {endpoints}, "
                f"expected {{{v_curr}, {v_next}}} at position {k}"
            ), {}

    # 5. Vertices alternieren A-B-A-B-A-B
    NM = N * M
    is_A = [v < NM for v in vertices]
    if is_A != [True, False, True, False, True, False]:
        return False, f"vertices do not alternate A-B-A-B-A-B: A={is_A}", {}

    # 6. plus: no vertex touches more than 2 plaquette edges
    incidence = {v: 0 for v in vertices}
    for e in edges:
        for v in get_edge_endpoints(e, N, M):
            if v in incidence:
                incidence[v] += 1
    max_inc = max(incidence.values())
    if max_inc != 2:
        return False, (
            f"vertex incidence not uniformly 2: max={max_inc}, "
            f"distribution={dict(incidence)}"
        ), {}

    return True, "OK", {"edges": edges, "vertices": vertices,
                        "incidence": incidence}


def check_all_plaquettes(N, M):
    """
    Checks all N*M plaquettes for cleanliness.
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
                    break  # only one example diagnosis

    summary["smallest_clean_NxM"] = smallest_clean

    print()
    if smallest_clean is None:
        print("STOP CRITERION: no NxM ≤ 5x5 is clean.")
        print("Diagnosis required, no further code.")
    else:
        print(f"Kleinstes sauberes Gitter: {smallest_clean}")
        # visualize one plaquette of the smallest clean lattice
        N_s, M_s = map(int, smallest_clean.split("x"))
        print(f"\nExample plaquette (0,0) on {smallest_clean}:")
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
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
