#!/usr/bin/env python3
"""Bounded graph cover review: n<=6 atlas, one optimal packing/all spanning trees.

Read-only import of the existing exact engine. Outputs are written to the requested file.
"""
import hashlib
import argparse
import importlib.util
import itertools
import json
from pathlib import Path
import signal
import sys
import time

sys.dont_write_bytecode = True
import networkx as nx

HERE = Path(__file__).resolve().parent
ENGINE = HERE.parent / "round4" / "r31_bond_rules.py"
spec = importlib.util.spec_from_file_location("r31_readonly", ENGINE)
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)


def vertex_cover_number(graph):
    nodes = list(graph)
    for size in range(len(nodes) + 1):
        for chosen in itertools.combinations(nodes, size):
            cover = set(chosen)
            if all(u in cover or v in cover for u, v in graph.edges()):
                return size


def check_graph(graph):
    bonds = engine.enumerate_bonds(graph)
    nu, packing, method = engine.max_bond_packing(graph, bonds)
    delta = len(graph) - 1 - nu
    tau = vertex_cover_number(graph)
    assert tau <= 4 * delta
    records = engine.edge_records(graph)
    trees = 0
    for selected in itertools.combinations(range(len(records)), len(graph) - 1):
        tree = nx.Graph()
        tree.add_nodes_from(graph)
        tree.add_edges_from((records[i][0], records[i][1]) for i in selected)
        if not nx.is_tree(tree):
            continue
        tree_mask = sum(1 << i for i in selected)
        good = 0
        q = 0
        for bond in packing:
            intersection = bond & tree_mask
            assert intersection
            if engine.popcount(intersection) == 1:
                good |= intersection
                i = intersection.bit_length() - 1
                a, b, _ = records[i]
                cut_tree = tree.copy()
                cut_tree.remove_edge(a, b)
                shore = nx.node_connected_component(cut_tree, a)
                fundamental = sum(1 << j for j, (u, v, _) in enumerate(records)
                                  if (u in shore) != (v in shore))
                assert fundamental == bond
            else:
                q += 1
        fmask = tree_mask & ~good
        endpoints = {v for i in selected if (fmask >> i) & 1
                     for v in records[i][:2]}
        assert q <= delta and engine.popcount(fmask) == delta + q
        assert len(endpoints) <= 4 * delta
        assert all(u in endpoints or v in endpoints for u, v in graph.edges())
        trees += 1
    classes = {}
    for v in graph:
        classes.setdefault(frozenset(graph[v]), []).append(v)
    fixed = all(len(vertices) <= len(neighbors)
                for neighbors, vertices in classes.items())
    if fixed:
        assert len(graph) <= tau * 2 ** (tau - 1)
    return trees, fixed, method


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError()))
    signal.alarm(240)
    started = time.monotonic()
    counts, total_trees, fixed_points, methods = {}, 0, 0, set()
    for graph in nx.graph_atlas_g():
        if not 3 <= len(graph) <= 6 or not nx.is_biconnected(graph):
            continue
        trees, fixed, method = check_graph(graph)
        counts[len(graph)] = counts.get(len(graph), 0) + 1
        total_trees += trees
        fixed_points += fixed
        methods.add(method)
    # Local block twins are NOT necessarily twins in the glued input graph.
    main_block = nx.complete_bipartite_graph(2, 3)
    original = main_block.copy()
    original.add_edges_from([(2, 5), (5, 6), (6, 2)])
    unsafe = original.copy()
    unsafe.remove_node(2)
    reduced_block = main_block.copy()
    reduced_block.remove_node(2)
    triangle = original.subgraph([2, 5, 6]).copy()
    invariants = {name: engine.invariant(g) for name, g in [
        ("original", original), ("unsafe_global_deletion", unsafe),
        ("original_main_block", main_block),
        ("safe_reduced_main_block", reduced_block), ("retained_triangle", triangle)]}
    assert invariants["original"]["delta"] == 2
    assert invariants["unsafe_global_deletion"]["delta"] == 1
    assert invariants["safe_reduced_main_block"]["delta"] + invariants["retained_triangle"]["delta"] == 2
    result = {
        "status": "PASS", "atlas_counts_by_n": counts,
        "graphs_checked": sum(counts.values()), "spanning_trees_checked": total_trees,
        "packing_scope": "one exact optimal packing per graph; all spanning trees",
        "fixed_points_checked": fixed_points, "packing_methods": sorted(methods),
        "unsafe_interpretation_example": {"edges": sorted(map(list, original.edges())),
            "deleted_vertex": 2, "invariants": invariants},
        "engine_sha256": hashlib.sha256(ENGINE.read_bytes()).hexdigest(),
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }
    output = args.output
    with output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k != "unsafe_interpretation_example"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
