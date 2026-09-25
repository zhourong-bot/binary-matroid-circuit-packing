#!/usr/bin/env python3
"""Exact bond packing for simple graphs.

The invariant is the cographic-matroid gap

    delta(G) = |V(G)| - c(G) - nu_bond(G),

where nu_bond is the maximum number of pairwise edge-disjoint bonds.  Bonds
are enumerated exactly as connected bipartitions of each connected component;
packing is solved by exact dynamic programming for small instances and by a
zero-gap HiGHS MILP for larger instances.


"""

from __future__ import print_function

import argparse
import copy
import itertools
import json
import math
import random
import time
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

import networkx as nx


SEED = 2801980
HERE = Path(__file__).resolve().parent


def popcount(x):
    return bin(int(x)).count("1")


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def components_count(graph):
    return nx.number_connected_components(graph) if len(graph) else 0


def simple_normalize(graph):
    out = nx.Graph()
    out.add_nodes_from(range(len(graph)))
    relabel = {v: i for i, v in enumerate(sorted(graph.nodes()))}
    out = nx.relabel_nodes(graph, relabel, copy=True)
    if isinstance(out, nx.MultiGraph):
        out = nx.Graph(out)
    out.remove_edges_from(nx.selfloop_edges(out))
    return out


def edge_records(graph):
    """Stable non-loop edge identities; parallel copies have distinct keys."""
    rows = []
    if graph.is_multigraph():
        for u, v, k in graph.edges(keys=True):
            if u != v:
                rows.append((min(u, v), max(u, v), str(k)))
    else:
        for u, v in graph.edges():
            if u != v:
                rows.append((min(u, v), max(u, v), None))
    return sorted(rows)


def enumerate_bonds(graph):
    """Return all bonds as edge bitmasks, without complementary duplication."""
    records = edge_records(graph)
    bonds = set()
    for component in nx.connected_components(graph):
        vertices = sorted(component)
        if len(vertices) < 2:
            continue
        root = vertices[0]
        rest = vertices[1:]
        for bits in range(1 << len(rest)):
            shore = {root}
            shore.update(rest[i] for i in range(len(rest)) if (bits >> i) & 1)
            if len(shore) == len(vertices):
                continue
            other = set(vertices) - shore
            if not nx.is_connected(graph.subgraph(shore)):
                continue
            if not nx.is_connected(graph.subgraph(other)):
                continue
            mask = 0
            for i, row in enumerate(records):
                u, v = row[0], row[1]
                if (u in shore) != (v in shore):
                    mask |= 1 << i
            if mask:
                bonds.add(mask)
    return tuple(sorted(bonds, key=lambda x: (popcount(x), x)))


def packing_dp(bonds):
    ordered = tuple(sorted(set(bonds), key=lambda x: (popcount(x), x)))

    @lru_cache(None)
    def solve(index, used):
        if index == len(ordered):
            return 0
        best = solve(index + 1, used)
        bond = ordered[index]
        if not (bond & used):
            best = max(best, 1 + solve(index + 1, used | bond))
        return best

    optimum = solve(0, 0)
    witness = []
    index, used = 0, 0
    while index < len(ordered):
        skip = solve(index + 1, used)
        bond = ordered[index]
        take = -1
        if not (bond & used):
            take = 1 + solve(index + 1, used | bond)
        if take >= skip and take == solve(index, used):
            witness.append(bond)
            used |= bond
        index += 1
    return optimum, witness, "exact_bitmask_dp"


def packing_milp(bonds, edge_count):
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import csc_matrix

    if not bonds:
        return 0, [], "exact_trivial"
    rows, cols = [], []
    for j, bond in enumerate(bonds):
        for e in range(edge_count):
            if (bond >> e) & 1:
                rows.append(e)
                cols.append(j)
    data = np.ones(len(rows), dtype=float)
    matrix = csc_matrix((data, (rows, cols)), shape=(edge_count, len(bonds)))
    result = milp(
        c=-np.ones(len(bonds), dtype=float),
        integrality=np.ones(len(bonds), dtype=int),
        bounds=Bounds(np.zeros(len(bonds)), np.ones(len(bonds))),
        constraints=LinearConstraint(matrix, -np.inf, np.ones(edge_count)),
        options={"mip_rel_gap": 0.0, "presolve": True},
    )
    if result.status != 0 or result.x is None:
        raise RuntimeError("HiGHS did not certify optimality: status=%r message=%r" %
                           (result.status, result.message))
    witness = [bonds[i] for i, value in enumerate(result.x) if value > 0.5]
    if any(a & b for i, a in enumerate(witness) for b in witness[i + 1:]):
        raise AssertionError("MILP witness is not edge-disjoint")
    optimum = int(round(-float(result.fun)))
    if optimum != len(witness):
        raise AssertionError("MILP objective/witness mismatch")
    return optimum, witness, "exact_highs_milp_zero_gap"


def max_bond_packing(graph, bonds=None):
    bonds = enumerate_bonds(graph) if bonds is None else tuple(bonds)
    m = len(edge_records(graph))
    if len(bonds) <= 80 and m <= 24:
        return packing_dp(bonds)
    return packing_milp(bonds, m)


def invariant(graph):
    bonds = enumerate_bonds(graph)
    nu, witness, method = max_bond_packing(graph, bonds)
    rank = len(graph) - components_count(graph)
    return {
        "vertices": len(graph),
        "edges": len(edge_records(graph)),
        "components": components_count(graph),
        "cographic_corank": rank,
        "bond_count": len(bonds),
        "nu_bond": nu,
        "delta": rank - nu,
        "packing_method": method,
        "packing_witness_masks": witness,
    }
