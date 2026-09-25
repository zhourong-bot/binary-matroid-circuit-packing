#!/usr/bin/env python3
"""Exact small-case checks for binary-matroid circuit packing on binary matroids.

A binary matroid on n labelled elements is stored by a row basis of its cycle
space C <= GF(2)^n.  Isomorphism classes are generated inductively from binary
linear codes, with exact coloured-graph isomorphism used for deduplication.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import time
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx
from networkx.algorithms import isomorphism


KNOWN_BINARY_MATROIDS = [1, 2, 4, 8, 16, 32, 68, 148, 342, 848, 2297]


def popcount(x: int) -> int:
    return bin(x).count("1")


def rref(rows: list[int] | tuple[int, ...], n: int) -> tuple[int, ...]:
    a = sorted(set(x for x in rows if x))
    out: list[int] = []
    for col in range(n):
        try:
            k = next(i for i, x in enumerate(a) if (x >> col) & 1)
        except StopIteration:
            continue
        pivot = a.pop(k)
        a = [x ^ pivot if (x >> col) & 1 else x for x in a]
        out = [x ^ pivot if (x >> col) & 1 else x for x in out]
        out.append(pivot)
    return tuple(sorted(out, key=lambda x: (x & -x).bit_length()))


def span_words(basis: tuple[int, ...]) -> list[int]:
    words = [0]
    for row in basis:
        words += [x ^ row for x in words]
    return words


def dual_basis(basis: tuple[int, ...], n: int) -> tuple[int, ...]:
    b = rref(basis, n)
    pivots = [(row & -row).bit_length() - 1 for row in b]
    free = [j for j in range(n) if j not in set(pivots)]
    ans = []
    for f in free:
        x = 1 << f
        for row, p in zip(b, pivots):
            if (row >> f) & 1:
                x |= 1 << p
        ans.append(x)
    return rref(ans, n)


def code_graph(basis: tuple[int, ...], n: int) -> nx.Graph:
    low = basis if len(basis) <= n - len(basis) else dual_basis(basis, n)
    g = nx.Graph()
    for j in range(n):
        g.add_node(("c", j), kind="coordinate")
    for idx, word in enumerate(sorted(span_words(low))[1:]):
        w = ("w", idx)
        g.add_node(w, kind="word")
        for j in range(n):
            if (word >> j) & 1:
                g.add_edge(w, ("c", j))
    return g


def graph_bucket(basis: tuple[int, ...], n: int, g: nx.Graph) -> tuple:
    low = basis if len(basis) <= n - len(basis) else dual_basis(basis, n)
    weights = tuple(sorted(popcount(x) for x in span_words(low)))
    h = nx.weisfeiler_lehman_graph_hash(g, node_attr="kind", iterations=5)
    return len(basis), weights, h


def extend_catalog(reps: list[tuple[int, ...]], old_n: int,
                   existing: list[tuple[int, ...]] | None = None) -> list[tuple[int, ...]]:
    n = old_n + 1
    buckets: dict[tuple, list[tuple[tuple[int, ...], nx.Graph]]] = defaultdict(list)
    node_match = isomorphism.categorical_node_match("kind", "")
    for cand in existing or []:
        g = code_graph(cand, n)
        buckets[graph_bucket(cand, n, g)].append((cand, g))
    for basis in reps:
        d = len(basis)
        candidates = []
        for functional in range(1 << d):
            rows = [row | (((functional >> i) & 1) << old_n) for i, row in enumerate(basis)]
            candidates.append(rref(rows, n))
        candidates.append(rref(list(basis) + [1 << old_n], n))
        for cand in candidates:
            g = code_graph(cand, n)
            key = graph_bucket(cand, n, g)
            if any(nx.is_isomorphic(g, old_g, node_match=node_match) for _, old_g in buckets[key]):
                continue
            buckets[key].append((cand, g))
    return [basis for group in buckets.values() for basis, _ in group]


def load_catalog(path: Path) -> tuple[int, list[tuple[int, ...]], list[int]]:
    if not path.exists():
        return 0, [tuple()], [1]
    obj = json.loads(path.read_text())
    return obj["n"], [tuple(x) for x in obj["reps"]], obj["counts"]


def save_catalog(path: Path, n: int, reps: list[tuple[int, ...]], counts: list[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    obj = {"n": n, "counts": counts, "reps": [list(x) for x in reps]}
    path.write_text(json.dumps(obj, separators=(",", ":")) + "\n")


def circuits_of(basis: tuple[int, ...]) -> list[int]:
    nonzero = sorted(set(span_words(basis)) - {0}, key=lambda x: (popcount(x), x))
    circuits = []
    for word in nonzero:
        if not any((c & word) == c for c in circuits):
            circuits.append(word)
    return circuits


def packing_number(circuits: list[int], ground: int) -> int:
    by_e = defaultdict(list)
    for c in circuits:
        for e in bits(c):
            by_e[e].append(c)
    memo = {}

    def go(avail: int) -> int:
        if avail in memo:
            return memo[avail]
        usable_e = next((e for e in bits(avail) if any((c & avail) == c for c in by_e[e])), None)
        if usable_e is None:
            return 0
        best = go(avail & ~(1 << usable_e))
        for c in by_e[usable_e]:
            if (c & avail) == c:
                best = max(best, 1 + go(avail & ~c))
        memo[avail] = best
        return best

    return go(ground)


def maximum_packings(circuits: list[int], ground: int):
    by_e = defaultdict(list)
    for c in circuits:
        for e in bits(c):
            by_e[e].append(c)
    memo = {}

    def val(avail: int) -> int:
        if avail in memo:
            return memo[avail]
        usable_e = next((e for e in bits(avail) if any((c & avail) == c for c in by_e[e])), None)
        if usable_e is None:
            memo[avail] = 0
            return 0
        best = val(avail & ~(1 << usable_e))
        for c in by_e[usable_e]:
            if (c & avail) == c:
                best = max(best, 1 + val(avail & ~c))
        memo[avail] = best
        return best

    def gen(avail: int):
        target = val(avail)
        usable_e = next((e for e in bits(avail) if any((c & avail) == c for c in by_e[e])), None)
        if usable_e is None:
            yield ()
            return
        skipped = avail & ~(1 << usable_e)
        if val(skipped) == target:
            yield from gen(skipped)
        for c in by_e[usable_e]:
            if (c & avail) == c and 1 + val(avail & ~c) == target:
                for rest in gen(avail & ~c):
                    yield (c,) + rest

    yield from gen(ground)


def bits(mask: int):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low
