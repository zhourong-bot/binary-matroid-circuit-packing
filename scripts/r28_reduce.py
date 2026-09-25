#!/usr/bin/env python3
"""Exact gap-preserving reductions for binary matroids.

The rows are a basis of the cycle space over GF(2).  Hence primal deletion is
shortening, primal contraction is puncturing, and equal nonzero columns are a
coparallel pair in the primal matroid.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import random
import sys
import time
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import check_proposition_p as core  # noqa: E402


SEED = 2801980


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def compress_word(word: int, n: int, removed: set[int]) -> int:
    out = 0
    j = 0
    for i in range(n):
        if i not in removed:
            out |= ((word >> i) & 1) << j
            j += 1
    return out


def puncture(basis: tuple[int, ...], n: int, removed: set[int]) -> tuple[int, ...]:
    """Cycle space of primal contraction: project away removed coordinates."""
    return core.rref([compress_word(x, n, removed) for x in basis], n - len(removed))


def shorten(basis: tuple[int, ...], n: int, removed: set[int]) -> tuple[int, ...]:
    """Cycle space of primal deletion: retain words zero on removed coordinates."""
    mask = sum(1 << e for e in removed)
    words = [compress_word(x, n, removed) for x in core.span_words(basis) if not (x & mask)]
    return core.rref(words, n - len(removed))


def column_patterns(basis: tuple[int, ...], n: int) -> list[int]:
    return [sum(((row >> e) & 1) << i for i, row in enumerate(basis)) for e in range(n)]


@lru_cache(maxsize=None)
def invariant(basis: tuple[int, ...], n: int) -> tuple[int, int, int, tuple[int, ...]]:
    circuits = tuple(core.circuits_of(basis))
    nu = core.packing_number(list(circuits), (1 << n) - 1)
    return len(basis), nu, len(basis) - nu, circuits


def rule_targets(basis: tuple[int, ...], n: int):
    pats = column_patterns(basis, n)
    coloops = [i for i, p in enumerate(pats) if p == 0]
    groups = defaultdict(list)
    for i, p in enumerate(pats):
        if p:
            groups[p].append(i)
    copairs = [(a, b) for g in groups.values() for a, b in itertools.combinations(g, 2)]
    circuits = list(invariant(basis, n)[3])
    circuit_set = set(circuits)
    circuit_components = []
    for comp in core.components(circuits, (1 << n) - 1):
        mask = sum(1 << e for e in comp)
        if mask in circuit_set:
            circuit_components.append(tuple(comp))
    return coloops, copairs, circuit_components


def reduce_matroid(basis: tuple[int, ...], n: int) -> tuple[tuple[int, ...], int, list[dict]]:
    trace = []
    while n:
        coloops, copairs, circuit_components = rule_targets(basis, n)
        if coloops:
            removed = set(coloops)
            trace.append({"rule": "R1", "count": len(removed)})
            basis, n = shorten(basis, n, removed), n - len(removed)
        elif circuit_components:
            removed = set(circuit_components[0])
            trace.append({"rule": "R3", "count": len(removed)})
            basis, n = shorten(basis, n, removed), n - len(removed)
        elif copairs:
            removed = {copairs[0][0]}
            trace.append({"rule": "R2", "count": 1})
            basis, n = puncture(basis, n, removed), n - 1
        else:
            break
    return basis, n, trace


def exact_n_reps(catalog: Path, n: int) -> list[tuple[int, ...]]:
    nmax, reps, _ = core.load_catalog(catalog)
    if n > nmax:
        raise ValueError("requested n exceeds complete catalogue")
    need = nmax - n
    out = []
    for basis in reps:
        zeros = [i for i, p in enumerate(column_patterns(basis, nmax)) if p == 0]
        if len(zeros) >= need:
            out.append(puncture(basis, nmax, set(zeros[:need])))
    expected = core.KNOWN_BINARY_MATROIDS[n]
    if len(out) != expected:
        raise AssertionError(f"n={n}: recovered {len(out)} classes, expected {expected}")
    return out
