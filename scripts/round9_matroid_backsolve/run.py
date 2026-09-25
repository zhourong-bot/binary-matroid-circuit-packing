#!/usr/bin/env python3
"""Exact small binary-matroid backsolve and candidate reductions.

Cycle-space rows use the existing binary representation convention.  No stochastic sampling.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import check_proposition_p as core  # noqa: E402
import r28_reduce as red  # noqa: E402


def bits(x):
    return tuple(core.bits(x))


def witness(circuits, n):
    by_e = [[] for _ in range(n)]
    for c in circuits:
        for e in bits(c):
            by_e[e].append(c)

    @lru_cache(None)
    def go(avail):
        e = next((e for e in bits(avail) if any(c & avail == c for c in by_e[e])), None)
        if e is None:
            return ()
        best = go(avail & ~(1 << e))
        for c in by_e[e]:
            if c & avail == c:
                z = (c,) + go(avail & ~c)
                if len(z) > len(best):
                    best = z
        return best

    return go((1 << n) - 1)


def bases_data(rows, n, selected_star=None):
    d = len(rows)
    words = core.span_words(rows)
    stars = ([selected_star] if selected_star is not None else
             (sum(1 << e for e in subset) for subset in itertools.combinations(range(n), d)))
    for star in stars:
        fund = {}
        for w in words[1:]:
            z = w & star
            if z and core.popcount(z) == 1:
                fund[z.bit_length() - 1] = w
        if len(fund) == d:
            yield star, fund


def backsolve(rows, n, star, fund, pack):
    d = len(rows)
    nu = len(pack)
    delta = d - nu
    good = {bits(c & star)[0]: c for c in pack if core.popcount(c & star) == 1}
    fset = tuple(e for e in bits(star) if e not in good)
    q = sum(core.popcount(c & star) >= 2 for c in pack)
    assert q <= delta and len(fset) == delta + q <= 2 * delta
    assert all(fund[e] == c for e, c in good.items())
    assert all(x & y == 0 for i, x in enumerate(good.values()) for y in list(good.values())[i + 1:])
    return good, fset, q


def type_of(c, fset, fund, star):
    # Multiset of F-incidence signatures of the B-elements in one good circuit.
    # The unique B* element is excluded; F labels stay fixed.
    signatures = []
    for x in bits(c & ~star):
        signatures.append(sum((1 << j) for j, f in enumerate(fset) if fund[f] & (1 << x)))
    return tuple(sorted(signatures))


def candidates(rows, n, star, fund, good, fset):
    groups = defaultdict(list)
    for e, c in good.items():
        groups[type_of(c, fset, fund, star)].append((e, c))
    for typ, members in groups.items():
        if len(members) <= len(fset):
            continue
        for e, c in members:
            yield typ, e, c, len(members)


def graph_rows(vertices, edges, dual=False):
    rows = []
    for v in range(vertices - 1):
        rows.append(sum(1 << i for i, (u, w) in enumerate(edges) if v in (u, w)))
    cut = core.rref(rows, len(edges))
    return cut if dual else core.dual_basis(cut, len(edges))


def kab(a, b):
    edges = [(i, a + j) for i in range(a) for j in range(b)]
    return graph_rows(a + b, edges, True), edges


def graph_cases(max_vertices=8, max_edges=10):
    # geng emits one graph per isomorphism class, including disconnected graphs.
    import subprocess
    import networkx as nx
    for v in range(1, max_vertices + 1):
        proc = subprocess.Popen(['geng', '-q', str(v)], stdout=subprocess.PIPE, text=True)
        for line in proc.stdout:
            g = nx.from_graph6_bytes(line.strip().encode())
            edges = tuple(sorted((min(u, w), max(u, w)) for u, w in g.edges()))
            if len(edges) <= max_edges:
                yield v, edges, graph_rows(v, edges, False)
        if proc.wait():
            raise RuntimeError('geng failed')


def check_one(rows, n, source, all_bases=True, selected_star=None, selected_pack=None,
              all_packings=False):
    circuits = core.circuits_of(rows)
    packs = ([selected_pack] if selected_pack is not None else
             list(core.maximum_packings(circuits, (1 << n) - 1)) if all_packings else
             [witness(circuits, n)])
    delta = len(rows) - len(packs[0])
    bases = 0
    triggers = 0
    failures = []
    tested_targets = set()
    for star, fund in bases_data(rows, n, selected_star):
        bases += 1
        for pack in packs:
            good, fset, _ = backsolve(rows, n, star, fund, pack)
            for typ, e, c, count in candidates(rows, n, star, fund, good, fset):
                # Operation: contract the entire good circuit, reducing r* by 1.
                key = c
                if key in tested_targets:
                    continue
                tested_targets.add(key)
                nr = red.puncture(rows, n, set(bits(c)))
                nn = n - core.popcount(c)
                new_delta = len(nr) - core.packing_number(core.circuits_of(nr), (1 << nn) - 1)
                triggers += 1
                if new_delta != delta:
                    failures.append({'source': source, 'n': n, 'rows': rows, 'Bstar': star,
                                     'pack': pack, 'F': fset, 'type': typ, 'type_count': count,
                                     'circuit': c, 'before': delta, 'after': new_delta})
        if not all_bases:
            break
    return bases, len(packs), triggers, failures, delta


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stage', choices=['catalog', 'family', 'graphs'], required=True)
    p.add_argument('--version', type=int, choices=[1], default=1)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--max-vertices', type=int, default=8)
    p.add_argument('--max-edges', type=int, default=10)
    p.add_argument('--all-packings', action='store_true')
    args = p.parse_args()
    start = time.time()
    stats = {'stage': args.stage, 'version': args.version, 'classes': 0, 'bases': 0,
             'packing_certificates': 0,
             'triggers': 0, 'failures': 0, 'triggers_by_delta': {},
             'fixed_classes_by_delta': {}, 'max_n_fixed_by_delta': {}, 'counterexamples': [],
             'family': []}
    if args.stage == 'catalog':
        nmax, reps, counts = core.load_catalog(Path(__file__).resolve().parents[2] / 'data/binary_codes_catalog.json')
        assert nmax == 10 and sum(counts) == 3766
        cases = ((n, rows, f'catalog:{n}:{i}', None, None)
                 for n in range(0, 11)
                 for i, rows in enumerate(red.exact_n_reps(Path(__file__).resolve().parents[2] / 'data/binary_codes_catalog.json', n)))
    elif args.stage == 'family':
        def families():
            seen = set()
            for a in (3, 4):
                for b in range(3, 9):
                    aa, bb = min(a, b), max(a, b)
                    if (aa, bb) in seen:
                        continue  # K(4,3) is isomorphic to K(3,4).
                    seen.add((aa, bb))
                    rows, edges = kab(aa, bb)
                    # Spanning tree: A0--all right vertices; R0--other left vertices.
                    tree = sum(1 << i for i, (u, v) in enumerate(edges)
                               if u == 0 or v == aa)
                    pack = tuple(sum(1 << i for i, (_, v) in enumerate(edges) if v == aa + j)
                                 for j in range(bb))
                    yield a * b, rows, f'K:{a}:{b}', tree, pack
        cases = families()
    else:
        cases = ((len(edges), rows, f'graph:{v}:{edges}', None, None)
                 for v, edges, rows in graph_cases(args.max_vertices, args.max_edges))
    for n, rows, source, selected_star, selected_pack in cases:
        bc, pc, tc, failures, delta = check_one(rows, n, source,
                                            all_bases=args.stage != 'family',
                                            selected_star=selected_star, selected_pack=selected_pack,
                                            all_packings=args.all_packings)
        stats['classes'] += 1
        stats['bases'] += bc
        stats['packing_certificates'] += pc
        stats['triggers'] += tc
        if tc:
            key = str(delta)
            stats['triggers_by_delta'][key] = stats['triggers_by_delta'].get(key, 0) + tc
        stats['failures'] += len(failures)
        stats['counterexamples'].extend(failures[:max(0, 12 - len(stats['counterexamples']))])
        if tc == 0 and delta in (1, 2, 3):
            key = str(delta)
            stats['fixed_classes_by_delta'][key] = stats['fixed_classes_by_delta'].get(key, 0) + 1
            stats['max_n_fixed_by_delta'][key] = max(stats['max_n_fixed_by_delta'].get(key, 0), n)
        if args.stage == 'family':
            stats['family'].append({'source': source, 'n': n, 'delta': delta, 'bases': bc,
                                    'triggers': tc, 'failures': len(failures)})
    stats['seconds'] = round(time.time() - start, 3)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(stats, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({k: v for k, v in stats.items() if k not in ('counterexamples', 'family')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
