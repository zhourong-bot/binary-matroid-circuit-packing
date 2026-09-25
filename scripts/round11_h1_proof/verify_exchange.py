#!/usr/bin/env python3.11
"""Finite audit of the exchange step; stdlib only, no network/subprocesses.

Requires Python >= 3.10; this host uses /opt/homebrew/bin/python3.11.
Usage: /opt/homebrew/bin/python3.11 -B verify_exchange.py --output NEW_FILE.json
Keeps original element labels after contraction. Does not import old scripts.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import itertools
import json
from pathlib import Path
import random
import time


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def xor_all(values):
    result = 0
    for value in values:
        result ^= value
    return result


def words(rows):
    span = {0}
    for row in rows:
        span |= {word ^ row for word in span}
    return span


def circuits(span):
    result = []
    for word in sorted(span - {0}, key=lambda x: (x.bit_count(), x)):
        if not any(word & c == c for c in result):
            result.append(word)
    return tuple(result)


def packing(circs, ground):
    by_element = {e: tuple(c for c in circs if c & (1 << e)) for e in bits(ground)}

    @lru_cache(None)
    def solve(avail):
        if not avail:
            return ()
        low = avail & -avail
        e = low.bit_length() - 1
        best = solve(avail ^ low)
        for c in by_element[e]:
            if c & avail == c:
                candidate = (c,) + solve(avail ^ c)
                if len(candidate) > len(best):
                    best = candidate
        return best

    return solve(ground)


def all_packings(circs, start=0, occupied=0, current=()):
    yield current
    for i in range(start, len(circs)):
        c = circs[i]
        if not c & occupied:
            yield from all_packings(circs, i + 1, occupied | c, current + (c,))


def rank(vectors):
    pivots = {}
    for x in vectors:
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                break
    return len(pivots)


def primal_columns(n, star, fund):
    """A=[I|fundamental-incidences], independent rank-oracle cross-check."""
    basis = tuple(bits(((1 << n) - 1) ^ star))
    columns = [0] * n
    for j, e in enumerate(basis):
        columns[e] = 1 << j
    for e in bits(star):
        columns[e] = sum(1 << j for j, b in enumerate(basis) if fund[e] & (1 << b))
    return columns


def rank_circuits(columns, contracted=0):
    available = tuple(e for e in range(len(columns)) if not contracted & (1 << e))
    fixed = [columns[e] for e in bits(contracted)]
    base = rank(fixed)
    result = []
    for size in range(1, len(available) + 1):
        for subset in itertools.combinations(available, size):
            mask = sum(1 << e for e in subset)
            if any(mask & c == c for c in result):
                continue
            if rank(fixed + [columns[e] for e in subset]) - base < size:
                result.append(mask)
    return set(result)


def audit_case(case, totals, exhaustive=False, independent=False):
    n, fund, good, fset, selected = (case[x] for x in ('n', 'fund', 'good', 'F', 'selected'))
    ground = (1 << n) - 1
    star = sum(1 << e for e in fund)
    fmask = sum(1 << f for f in fset)
    block_union = sum(good[e] for e in selected)
    all_good_union = sum(good.values())
    assert all(fund[e] & star == 1 << e for e in fund)
    assert all(good[e] == fund[e] for e in good)
    assert sum(c.bit_count() for c in good.values()) == all_good_union.bit_count()
    assert set(fund) == set(good) | set(fset)
    assert not set(good) & set(fset)
    k, t = len(fset), len(selected)
    assert t > k
    sigs = []
    for e in selected:
        sigs.append(sorted(sum(1 << j for j, f in enumerate(fset) if fund[f] & (1 << x))
                           for x in bits(good[e] & ~star)))
    assert all(s == sigs[0] for s in sigs)
    span = words(fund.values())
    assert len(span) == 1 << len(fund)
    original = circuits(span)
    assert all(c in original for c in good.values())
    original_set = set(original)
    original_pack = packing(original, ground)
    neutral_original = {c for c in original if not c & block_union}
    neutral_number = len(packing(tuple(neutral_original), ground & ~block_union))
    assert len(original_pack) == t + neutral_number
    totals['cases'] += 1
    totals['max_n'] = max(totals['max_n'], n)
    totals['cases_with_nonzero_R'] += int(any(fund[f] & ~(star | all_good_union) for f in fset))
    totals['cases_with_other_good_blocks'] += int(len(good) > t)
    totals['cases_with_repeated_signatures'] += int(len(sigs[0]) != len(set(sigs[0])))
    totals['cases_with_empty_signature'] += int(0 in sigs[0])
    totals['cases_k_zero'] += int(k == 0)
    totals['original_circuits'] += len(original)
    if independent:
        columns = primal_columns(n, star, fund)
        assert rank_circuits(columns) == original_set
        totals['independent_rank_originals'] += 1

    def local(z, e):
        return xor_all(fund[f] for f in fset if z & (1 << f)) & good[e]

    for removed_e in selected:
        removed = good[removed_e]
        remaining = ground & ~removed
        survivor_keys = tuple(e for e in selected if e != removed_e)
        stars = {good[e] for e in survivor_keys}
        surviving_union = block_union & ~removed
        projected_span = {word & ~removed for word in span}
        minor_circuits = circuits(projected_span)
        minor_set = set(minor_circuits)
        assert len(projected_span) * 2 == len(span)
        assert stars <= minor_set
        neutral_minor = set()
        for c in minor_circuits:
            z = c & fmask
            a = local(z, removed_e)
            if c in stars:
                assert not z
                totals['star_circuits'] += 1
            elif a:
                assert z
                assert all(c & good[e] for e in survivor_keys)
                assert all((c & good[e]) in (local(z, e), good[e] ^ local(z, e))
                           for e in survivor_keys)
                totals['crossing_circuits'] += 1
            else:
                assert not c & surviving_union
                assert c in original_set
                neutral_minor.add(c)
                totals['neutral_circuits'] += 1
        assert neutral_minor == neutral_original
        witness = packing(minor_circuits, remaining)
        assert len(witness) == len(original_pack) - 1
        totals['contractions'] += 1
        if independent:
            assert rank_circuits(columns, removed) == minor_set
            totals['independent_rank_contractions'] += 1
        tested = all_packings(minor_circuits) if exhaustive else [witness]
        for pack in tested:
            q = [c for c in pack if local(c & fmask, removed_e)]
            m = [c for c in pack if c in stars]
            keep = tuple(c for c in pack if c not in q and c not in m)
            assert len(q) <= k
            assert not (q and m)
            assert len(q) + len(m) <= t - 1
            lifted = keep + tuple(good[e] for e in selected)
            assert len(lifted) >= len(pack) + 1
            assert all(c in original_set for c in lifted)
            assert sum(c.bit_count() for c in lifted) == xor_all(lifted).bit_count()
            totals['exchanged_packings'] += 1
            if q and len(q) == k == t - 1:
                totals['tight_boundary_exchanges'] += 1
    return {'name': case['name'], 'n': n, 'k': k, 't': t,
            'nu': len(original_pack), 'neutral_nu': neutral_number,
            'all_minor_packings': exhaustive, 'rank_oracle': independent}


def make_case(rng, k, t, p, residual, extras, name):
    """Systematic fundamental basis with independently permuted signature copies."""
    good, parts, cursor = {}, {}, 0
    selected = []
    for i in range(t + extras):
        length = p if i < t else rng.randint(0, 2)
        e = cursor
        cursor += length + 1
        good[e] = ((1 << (length + 1)) - 1) << e
        parts[e] = list(range(e + 1, cursor))
        if i < t:
            selected.append(e)
    p0 = list(range(cursor, cursor + residual))
    cursor += residual
    fset = list(range(cursor, cursor + k))
    cursor += k
    if cursor > 16:
        return None
    fund = dict(good)
    signatures = [rng.randrange(1 << k) for _ in range(p)]
    copies = {}
    for e in selected:
        perm = list(signatures)
        rng.shuffle(perm)
        copies[e] = perm
    for j, f in enumerate(fset):
        c = 1 << f
        for e in selected:
            for x, sig in zip(parts[e], copies[e]):
                if sig & (1 << j):
                    c |= 1 << x
        for e in good:
            if e not in selected:
                for x in parts[e]:
                    if rng.getrandbits(1):
                        c |= 1 << x
        for x in p0:
            if rng.getrandbits(1):
                c |= 1 << x
        fund[f] = c
    return {'name': name, 'n': cursor, 'fund': fund, 'good': good,
            'F': fset, 'selected': selected}


def equality_counterexample(path):
    example = json.loads(path.read_text())['examples']['t2']
    n, rows, star, removed = (example[x] for x in ('n', 'rows', 'Bstar', 'D'))
    span = words(rows)
    fund = {e: next(w for w in span if w & star == 1 << e) for e in bits(star)}
    original = circuits(span)
    minor = circuits({w & ~removed for w in span})
    old_pack = packing(original, (1 << n) - 1)
    new_pack = packing(minor, ((1 << n) - 1) & ~removed)
    assert len(old_pack) == len(new_pack) == 2
    assert len(span) == 16 and len({w & ~removed for w in span}) == 8
    columns = primal_columns(n, star, fund)
    assert rank_circuits(columns) == set(original)
    assert rank_circuits(columns, removed) == set(minor)
    return {'n': n, 'fundamental_rows': fund, 'Bstar': star, 'F': example['F'],
            'original_good_pack': example['pack'], 'D': removed,
            'minor_packing': new_pack, 'before': example['before'], 'after': example['after']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must be a new file')
    start = time.monotonic()
    rng = random.Random(250925)
    totals = Counter()
    cases = []
    # A tight q=k=t-1 example. Third block reverses signature order; R is nonempty.
    cases.append({'name': 'tight_permuted_nonzero_R', 'n': 13,
                  'fund': {1: 70, 3: 152, 8: 1792,
                           0: 21 | (1 << 10) | (1 << 11),
                           5: 224 | (1 << 9) | (1 << 12)},
                  'good': {1: 70, 3: 152, 8: 1792}, 'F': [0, 5], 'selected': [1, 3, 8]})
    for k, p, residual, extras in itertools.product(range(5), range(4), range(3), range(2)):
        t = k + 1
        c = make_case(rng, k, t, p, residual, extras,
                      f'structured_k{k}_p{p}_r{residual}_extra{extras}')
        if c is not None:
            cases.append(c)
    results = []
    for i, case in enumerate(cases):
        if time.monotonic() - start > 180:
            raise RuntimeError('front-end time bound reached; no partial PASS written')
        results.append(audit_case(case, totals, exhaustive=case['n'] <= 10 or i == 0,
                                  independent=i == 0 or (case['n'] <= 8 and i % 11 == 0)))
    root = Path(__file__).resolve().parents[2]
    boundary = equality_counterexample(root / 'data/round10/threshold_minima.json')
    assert totals['tight_boundary_exchanges'] > 0
    output = {'status': 'PASS', 'seed': 250925,
              'scope': 'Proof-step audit of constructed fundamental representations, not a matroid census; no optimality assumption on input good family.',
              'totals': dict(totals), 'cases': results,
              'equal_threshold_counterexample': boundary,
              'elapsed_seconds': round(time.monotonic() - start, 3)}
    with args.output.open('x') as stream:
        json.dump(output, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    print(json.dumps({key: output[key] for key in ('status', 'totals', 'elapsed_seconds')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
