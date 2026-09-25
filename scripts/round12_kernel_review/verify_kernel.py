#!/usr/bin/env python3
"""Kernel pipeline verification. New outputs only; run with python3 -B.

Exact invariants are verification oracles, never reduction-selection inputs.
Usage: python3 -B verify_kernel.py --output NEW_DIRECTORY --seconds 1800
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import random
import signal
import sys
import time

SCRIPTS = Path(__file__).resolve().parents[1]
PROJECT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))
import check_proposition_p as core
import r28_reduce as red

def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / relative)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

old = module('backsolve', 'round9_matroid_backsolve/run.py')
t25 = module('exchange', 'round11_h1_proof/verify_exchange.py')
pc = core.popcount
bits = lambda x: tuple(core.bits(x))

@lru_cache(maxsize=30000)
def exact(rows, n):
    cs = tuple(core.circuits_of(rows))
    nu = core.packing_number(list(cs), (1 << n) - 1)
    return (len(rows), nu, len(rows) - nu)

def fundamental(rows):
    return {(r & -r).bit_length() - 1: r for r in rows}

def conflicts(fund):
    return tuple((a, b) for a, b in itertools.combinations(sorted(fund), 2)
                 if fund[a] & fund[b])

def cover_ok(fset, edges):
    return all(a in fset or b in fset for a, b in edges)

def minimum_cover(fund, k):
    if k < 0:
        return None
    es = conflicts(fund)
    for size in range(min(2*k, len(fund)) + 1):
        for v in itertools.combinations(sorted(fund), size):
            if cover_ok(set(v), es):
                return tuple(v)
    return None

def valid_state(rows, n, fund, fset):
    star = sum(1 << e for e in fund)
    assert len(fund) == len(rows)
    assert core.rref(tuple(fund.values()), n) == rows
    assert all(c & star == 1 << e for e, c in fund.items())
    good = {e: c for e, c in fund.items() if e not in fset}
    assert all(not (a & b) for a, b in itertools.combinations(good.values(), 2))
    return star, good

def h1_targets(fund, fset):
    star = sum(1 << e for e in fund)
    good = {e: c for e, c in fund.items() if e not in fset}
    return list(old.candidates((), 0, star, fund, good, tuple(sorted(fset))))

def contract(rows, n, fund, fset, removed):
    keep = [e for e in range(n) if e not in removed]
    ren = {e: i for i, e in enumerate(keep)}
    rows2 = red.puncture(rows, n, removed)
    fund2 = {ren[e]: red.compress_word(c, n, removed)
             for e, c in fund.items() if e not in removed}
    assert not set(fset) & removed
    fset2 = tuple(ren[e] for e in fset)
    valid_state(rows2, len(keep), fund2, fset2)
    return rows2, len(keep), fund2, fset2, keep

def bound_ok(n, k):
    # g(k) >= 2^(2k). Avoid constructing double-exponential-sized integers.
    if k < 0:
        return False
    if 2*k >= max(1, n.bit_length()):
        return True
    a = 1 << (2*k)
    if a >= max(1, n.bit_length()) and k:
        return True
    return n <= 2*k + a + 2*k*(1 << a)*a

def structural_bound(rows, n, fund, fset):
    star, good = valid_state(rows, n, fund, fset)
    pats = red.column_patterns(tuple(fund.values()), n)
    assert all(pats) and len(set(pats)) == n
    groups = Counter(old.type_of(c, fset, fund, star) for c in good.values())
    h = len(fset)
    assert all(count <= h for count in groups.values())
    p0 = ((1 << n)-1) & ~star
    for c in good.values():
        typ = old.type_of(c, fset, fund, star)
        assert len(set(typ)) == len(typ) and 0 not in typ
        assert pc(c) <= 1 << h
        p0 &= ~c
    assert pc(p0) <= (1 << h)-1
    # Stronger instance-sensitive bound; no gigantic exponent required.
    assert n == h + pc(p0) + sum(pc(c) for c in good.values())
    return {'F': h, 'P0': pc(p0), 'blocks': len(good),
            'types': len(groups), 'largest_class': max(groups.values(), default=0)}

def run_kernel(rows, n, k, label, trace, totals):
    initial = exact(rows, n)
    fund = fundamental(rows)
    labels = list(range(n))
    steps = 0
    rounds = 0
    first_F = None
    while True:
        rounds += 1
        fset = minimum_cover(fund, k)
        if fset is None:
            assert initial[2] > k
            return {'k': k, 'answer': False, 'rejected': True,
                    'reject_round': rounds, 'n': n, 'steps': steps}
        if first_F is None:
            first_F = len(fset)
        star, good = valid_state(rows, n, fund, fset)
        candidates = h1_targets(fund, fset)
        witness = None
        if candidates:
            typ, e, c, t = candidates[0]
            removed = set(bits(c))
            rule = 'H1'
            witness = {'F': [labels[f] for f in fset], 't': t,
                       'type': list(typ), 'e': labels[e]}
        else:
            pats = red.column_patterns(tuple(fund.values()), n)
            zeros = [e for e, p in enumerate(pats) if p == 0]
            if zeros:
                removed = {zeros[0]}
                rule = 'R1'
                assert not (star & (1 << zeros[0]))
            else:
                groups = defaultdict(list)
                for e, p in enumerate(pats):
                    groups[p].append(e)
                pair = next((v[:2] for v in groups.values() if len(v)>1), None)
                if pair:
                    e = next(x for x in pair if not (star & (1 << x)))
                    removed = {e}
                    rule = 'R2'
                    witness = {'pair': [labels[x] for x in pair]}
                else:
                    structure = structural_bound(rows, n, fund, fset)
                    assert bound_ok(n, k)
                    assert exact(rows, n)[2] == initial[2]
                    return {'k': k, 'answer': initial[2] <= k, 'rejected': False,
                            'n': n, 'steps': steps, 'rows': rows,
                            'labels': labels, 'initial_F': first_F,
                            'structure': structure}
        before = exact(rows, n)
        nextrows, nextn, nextfund, nextF, keep = contract(rows, n, fund, fset, removed)
        after = exact(nextrows, nextn)
        event = {'case': label, 'k': k, 'step': steps+1, 'rule': rule,
                 'removed': [labels[x] for x in sorted(removed)],
                 'n_before': n, 'n_after': nextn, 'before': before,
                 'after': after, 'witness': witness}
        trace.write(json.dumps(event)+'\n')
        totals['reductions'] += 1
        totals['rule_'+rule] += 1
        if before[2] != after[2]:
            totals['delta_failures'] += 1
            raise AssertionError(event)
        shift = 1 if rule == 'H1' else 0
        assert after[0] == before[0]-shift and after[1] == before[1]-shift
        assert nextn < n and len(nextF) == len(fset)
        labels = [labels[x] for x in keep]
        rows, n, fund = nextrows, nextn, nextfund
        steps += 1
        # Fresh minimum cover on the maintained basis at each round.

def other_covers(fund, k, mincover, rng, limit=4):
    if mincover is None or len(mincover) >= min(2*k, len(fund)):
        return []
    es = conflicts(fund)
    # Reservoir sampling from ALL eligible nonminimum covers, not only
    # supersets of the particular minimum cover chosen by the pipeline.
    chosen = []
    seen = 0
    for size in range(len(mincover)+1, min(2*k, len(fund))+1):
        for cand in itertools.combinations(sorted(fund), size):
            if cover_ok(set(cand), es):
                seen += 1
                if len(chosen) < limit:
                    chosen.append(cand)
                else:
                    j = rng.randrange(seen)
                    if j < limit:
                        chosen[j] = cand
    return chosen

def audit_nonminimum(rows, n, delta, label, rng, log, totals):
    fund = fundamental(rows)
    before = exact(rows, n)
    result = Counter()
    for k in (delta, delta-1):
        mincover = minimum_cover(fund, k)
        covers = other_covers(fund, k, mincover, rng)
        for fset in covers:
            result['covers'] += 1
            valid_state(rows, n, fund, fset)
            nonoptimal = len(fund)-len(fset) < before[1]
            result['nonoptimal_families'] += nonoptimal
            targets = h1_targets(fund, fset)
            result['triggering_covers'] += bool(targets)
            for typ, e, c, t in targets:
                nextrows, nextn, _, _, _ = contract(rows, n, fund, fset, set(bits(c)))
                after = exact(nextrows, nextn)
                event = {'case': label, 'k': k, 'F': fset,
                         'minimum_cover_size': len(mincover),
                         'family_size': len(fund)-len(fset), 'nu': before[1],
                         'nonoptimal': nonoptimal, 'e': e, 't': t,
                         'type': typ, 'before': before, 'after': after}
                log.write(json.dumps(event)+'\n')
                result['h1_tests'] += 1
                result['nonoptimal_h1_tests'] += nonoptimal
                totals['reductions'] += 1
                totals['nonminimum_H1'] += 1
                if before[2] != after[2]:
                    totals['delta_failures'] += 1
                    raise AssertionError(event)
                assert after[:2] == (before[0]-1, before[1]-1)
    totals.update({'nonminimum_'+key: val for key, val in result.items()})
    return dict(result)

def cases(rng):
    catalog = PROJECT/'data/binary_codes_catalog.json'
    for n in range(11):
        for i, rows in enumerate(red.exact_n_reps(catalog, n)):
            yield 'catalog:%d:%d'%(n,i), n, rows
    for n in range(11,17):
        for i in range(32):
            d = i % (n+1)
            # A random systematic generator, then random coordinate labels.
            rows = [((1 << j) | (rng.getrandbits(n-d) << d)) for j in range(d)]
            perm = list(range(n)); rng.shuffle(perm)
            rows = core.rref([sum(1 << perm[e] for e in bits(r)) for r in rows],n)
            yield 'random:%d:%d'%(n,i), n, rows
    # Exchange-step generators with duplicate/empty signatures and residual base elements.
    j=0
    for h in range(4):
        for t in range(h+1,h+4):
            for p in range(4):
                for residual in (0,2):
                    c = t25.make_case(rng,h,t,p,residual,1,'structured:%d'%j)
                    if c and 11 <= c['n'] <=16:
                        yield c['name'], c['n'], core.rref(tuple(c['fund'].values()),c['n'])
                        j += 1
    # Nonoptimal independent families with many identical loop blocks.
    for n in range(11,17):
        for d in (1,2,3):
            rows=[]
            for i in range(d):
                rows.extend([3 << (3*i),5 << (3*i)])
            rows.extend(1 << j for j in range(3*d,n))
            yield 'loops_plus_parallel:%d:%d'%(n,d), n, core.rref(rows,n)
    for a in (3,4):
        for b in range(1,9):
            rows, edges = old.kab(a,b)
            yield 'cographic:K%d,%d'%(a,b),len(edges),rows

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--seconds',type=int,default=1800)
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    def timeout(signum, frame):
        raise TimeoutError('Foreground computation wall-clock limit')
    signal.signal(signal.SIGALRM,timeout)
    signal.alarm(args.seconds)
    totals=Counter(instances=0,reductions=0,delta_failures=0)
    maxima={str(d): {'n':0,'case':None} for d in (1,2,3)}
    counts=Counter()
    status='PASS'
    error=None
    current=None
    seed=26092512
    rng=random.Random(seed)
    source_rng=random.Random(seed+1)
    sources=[PROJECT/'data/binary_codes_catalog.json',
             SCRIPTS/'check_proposition_p.py',SCRIPTS/'r28_reduce.py',
             SCRIPTS/'round9_matroid_backsolve/run.py',
             SCRIPTS/'round11_h1_proof/verify_exchange.py',Path(__file__)]
    hashes={p.relative_to(PROJECT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    try:
        with (args.output/'instances.jsonl').open('x') as cases_out, \
             (args.output/'reductions.jsonl').open('x') as trace, \
             (args.output/'nonminimum.jsonl').open('x') as other:
            for label,n,rows in cases(source_rng):
                current=label
                rows=core.rref(rows,n)
                before=exact(rows,n)
                delta=before[2]
                totals['instances']+=1
                counts[label.split(':')[0]]+=1
                if label.startswith('catalog:'):
                    totals['catalog_n_'+str(n)]+=1
                if n<=6 or label.startswith('random:') and int(label.split(':')[-1])<3:
                    cs=core.circuits_of(rows)
                    assert len(t25.packing(cs,(1<<n)-1))==before[1]
                    totals['independent_packing_checks']+=1
                outcomes=[]
                for k in (delta,delta-1):
                    res=run_kernel(rows,n,k,label,trace,totals)
                    assert res['answer']==(delta<=k)
                    if k==delta:
                        assert not res['rejected']
                        if str(delta) in maxima and res['n']>maxima[str(delta)]['n']:
                            maxima[str(delta)]={'n':res['n'],'case':label,'rows':res['rows']}
                    elif res['rejected']:
                        totals['NO_rejected']+=1
                        totals['NO_rejected_first']+=res['reject_round']==1
                    else:
                        totals['NO_retained']+=1
                    totals['kernel_runs']+=1
                    outcomes.append(res)
                nonminimum=audit_nonminimum(rows,n,delta,label,rng,other,totals)
                cases_out.write(json.dumps({'case':label,'n':n,'rows':rows,
                    'invariants':before,'runs':outcomes,'nonminimum':nonminimum})+'\n')
                cases_out.flush()
                if totals['instances']%250==0 or label.startswith('cographic:'):
                    print(json.dumps({'progress':totals['instances'],'case':label,
                        'seconds':round(time.monotonic()-started,2)}),flush=True)
    except Exception as ex:
        status='FAIL' if not isinstance(ex,TimeoutError) else 'TIMEOUT'
        error={'type':type(ex).__name__,'message':str(ex),'case':current}
    finally:
        signal.alarm(0)
    summary={'status':status,'seed':seed,'seconds':time.monotonic()-started,
             'totals':dict(totals),'counts':dict(counts),'max_E_by_delta':maxima,
             'error':error,'source_sha256':hashes,'exact_cache':str(exact.cache_info())}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,ensure_ascii=False),flush=True)
    return 0 if status=='PASS' else 1

if __name__=='__main__':
    sys.exit(main())
