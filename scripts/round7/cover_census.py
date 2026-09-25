#!/usr/bin/env python3
"""Exact cover/BondGap census. Each n=9 shard is a separate foreground unit."""
import argparse
import json
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parents[2]
OUT = None
sys.path.insert(0, str(ROOT / "scripts" / "round4"))
import r31_bond_rules as engine

GENG = "geng"


def tau(g):
    n = len(g)
    adj = [sum(1 << j for j in g.neighbors(i)) for i in range(n)]
    independent = [False] * (1 << n)
    independent[0] = True
    best = 0
    for mask in range(1, 1 << n):
        bit = mask & -mask
        v = bit.bit_length() - 1
        rest = mask ^ bit
        independent[mask] = independent[rest] and not (adj[v] & rest)
        if independent[mask]:
            best = max(best, bin(mask).count("1"))
    return n - best


def process_shard(arg):
    n, residue, modulus = arg
    OUT.mkdir(parents=True, exist_ok=True)
    name = f"biconnected_n{n}_r{residue:03d}of{modulus}.json"
    path = OUT / name
    if path.exists():
        return (name, "existing")
    cmd = [GENG, "-Cq", str(n)]
    if modulus > 1:
        cmd.append(f"{residue}/{modulus}")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    count = 0
    by_tau = defaultdict(lambda: {"count": 0, "min_delta": None, "example": None})
    by_pair = defaultdict(int)
    begin = time.monotonic()
    for line in proc.stdout:
        code = line.strip()
        if not code:
            continue
        g = nx.from_graph6_bytes(code.encode("ascii"))
        t = tau(g)
        d = engine.invariant(g)["delta"]
        count += 1
        row = by_tau[t]
        row["count"] += 1
        if row["min_delta"] is None or d < row["min_delta"]:
            row["min_delta"] = d
            row["example"] = code
        by_pair[(t, d)] += 1
    proc.stdout.close()
    stderr = proc.stderr.read()
    proc.stderr.close()
    if proc.wait() != 0:
        raise RuntimeError(stderr[-1000:])
    data = {"n": n, "residue": residue, "modulus": modulus, "graphs": count,
            "elapsed_seconds": round(time.monotonic() - begin, 2),
            "by_tau": {str(k): v for k, v in sorted(by_tau.items())},
            "by_tau_delta": {f"{t},{d}": c for (t, d), c in sorted(by_pair.items())}}
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
    tmp.replace(path)
    return (name, count, data["elapsed_seconds"])


def aggregate():
    required = [OUT / f"biconnected_n{n}_r000of1.json" for n in range(3, 9)]
    required += [OUT / f"biconnected_n9_r{r:03d}of256.json" for r in range(256)]
    missing = [p.name for p in required if not p.exists()]
    if missing:
        raise RuntimeError(f"missing {len(missing)} shards: {missing[:5]}")
    pair = defaultdict(int)
    by_tau = defaultdict(lambda: {"count": 0, "min_delta": None, "example": None})
    by_n = defaultdict(int)
    for path in required:
        obj = json.loads(path.read_text())
        by_n[obj["n"]] += obj["graphs"]
        for key, val in obj["by_tau_delta"].items():
            pair[key] += val
        for key, row in obj["by_tau"].items():
            t = by_tau[int(key)]
            t["count"] += row["count"]
            if t["min_delta"] is None or row["min_delta"] < t["min_delta"]:
                t["min_delta"] = row["min_delta"]
                t["example"] = row["example"]
    total = sum(by_n.values())
    if total != 201727 or by_n[9] != 194066:
        raise AssertionError((total, by_n[9]))
    result = {"graphs": total, "by_n": dict(sorted(by_n.items())),
              "by_tau": {str(t): row for t, row in sorted(by_tau.items())},
              "by_tau_delta": dict(sorted(pair.items(), key=lambda kv: tuple(map(int, kv[0].split(','))))) }
    assert all(int(t) <= 2 * int(d) and int(t) <= 4 * int(d)
               for t, d in (key.split(",") for key, count in pair.items() if count))
    expected = json.loads((ROOT / "data/round7/biconnected_summary.json").read_text())
    assert result["graphs"] == expected["graphs"]
    assert result["by_n"] == {int(k): v for k, v in expected["by_n"].items()}
    assert result["by_tau_delta"] == expected["by_tau_delta"]
    (OUT / "biconnected_summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"graphs": total, "by_n": result["by_n"], "by_tau": result["by_tau"]}))


def verify_small():
    """Recompute the n=3..6 part of the census directly from geng."""
    expected = json.loads((ROOT / "data/round7/biconnected_summary.json").read_text())
    counts = {}
    for n in range(3, 7):
        proc = subprocess.Popen([GENG, "-Cq", str(n)], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True)
        count = 0
        for line in proc.stdout:
            if not line.strip():
                continue
            graph = nx.from_graph6_bytes(line.strip().encode("ascii"))
            deficit = engine.invariant(graph)["delta"]
            cover = tau(graph)
            assert cover <= 2 * deficit
            assert cover <= 4 * deficit
            count += 1
        proc.stdout.close()
        stderr = proc.stderr.read()
        proc.stderr.close()
        assert proc.wait() == 0, stderr
        assert count == expected["by_n"][str(n)]
        counts[n] = count
    print(json.dumps({"status": "PASS", "small_graphs": sum(counts.values()),
                      "by_n": counts, "bounds": ["tau<=2delta", "tau<=4delta"]}))

def init_worker(output_dir):
    global OUT
    OUT = output_dir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["census", "aggregate", "small"])
    ap.add_argument("--output-dir", type=Path)
    ap.add_argument("--n", type=int)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    global OUT
    OUT = args.output_dir
    if args.action == "small":
        verify_small()
    elif OUT is None:
        ap.error("--output-dir is required for census and aggregate")
    elif args.action == "aggregate":
        aggregate()
    else:
        from multiprocessing import Pool
        if args.n is None:
            raise SystemExit("--n required")
        modulus = 256 if args.n == 9 else 1
        stop = args.stop if args.stop is not None else modulus
        work = [(args.n, r, modulus) for r in range(args.start, stop)]
        if args.jobs == 1:
            for item in work:
                print(json.dumps(process_shard(item)), flush=True)
        else:
            with Pool(args.jobs, initializer=init_worker, initargs=(OUT,)) as pool:
                for row in pool.imap_unordered(process_shard, work):
                    print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
