#!/usr/bin/env python3
"""Check the constructive tau <= 4 delta certificate on geng graphs."""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "round4"))
import r31_bond_rules as engine


def check(g):
    inv = engine.invariant(g)
    records = engine.edge_records(g)
    tree = nx.minimum_spanning_tree(g)
    tree_indices = {i for i, (u, v, _) in enumerate(records) if tree.has_edge(u, v)}
    good = set()
    bad = 0
    for bond in inv["packing_witness_masks"]:
        hits = [i for i in tree_indices if bond >> i & 1]
        assert hits
        if len(hits) == 1:
            good.add(hits[0])
        else:
            bad += 1
    assert len(good) == inv["nu_bond"] - bad
    leftover = tree_indices - good
    cover = {v for i in leftover for v in records[i][:2]}
    assert len(leftover) == inv["delta"] + bad
    assert bad <= inv["delta"]
    assert len(cover) <= 4 * inv["delta"]
    assert all(u in cover or v in cover for u, v in g.edges())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--res", type=int, default=0)
    ap.add_argument("--mod", type=int, default=1)
    args = ap.parse_args()
    cmd = ["geng", "-Cq", str(args.n)]
    if args.mod > 1:
        cmd.append(f"{args.res}/{args.mod}")
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    count = 0
    for line in p.stdout:
        if line.strip():
            check(nx.from_graph6_bytes(line.strip().encode()))
            count += 1
    p.stdout.close()
    stderr = p.stderr.read()
    p.stderr.close()
    assert p.wait() == 0, stderr
    print(json.dumps({"n": args.n, "res": args.res, "mod": args.mod,
                      "certificates_checked": count}))


if __name__ == "__main__":
    main()
