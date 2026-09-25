#!/usr/bin/env python3
"""Compare freshly computed results and supplied summaries with the paper."""

import argparse
import json
from pathlib import Path


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    results = args.results

    catalog = read(results / "catalog.json")
    assert (catalog["classes"], catalog["bases"], catalog["failures"]) == (3766, 134860, 0)
    family = read(results / "family.json")
    assert family["failures"] == 0 and max(item["n"] for item in family["family"]) == 32
    graphs_path = results / "matroid_graphs.json"
    if graphs_path.exists():
        graphs = read(graphs_path)
        assert (graphs["classes"], graphs["bases"], graphs["failures"]) == (8166, 2462020, 0)
    exchange = read(results / "exchange.json")
    assert exchange["status"] == "PASS" and exchange["totals"]["cases"] == 95
    kernel = read(results / "kernel/summary.json")
    totals = kernel["totals"]
    assert kernel["status"] == "PASS"
    assert (totals["instances"], totals["kernel_runs"], totals["reductions"],
            totals["delta_failures"], totals["nonminimum_H1"]) == (4021, 8042, 15716, 0, 1372)
    assert totals["reductions"] - totals["nonminimum_H1"] == 14344
    assert {int(k): item["n"] for k, item in kernel["max_E_by_delta"].items()} == {1: 9, 2: 13, 3: 20}
    cover = read(results / "cover.json")
    assert (cover["status"], cover["graphs_checked"], cover["spanning_trees_checked"]) == (
        "PASS", 70, 9588)

    census = read(root / "data/round7/biconnected_summary.json")
    assert census["graphs"] == 201727
    assert sum(census["by_n"].values()) == 201727
    assert all(2 * int(d) >= int(t) and 4 * int(d) >= int(t)
               for (t, d), count in ((key.split(","), value)
                                      for key, value in census["by_tau_delta"].items()) if count)
    baseline = read(root / "data/round12/summary.json")
    assert baseline["status"] == "PASS" and baseline["totals"]["instances"] == 4021
    print("PASS catalog=3766 classes/134860 bases; kernel=4021 inputs/8042 parameters")
    print("PASS local=15716 (14344 pipeline+1372 exchange), deficit failures=0")
    print("PASS graph census=201727, tau<=2delta and tau<=4delta in supplied summary")
    print("PASS cover certificates=70 graphs/9588 trees; terminal maxima=9,13,20")


if __name__ == "__main__":
    main()
