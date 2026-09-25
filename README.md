# Binary matroid circuit packing: computational checks

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22952617.svg)](https://doi.org/10.5281/zenodo.22952617)


This repository accompanies the **Computational verification** section of
[`paper/main.tex`](paper/main.tex) and the corresponding
[`paper/main.pdf`](paper/main.pdf). The paper's bibliography is in
[`paper/refs.bib`](paper/refs.bib). The scripts test finite instances of the
binary matroid and graph arguments. **These implementation checks are not
proofs of the lemmas, bounds, or algorithm.** In particular, the observed
graph inequality `tau <= 2*delta` is not claimed as a theorem.

## Environment

Use Python **3.10 or newer**; the code uses `int.bit_count()`. Python 3.11 is
recommended. Install the packages in `requirements.txt` in your Python 3.11
environment. `networkx` is used throughout; `numpy` and `scipy` provide the
exact integer optimization backend for larger graph instances. The full graph
enumeration and graph-based matroid check also require `geng` from nauty on
`PATH`. Install nauty separately; this repository contains no executables.
The census aggregation uses no `labelg` or external data beyond its own
generated shards.

Run the standard suite from the repository root:

```sh
./run_all.sh --quick
```

This runs the complete 3,766-class matroid catalogue, the selected cographic
family, the exchange check, all 4,021 kernel inputs at both parameter values,
the 70-graph/9,588-tree cover check, and a fresh `geng` census for graphs with
3–6 vertices. It also verifies the supplied 201,727-graph aggregate and its
`tau <= 2*delta` and `tau <= 4*delta` counts. On the development machine,
this took **13 seconds** with Python 3.11 and the packages already installed.
`./run_all.sh` additionally checks graphic matroid examples through eight
vertices and at most 14 edges. It took **171 seconds** on the development
machine, within the 15-minute target.
Temporary logs and detailed output are removed after a successful run.

Expected final lines of the quick suite:

```text
PASS catalog=3766 classes/134860 bases; kernel=4021 inputs/8042 parameters
PASS local=15716 (14344 pipeline+1372 exchange), deficit failures=0
PASS graph census=201727, tau<=2delta and tau<=4delta in supplied summary
PASS cover certificates=70 graphs/9588 trees; terminal maxima=9,13,20
```

## Individual checks

The commands below write to a new `results/` directory, which is ignored by
Git. Run each command only once per named output file, or choose a new output
name. Reported times are from the quick suite on the development machine and
will vary by hardware.

| Check | Reproduction command from the repository root | Expected result | Approximate time |
| --- | --- | --- | --- |
| Catalogue and counting/backsolve | `python3.11 -B scripts/round9_matroid_backsolve/run.py --stage catalog --output results/catalog.json` | 3,766 classes, 134,860 bases, zero failures | 4 s |
| Selected cographic family | `python3.11 -B scripts/round9_matroid_backsolve/run.py --stage family --output results/family.json` | 11 cases, up to 32 elements, zero failures | 3 s |
| Exchange step | `python3.11 -B scripts/round11_h1_proof/verify_exchange.py --output results/exchange.json` | PASS, 95 constructed cases | <1 s |
| Kernel pipeline | `python3.11 -B scripts/round12_kernel_review/verify_kernel.py --output results/kernel --seconds 600` | PASS, 4,021 inputs, 8,042 parameter choices, 15,716 local checks, zero deficit failures; terminal maxima 9, 13, 20 | 5 s |
| Cover argument | `python3.11 -B scripts/round8_review/t21_check.py --output results/cover.json` | PASS, 70 graphs and 9,588 spanning trees | 1 s |
| Small graph census | `python3.11 -B scripts/round7/cover_census.py small` | PASS, 70 graphs, both inequalities | <1 s |
| Graphic matroid cases | `python3.11 -B scripts/round9_matroid_backsolve/run.py --stage graphs --max-vertices 8 --max-edges 14 --output results/matroid_graphs.json` | 8,166 cases, 2,462,020 bases, zero failures | 157 s |

Create `results/` before running individual commands. The kernel command
creates its own output directory and writes detailed JSON Lines there; those
large per-instance logs are not distributed. The small JSON summaries under
`data/` are provided so the paper's reported totals can be checked without
repeating the full graph enumeration. The 47 KB catalogue is provided at
`data/binary_codes_catalog.json`.

## Full graph census

`scripts/round7/cover_census.py` is the exact `geng` census used for the graph
statistics. To independently regenerate all 201,727 two-connected simple
graphs on 3–9 vertices, run the following. Each output shard contains only
counts and examples, not per-graph records. The nine-vertex census is much
heavier than the standard suite; no sub-15-minute runtime is claimed for it.

```sh
mkdir -p results/graph
for n in 3 4 5 6 7 8; do
  python3.11 -B scripts/round7/cover_census.py census --n "$n" --jobs 1 --output-dir results/graph
done
python3.11 -B scripts/round7/cover_census.py census --n 9 --jobs 8 --output-dir results/graph
python3.11 -B scripts/round7/cover_census.py aggregate --output-dir results/graph
```

Aggregation checks the total and the complete `(tau, delta)` count table
against `data/round7/biconnected_summary.json`. It fails if any graph violates
either recorded inequality. The cover-certificate check is separate:
`python3.11 -B scripts/round7/check_tree_certificate.py --n 6` checks its
constructive bound on all 56 two-connected six-vertex graphs.

The scripts represent binary matroids by a basis of the circuit space over
GF(2), use row reduction to canonicalize it, enumerate minimal nonzero words
as circuits, and compute exact disjoint-circuit packing numbers. `puncture`
implements contraction in this representation. The graph code enumerates
bonds as connected bipartitions and computes exact disjoint-bond packings.

## Citation

Rong Zhou, *Packing circuits below the corank: a kernel for binary matroids and edge-disjoint bond packing*, preprint, 2026. Zenodo, doi:[10.5281/zenodo.22952617](https://doi.org/10.5281/zenodo.22952617). ORCID [0009-0004-8825-6103](https://orcid.org/0009-0004-8825-6103).
