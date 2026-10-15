"""Command line: mutation-hunter fetch | phewas | gwas | fdr."""

import argparse
import sys
from pathlib import Path

import pandas as pd

from . import __version__, analysis, data, report, stats


def _phenotype(arg: str) -> pd.Series:
    """A condition name from the panel, or a TSV file: strain<TAB>value with a header line."""
    if Path(arg).is_file():
        y = pd.read_csv(arg, sep="\t", index_col=0).iloc[:, 0].rename(Path(arg).stem)
        if y.index.intersection(data.load_distance().index).empty:
            raise ValueError(f"{arg}: no strain ID in common with the 1011 panel (expected IDs like AAA, SACE_YBA)")
        return y
    Y = data.load_phenotypes()
    if arg not in Y:
        raise ValueError(f"unknown phenotype {arg!r}; run 'mutation-hunter fetch' to list the 35 conditions")
    return Y[arg]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="mutation-hunter", description=__doc__)
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("fetch", help="download and verify the public data; print the 35 conditions")
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-o", "--out", required=True, help="output prefix (writes .tsv, .json, .html)")
    common.add_argument("--max-freq", type=float, default=0.02, help="rare = allele frequency <= this (default 0.02)")
    common.add_argument("--min-carriers", type=int, default=3, help="minimum carriers to test (default 3)")
    ph = sub.add_parser("phewas", parents=[common], help="one variant or gene vs. all 35 phenotypes")
    ph.add_argument("query", help="GENE (burden of rare variants) or GENE:REFposALT, e.g. YNL079C:A123V")
    gw = sub.add_parser("gwas", parents=[common], help="one phenotype vs. all genes")
    gw.add_argument("phenotype", help="one of the 35 conditions, or a TSV file strain<TAB>value")
    gw.add_argument("--level", choices=["gene", "variant"], default="gene")
    fd = sub.add_parser("fdr", help="merge result TSVs and recompute q over all tests")
    fd.add_argument("tsv", nargs="+")
    fd.add_argument("-o", "--out", required=True, help="output prefix")
    a = p.parse_args(argv)

    try:
        if a.command == "fetch":
            for name in data.FILES:
                data.fetch(name)
            data.proteome_dir()
            print("\n".join(data.load_phenotypes().columns))
        elif a.command == "phewas":
            res, rare = analysis.phewas(a.query, a.max_freq, a.min_carriers)
            report.write(a.out, res, vars(a), rare)
        elif a.command == "gwas":
            res = analysis.gwas(_phenotype(a.phenotype), a.level, a.max_freq, a.min_carriers)
            report.write(a.out, res, vars(a))
        else:
            res = pd.concat([pd.read_csv(f, sep="\t") for f in a.tsv], ignore_index=True)
            res["q"] = stats.bh(res["p"])
            res.to_csv(f"{a.out}.tsv", sep="\t", index=False)
    except (ValueError, KeyError, RuntimeError) as e:
        print(f"mutation-hunter: error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
