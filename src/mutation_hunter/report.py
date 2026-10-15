"""Write results as TSV (the product), JSON (provenance) and a single-file HTML report."""

import base64
import html
import io
import json
import sys
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy  # noqa: E402
from scipy import stats  # noqa: E402

from . import __version__, data  # noqa: E402


def lambda_gc(p) -> float:
    """Genomic inflation: median observed chi2 / expected (0.455). ~1 is calibrated; >>1 means confounding."""
    p = np.asarray(p, float)
    p = p[~np.isnan(p)]
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)) if p.size else float("nan")


def provenance(params: dict) -> dict:
    return {
        "mutation_hunter": __version__,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": params,
        "data_source": data.BASE_URL,
        "data_md5": data.FILES,
        "python": sys.version.split()[0],
        "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
    }


def _png(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return f'<img src="data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}">'


def _forest(res: pd.DataFrame) -> str:
    r = res.dropna(subset=["beta"]).sort_values("beta")
    fig, ax = plt.subplots(figsize=(6, 0.25 * len(r) + 1))
    ax.errorbar(r["beta"], r["phenotype"], xerr=1.96 * r["se"], fmt="o")
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlabel("effect on growth ratio (beta ± 95% CI)")
    return _png(fig)


def _manhattan(res: pd.DataFrame) -> str:
    r = res.dropna(subset=["p"]).sort_values("gene")
    fig, ax = plt.subplots(figsize=(9, 3))
    ax.scatter(range(len(r)), -np.log10(r["p"]), s=4)
    ax.set_xlabel("tests, ordered by systematic gene name")
    ax.set_ylabel("-log10 p")
    return _png(fig)


def _lollipop(rv: pd.DataFrame) -> str:
    fig, ax = plt.subplots(figsize=(9, 2.5))
    ax.stem(rv["pos"], rv["n_carriers"])
    ax.set_xlabel("protein position")
    ax.set_ylabel("carriers")
    return _png(fig)


def write(prefix: str, res: pd.DataFrame, params: dict, rare: pd.DataFrame | None = None) -> None:
    """<prefix>.tsv, <prefix>.json and <prefix>.html. rare is given for PheWAS, None for GWAS."""
    res.to_csv(f"{prefix}.tsv", sep="\t", index=False)
    prov = provenance(params) | {"lambda_gc": lambda_gc(res["p"])}
    with open(f"{prefix}.json", "w") as f:
        json.dump(prov, f, indent=2)
    figs = [_forest(res)] if rare is not None else [_manhattan(res)]
    if rare is not None and not rare.empty:
        figs.append(_lollipop(rare))
    top = res.nsmallest(50, "p").to_html(index=False, float_format="%.3g")
    title = html.escape(f"mutation-hunter {params.get('command', '')}: {params.get('query', params.get('phenotype', ''))}")
    with open(f"{prefix}.html", "w") as f:
        f.write(f"<!doctype html><meta charset='utf-8'><title>{title}</title><h1>{title}</h1>"
                f"<p>Genomic inflation &lambda;<sub>GC</sub> = {prov['lambda_gc']:.2f} (~1 is calibrated)</p>"
                f"{''.join(figs)}<h2>Top tests</h2>{top}"
                f"<h2>Provenance</h2><pre>{html.escape(json.dumps(prov, indent=2))}</pre>")
