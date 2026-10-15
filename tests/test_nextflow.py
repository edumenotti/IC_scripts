import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from mutation_hunter import data

MAIN = Path(__file__).parents[1] / "nextflow" / "main.nf"


@pytest.mark.skipif(not shutil.which("nextflow"), reason="nextflow not installed")
def test_pipeline_runs_and_merges(panel, tmp_path):
    for name in data.FILES:  # warm the cache here: the subprocesses do not see the patched MD5s
        data.fetch(name)
    env = os.environ | {"PATH": f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}"}
    subprocess.run(["nextflow", "run", str(MAIN), "--outdir", str(tmp_path / "res"), "--max_freq", "0.1"],
                   cwd=tmp_path, env=env, check=True, capture_output=True)
    merged = pd.read_csv(tmp_path / "res" / "gwas_all.tsv", sep="\t")
    assert set(merged["phenotype"]) == {"COND1", "COND2"} and merged["q"].notna().any()
