"""Manual check of a built panel against the public SNP-only proteome (1011proteome.tar.gz, extracted).

For every missense row on one chromosome, compare its carriers with the strains whose proteome sequence has the
alternative amino acid at that position. Disagreements come mostly from codons with two changes: bcftools csq -l
annotates each variant alone, while the proteome applies both (see README.md in this folder).

usage: python compare_proteome.py <panel folder> <Proteome_1011 folder> [chromosome1]
"""

import sys
from pathlib import Path

import numpy as np

from mutation_hunter.store import Panel


def load_gene(folder: Path, gene: str) -> dict[str, str]:
    seqs, strain = {}, None
    for line in (folder / f"{gene}.fasta").read_text().splitlines():
        if line.startswith(">"):
            strain = line[1:].split(f"_{gene}_")[0]
            seqs[strain] = ""
        else:
            seqs[strain] += line.strip()
    return seqs


def main(panel_dir, proteome_dir, chrom="chromosome1") -> None:
    p, prot = Panel(panel_dir), Path(proteome_dir)
    strains = np.array(p.strains.index)
    hom = het = n_rows = 0
    hom_ok = het_ok = prot_found = prot_total = 0
    for gene in p.genes.index[p.genes["chrom"] == chrom]:
        if not (prot / f"{gene}.fasta").exists():
            continue
        var, geno = p.gene(gene)
        seqs = load_gene(prot, gene)
        for i, (cls, aa) in enumerate(zip(var["class"], var["aa_change"])):
            ch = aa.split(">")
            if cls != "missense" or len(ch) != 2 or not ch[1][:-1].isdigit():
                continue
            pos, alt = int(ch[1][:-1]), ch[1][-1]
            in_prot = np.array([s in seqs for s in strains])
            prot_alt = np.array([s in seqs and len(seqs[s]) >= pos and seqs[s][pos - 1] == alt for s in strains])
            g = geno[i]
            n_rows += 1
            hom += ((g == 2) & in_prot).sum()
            hom_ok += ((g == 2) & prot_alt).sum()
            het += ((g == 1) & in_prot).sum()
            het_ok += ((g == 1) & prot_alt).sum()
            prot_total += (prot_alt & (g >= 0)).sum()
            prot_found += (prot_alt & (g > 0)).sum()
    print(f"{chrom}: {n_rows} single-codon missense rows")
    print(f"homozygous carriers that the proteome agrees with: {hom_ok}/{hom} = {hom_ok / hom:.3f}")
    print(f"heterozygous carriers that the proteome agrees with: {het_ok}/{het} = {het_ok / het:.3f}")
    print(f"proteome carriers (called in panel) that the panel also calls: {prot_found}/{prot_total} = {prot_found / prot_total:.3f}")


if __name__ == "__main__":
    main(*sys.argv[1:])
