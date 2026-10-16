# Build check against the public proteome

`compare_proteome.py` compares the missense carriers of a built panel with `1011proteome.tar.gz`
(SNP-only proteins). Run on panel draft built 2026-10-15 and reproduced exactly on panel-0.1.0 (pipeline commit 62c724b, 2026-10-16;
bcftools 1.20, min GQ 20, min DP 5):

| chromosome1, 9,972 single-codon missense rows | agreement |
|---|---|
| homozygous carriers the proteome agrees with | 204,428 / 217,364 = 0.940 |
| heterozygous carriers the proteome agrees with | 50,688 / 54,859 = 0.924 |
| proteome carriers the panel also calls | 255,116 / 258,466 = 0.987 |

Of the 12,936 homozygous disagreements, 9,855 (76%) are codons with a second change (the proteome shows a third
amino acid; `csq -l` annotates each variant alone, a documented limit) and 3,081 (24%) show the reference amino acid
in the proteome. Heterozygous calls agree far more than the 50% a random allele draw would give, so the proteome's
"random allele" is not uniform; we did not investigate further.
