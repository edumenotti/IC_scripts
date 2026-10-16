#!/usr/bin/env bash
# Normalise and annotate one region of the panel VCF. Same flags for the panel and for external strains.
# usage: annotate.sh <vcf.gz> <region> <ref.fa> <ref.gff3> <out.bcf>
set -euo pipefail
bcftools view -r "$2" "$1" -Ou \
  | bcftools norm -f "$3" -c x -m -any -a --atom-overlaps . --old-rec-tag OLD_REC -Ou 2>/dev/null \
  | bcftools csq -f "$3" -g "$4" -l -s - -v 0 -Ob -o "$5"
bcftools index -f "$5"
