# Utility script: Convert GenBank annotations to BED12 format
#
# Purpose:
#   Export annotated genomic features 
#   from a GenBank file into BED12 format for visualization and
#   coordinate-based comparison with RNA-seq coverage and ANNOgesic results.
#
import sys
from Bio import SeqIO

gbk, out = sys.argv[1], sys.argv[2]
rec = SeqIO.read(gbk, "genbank")
chrom = rec.id

def pick_name(q, fallback):
    for k in ("locus_tag", "gene", "product"):
        if k in q and q[k]:
            return q[k][0].replace(" ", "_")
    return fallback

with open(out, "w") as f:
    for feat in rec.features:
        if feat.type not in ("CDS", "rRNA", "tRNA", "ncRNA", "gene"):
            continue

        strand = "."
        if feat.location.strand == 1:
            strand = "+"
        elif feat.location.strand == -1:
            strand = "-"

        parts = list(getattr(feat.location, "parts", [feat.location]))
        blocks = [(int(p.start), int(p.end)) for p in parts]
        if not blocks:
            continue

        start0 = min(s for s, e in blocks)
        end = max(e for s, e in blocks)

        block_sizes = ",".join(str(e - s) for s, e in blocks) + ","
        block_starts = ",".join(str(s - start0) for s, e in blocks) + ","

        name = pick_name(feat.qualifiers, feat.type)

        f.write("\t".join(map(str, [
            chrom, start0, end, name, 0, strand,
            start0, end, 0,
            len(blocks), block_sizes, block_starts
        ])) + "\n")

print(f"Wrote BED12: {out}")
