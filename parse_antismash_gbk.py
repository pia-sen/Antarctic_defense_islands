#!/usr/bin/env python3
"""
parse_antismash_gbk.py
Parse antiSMASH GenBank (.gbk) files.

Tasks:
- Extracts BGC class/subclass
- Extracts sub-clusters/domains per gene
- Flags incomplete clusters
- Recursive parsing of folders
- Outputs single CSV
"""

import sys
import csv
from pathlib import Path
from Bio import SeqIO

def extract_domains(qualifiers):
    # antiSMASH domain key qualifers (multiple)
    domain_keys = ["Note", "note", "antiSMASH_subclass", "product"]
    domains = []
    for key in domain_keys:
        if key in qualifiers:
            domains += qualifiers[key]
    return "; ".join(domains)

def cds_overlaps_cluster(cds, cluster):
    return not (int(cds.location.end) < cluster.location.start or int(cds.location.start) > cluster.location.end)

def parse_gbk_file(gbk_path, writer):
    for record in SeqIO.parse(gbk_path, "genbank"):
        cluster_features = [f for f in record.features if f.type in ("region", "cluster")]

        for cluster in cluster_features:
            qual = cluster.qualifiers
            bgc_id = qual.get("region_name", qual.get("note", [""]))[0]
            bgc_type = qual.get("product", [""])[0]
            bgc_subclass = qual.get("subclass", [""])[0] if "subclass" in qual else ""
            bgc_start = int(cluster.location.start) + 1
            bgc_end = int(cluster.location.end)
            # Flag incomplete clusters if at contig edges or noted as incomplete
            cluster_complete = "Yes"
            if int(cluster.location.start) == 0 or "incomplete" in str(qual).lower():
                cluster_complete = "No"

            for feature in record.features:
                if feature.type == "CDS" and cds_overlaps_cluster(feature, cluster):
                    qualifiers = feature.qualifiers
                    gene_id = qualifiers.get("locus_tag", [""])[0]
                    product = qualifiers.get("product", [""])[0]
                    domains = extract_domains(qualifiers)
                    gene_start = int(feature.location.start) + 1
                    gene_end = int(feature.location.end)
                    strand = feature.location.strand

                    writer.writerow([
                        gbk_path.name, bgc_id, bgc_type, bgc_subclass, bgc_start, bgc_end, cluster_complete,
                        gene_id, product, gene_start, gene_end, strand, domains
                    ])

def parse_folder(input_folder, output_file):
    gbk_files = list(Path(input_folder).rglob("*.gbk"))
    if not gbk_files:
        print(f"[!] No .gbk files found in {input_folder}")
        return

    with open(output_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "Filename", "BGC_ID", "BGC_Type", "BGC_Subclass", "BGC_Start", "BGC_End",
            "Cluster_Complete", "Gene_ID", "Gene_Product", "Gene_Start", "Gene_End", "Strand", "Domains"
        ])

        for gbk in gbk_files:
            print(f"[+] Parsing {gbk}")
            parse_gbk_file(gbk, writer)

    print(f"[✓] Fancy antiSMASH data saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_antismash_gbk.py input_folder output_file.csv")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_file = sys.argv[2]
    parse_folder(input_folder, output_file)
```
Loop parse_antismash_gbk.py:   \
```{bash}
#pwd:/Users/pia.sen/Documents/Research/ATA_defense/BV_partial_output
while read -r j; do
python3 parse_antismash_gbk.py ./"${j}" "${j}_antismash.csv"
done < list.txt
