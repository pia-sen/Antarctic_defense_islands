import pandas as pd
import glob
import os

MAX_DIST_CONTIG = 10000  # 10 kb distance
EDGE_DIST = 10000        # 10 kb distance to contig end
all_pairs = []

# Process all FAA files (one per MAG)
for faa in glob.glob("*.faa"):
    prefix = faa.replace(".faa", "")
    print("Processing MAG:", prefix)

# Parse FAA headers for coordinates
    genes = []
    contig_lengths = {}
    with open(faa) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                gene_id = header.split("#")[0].strip()
                parts = header.split("#")
                start = int(parts[1].strip())
                end = int(parts[2].strip())
                strand = "+" if parts[3].strip() == "1" else "-"
                contig = "_".join(gene_id.split("_")[:-1])
                genes.append({
                    "gene_id": gene_id,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand
                })
                contig_lengths[contig] = max(contig_lengths.get(contig, 0), end)

    gene_df = pd.DataFrame(genes)

# Load all DIAMOND output toxin and antitoxin hit files

    cols = ["gene_id","ta_id","description","pident","length",
            "qstart","qend","sstart","send","evalue","bitscore"]

    hits_list = []

    for file in glob.glob(f"{prefix}_vs_*_exp_db_hits.tsv"):
        if os.path.getsize(file) == 0:
            continue
        base = os.path.basename(file)

        # Extract Type and role from filename
        parts = base.split("_vs_")[1].split("_")
        type_label = parts[0] + "_" + parts[1]   # e.g., type_II
        role = parts[2]                           # "T" or "AT"

        df = pd.read_csv(file, sep="\t", header=None, names=cols)
        df["type_label"] = type_label
        df["role"] = "toxin" if role.upper() == "T" else "antitoxin"

        hits_list.append(df)

    if not hits_list:
        continue

    hits = pd.concat(hits_list, ignore_index=True)
    hits = hits.sort_values("bitscore", ascending=False).drop_duplicates(["gene_id","role"])
    hits = hits.merge(gene_df, on="gene_id", how="left").dropna(subset=["contig"])

    toxins = hits[hits["role"]=="toxin"]
    antitoxins = hits[hits["role"]=="antitoxin"]

# Find TA pairs by Type
    for _, t in toxins.iterrows():
        same_type_ants = antitoxins[antitoxins["type_label"]==t["type_label"]]
        for _, a in same_type_ants.iterrows():
            distance = None
            pair_ok = False

            if t["contig"] == a["contig"]:
                distance = min(abs(t["start"] - a["end"]), abs(t["end"] - a["start"]))
                if distance <= MAX_DIST_CONTIG:
                    pair_ok = True
            else:
                # check contig edge for cross-contig
                t_near_edge = t["end"] >= contig_lengths[t["contig"]] - EDGE_DIST
                a_near_edge = a["end"] >= contig_lengths[a["contig"]] - EDGE_DIST
                if t_near_edge or a_near_edge:
                    pair_ok = True

            if pair_ok:
                all_pairs.append({
                    "MAG": prefix,
                    "type": t["type_label"],
                    "toxin_gene": t["gene_id"],
                    "antitoxin_gene": a["gene_id"],
                    "distance": distance if distance is not None else "NA",
                    "toxin_contig": t["contig"],
                    "antitoxin_contig": a["contig"],
                    "toxin_hit": t["ta_id"],
                    "antitoxin_hit": a["ta_id"]
                })

# Save results
pairs_df = pd.DataFrame(all_pairs)
pairs_df.to_csv("TA_pairs_by_type_10kb.tsv", sep="\t", index=False)
print("Total TA pairs found:", len(all_pairs))
