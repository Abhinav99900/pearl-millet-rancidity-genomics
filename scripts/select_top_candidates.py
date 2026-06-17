#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def clean_gene(x):
    return str(x).strip()


def clean_impact(row):
    impact_text = str(row.get("impact_category_for_top_candidate_selection", row.get("impact", ""))).lower()
    if "lof" in impact_text or "loss" in impact_text:
        return "LoF-associated"
    if "protein" in impact_text or "moderate" in impact_text:
        return "Protein-altering"
    if "regulatory" in impact_text or "modifier" in impact_text:
        return "Regulatory"
    if "low" in impact_text:
        return "LOW-impact coding"
    return str(row.get("impact_category_for_top_candidate_selection", "Unclassified"))


def impact_score(row):
    if "impact_score" in row and pd.notna(row["impact_score"]):
        try:
            return int(float(row["impact_score"]))
        except Exception:
            pass

    c = clean_impact(row)
    return {
        "LoF-associated": 4,
        "Protein-altering": 3,
        "Regulatory": 2,
        "LOW-impact coding": 1
    }.get(c, 0)


def relevance_score_and_group(function_text):
    t = str(function_text).lower()

    direct = [
        "lipase", "esterase", "gdsl", "phospholipase", "triacylglycerol",
        "tag lipase", "lipoxygenase", "lox", "carboxylesterase", "thioesterase"
    ]

    oxidative = [
        "peroxidase", "aldehyde dehydrogenase", "alcohol dehydrogenase",
        "superoxide dismutase", "catalase", "glutathione", "oxidoreductase",
        "peroxide", "detox", "ros", "oxidative", "redox", "dehydrogenase",
        "aminobutyraldehyde", "betaine aldehyde"
    ]

    lipid = [
        "fatty acid", "fatty-acid", "acyl-coa", "beta-oxidation",
        "desaturase", "reductase", "lipid metabolism", "glycerolipid",
        "membrane lipid", "hydratase", "thiolase"
    ]

    if any(x in t for x in direct):
        return 4, "Direct rancidity initiation"
    if any(x in t for x in oxidative):
        return 3, "Oxidative deterioration/detoxification"
    if any(x in t for x in lipid):
        return 2, "Supportive lipid metabolism"
    if any(x in t for x in ["stress", "metabolism", "transferase", "response", "enzyme"]):
        return 1, "General stress/metabolism"
    return 0, "Unclear/unknown"


def annotation_support(row):
    text = " ".join([str(row.get("source", "")), str(row.get("function", ""))]).lower()

    swiss = "swiss" in text or "swissprot" in text or "swiss-prot" in text
    trembl = "trembl" in text
    eggnog = "eggnog" in text
    domain = "pfam" in text or "cdd" in text or "domain" in text
    pathway = "kegg" in text or "ko" in text or "go" in text or "ec" in text

    if swiss and (domain or pathway):
        return "High"
    if swiss or (trembl and (domain or pathway)) or (eggnog and (domain or pathway)):
        return "Moderate"
    if trembl or eggnog or domain:
        return "Low"
    return "Limited"


def support_rank(x):
    return {"High": 3, "Moderate": 2, "Low": 1, "Limited": 0}.get(x, 0)


def prepare(df, label):
    out = df.copy()
    out["gene_id_clean"] = out["gene_id"].apply(clean_gene)
    out["candidate_gene_set"] = label
    out["Impact class"] = out.apply(clean_impact, axis=1)
    out["Impact score"] = out.apply(impact_score, axis=1)

    rel = out["function"].fillna("").apply(relevance_score_and_group)
    out["Rancidity relevance score"] = [x[0] for x in rel]
    out["Functional group"] = [x[1] for x in rel]

    out["Priority score"] = 2 * out["Impact score"] + 2 * out["Rancidity relevance score"]
    out["Annotation support"] = out.apply(annotation_support, axis=1)
    out["support_rank"] = out["Annotation support"].apply(support_rank)
    out["Putative function"] = out["function"].fillna("")

    return out.drop_duplicates("gene_id_clean")


def diversity_select(df, n=5):
    df = df.sort_values(
        ["Priority score", "Impact score", "Rancidity relevance score", "support_rank"],
        ascending=False
    ).reset_index(drop=True)

    selected = []
    impact_count = {}
    group_count = {}

    impact_cap = {
        "LoF-associated": 2,
        "Protein-altering": 2,
        "Regulatory": 1,
        "LOW-impact coding": 1,
        "Unclassified": 1
    }

    group_cap = 3

    for _, row in df.iterrows():
        impact = row["Impact class"]
        group = row["Functional group"]

        if impact_count.get(impact, 0) >= impact_cap.get(impact, 1):
            continue
        if group_count.get(group, 0) >= group_cap:
            continue

        selected.append(row)
        impact_count[impact] = impact_count.get(impact, 0) + 1
        group_count[group] = group_count.get(group, 0) + 1

        if len(selected) == n:
            break

    if len(selected) < n:
        used = set([x["gene_id_clean"] for x in selected])
        for _, row in df.iterrows():
            if row["gene_id_clean"] not in used:
                selected.append(row)
                used.add(row["gene_id_clean"])
            if len(selected) == n:
                break

    return pd.DataFrame(selected)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--g732", required=True)
    p.add_argument("--g538", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--n", type=int, default=5)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df732 = prepare(read_tsv(args.g732), "GHB-732(LR)-unique")
    df538 = prepare(read_tsv(args.g538), "GHB-538(HR)-unique")

    genes732 = set(df732["gene_id_clean"])
    genes538 = set(df538["gene_id_clean"])

    common_genes = genes732 & genes538
    unique732_genes = genes732 - genes538
    unique538_genes = genes538 - genes732

    u732 = df732[df732["gene_id_clean"].isin(unique732_genes)].copy()
    u538 = df538[df538["gene_id_clean"].isin(unique538_genes)].copy()

    c732 = df732[df732["gene_id_clean"].isin(common_genes)].copy()
    c538 = df538[df538["gene_id_clean"].isin(common_genes)].copy()

    merged = c732.merge(c538, on="gene_id_clean", suffixes=("_732", "_538"))

    common_rows = []
    for _, r in merged.iterrows():
        # Conservative rule: lower impact score for same common gene
        if r["Impact score_732"] <= r["Impact score_538"]:
            impact_class = r["Impact class_732"]
            impact_score_final = r["Impact score_732"]
        else:
            impact_class = r["Impact class_538"]
            impact_score_final = r["Impact score_538"]

        # Use stronger functional relevance if annotation differs
        if r["Rancidity relevance score_732"] >= r["Rancidity relevance score_538"]:
            func_group = r["Functional group_732"]
            rel_score = r["Rancidity relevance score_732"]
            put_func = r["Putative function_732"]
        else:
            func_group = r["Functional group_538"]
            rel_score = r["Rancidity relevance score_538"]
            put_func = r["Putative function_538"]

        support = r["Annotation support_732"] if r["support_rank_732"] >= r["support_rank_538"] else r["Annotation support_538"]
        support_value = max(r["support_rank_732"], r["support_rank_538"])

        common_rows.append({
            "candidate_gene_set": "Common",
            "gene_id_clean": r["gene_id_clean"],
            "Impact class": impact_class,
            "Impact score": impact_score_final,
            "Functional group": func_group,
            "Rancidity relevance score": rel_score,
            "Priority score": 2 * impact_score_final + 2 * rel_score,
            "Annotation support": support,
            "support_rank": support_value,
            "Putative function": put_func
        })

    common = pd.DataFrame(common_rows)

    top_common = diversity_select(common, args.n)
    top_732 = diversity_select(u732, args.n)
    top_538 = diversity_select(u538, args.n)

    final = pd.concat([top_common, top_732, top_538], ignore_index=True)

    final_table = final[[
        "candidate_gene_set",
        "gene_id_clean",
        "Impact class",
        "Functional group",
        "Putative function",
        "Annotation support",
        "Priority score"
    ]].rename(columns={
        "candidate_gene_set": "Candidate gene set",
        "gene_id_clean": "Gene ID"
    })

    common.to_csv(outdir / "COMMON_ALL_SCORED.tsv", sep="\t", index=False)
    u732.to_csv(outdir / "GHB732_UNIQUE_ALL_SCORED.tsv", sep="\t", index=False)
    u538.to_csv(outdir / "GHB538_UNIQUE_ALL_SCORED.tsv", sep="\t", index=False)
    final_table.to_csv(outdir / "TABLE6_READY_TOP_CANDIDATES.tsv", sep="\t", index=False)

    try:
        final_table.to_excel(outdir / "TABLE6_READY_TOP_CANDIDATES.xlsx", index=False)
    except Exception:
        pass

    with open(outdir / "TOP_CANDIDATE_SELECTION_LOG.txt", "w") as f:
        f.write(f"GHB-732 final candidates: {len(genes732)}\n")
        f.write(f"GHB-538 final candidates: {len(genes538)}\n")
        f.write(f"Common candidates: {len(common_genes)}\n")
        f.write(f"GHB-732 unique candidates: {len(unique732_genes)}\n")
        f.write(f"GHB-538 unique candidates: {len(unique538_genes)}\n")
        f.write("Formula: Priority score = 2 x impact score + 2 x rancidity relevance score\n")
        f.write("Annotation support was reported qualitatively and was not added to the score.\n")
        f.write("Diversity-aware selection was applied to reduce dominance of one impact class or functional group.\n")

    print("Done")
    print(f"Output: {outdir}")
    print(f"Common={len(common_genes)} | 732-unique={len(unique732_genes)} | 538-unique={len(unique538_genes)}")


if __name__ == "__main__":
    main()
