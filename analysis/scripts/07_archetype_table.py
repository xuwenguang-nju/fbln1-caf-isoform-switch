"""Step 6: assemble the archetype overview table from the retrieved evidence.

Inputs : data/processed/as_pubmed_query_counts.json
         data/processed/as_records_enriched.csv
         data/processed/archetype_component_prevalence.csv
         data/processed/as_year_trend.json
Outputs: results/tables/as_archetype_overview.csv
"""
import json
import pandas as pd

counts = json.load(open("data/processed/as_pubmed_query_counts.json", encoding="utf-8"))
recs = pd.read_csv("data/processed/as_records_enriched.csv")
trend = json.load(open("data/processed/as_year_trend.json", encoding="utf-8"))

rows = []
for key, meta in counts.items():
    sub = recs[recs["archetype_query"] == key]
    rows.append({
        "archetype": key,
        "pubmed_hits_2026_02": meta["count"],
        "sampled_records": len(sub),
        "median_citations_sample": round(sub["cited_by"].median(), 1) if sub["cited_by"].notna().any() else None,
        "max_citations_sample": sub["cited_by"].max(),
        "median_year_sample": round(pd.to_numeric(sub["year"], errors="coerce").median(), 0),
        "query": meta["query"],
    })
out = pd.DataFrame(rows).sort_values("pubmed_hits_2026_02", ascending=False)
out.to_csv("results/tables/as_archetype_overview.csv", index=False, encoding="utf-8-sig")
print(out.drop(columns=["query"]).to_string(index=False))
print("\nyear trend keys:", list(trend))
