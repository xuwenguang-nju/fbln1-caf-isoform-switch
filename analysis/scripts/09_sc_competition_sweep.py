"""Step 9 (this project): competition sweep for a single-cell splicing study centred
on the oral-cancer microenvironment.

Counts and retrieves PubMed records for four blocks of queries:
  A. single-cell splicing as a field (global volume + where it publishes)
  B. oral / head-and-neck single-cell literature (how crowded the tissue is)
  C. the specific niche: CAF, macrophage, OPMD/OSMF + splicing/isoform
  D. closest prior art: any oral-cancer paper already doing cell-type splicing

Inputs : none (live PubMed E-utilities)
Outputs: data/processed/sc_competition_counts.json
         data/processed/sc_competition_records.json
         results/tables/sc_competition_sweep.csv
"""
import json, os, time, csv, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT = "data/processed"
os.makedirs(OUT, exist_ok=True)

GROUPS = {
    "A_field": {
        "sc_as_all": 'alternative splicing[TIAB] AND ("single-cell"[TIAB] OR "single cell"[TIAB] OR scRNA-seq[TIAB] OR "single-nucleus"[TIAB])',
        "sc_as_cancer": 'alternative splicing[TIAB] AND ("single-cell"[TIAB] OR "single cell"[TIAB] OR scRNA-seq[TIAB]) AND (cancer[TIAB] OR tumor[TIAB] OR tumour[TIAB] OR carcinoma[TIAB])',
        "isoform_sc_cancer": '(isoform*[TIAB]) AND ("single-cell"[TIAB] OR scRNA-seq[TIAB]) AND (cancer[TIAB] OR tumor[TIAB] OR carcinoma[TIAB])',
        "celltype_specific_splicing": '("cell type-specific"[TIAB] OR "cell-type-specific"[TIAB] OR "cell type specific"[TIAB]) AND (splicing[TIAB] OR isoform*[TIAB])',
    },
    "B_tissue": {
        "hnsc_sc": '("head and neck"[TIAB] OR HNSCC[TIAB]) AND ("single-cell"[TIAB] OR "single cell"[TIAB] OR scRNA-seq[TIAB])',
        "oscc_sc": '(OSCC[TIAB] OR "oral squamous cell carcinoma"[TIAB] OR "oral cancer"[TIAB] OR "oral cavity"[TIAB]) AND ("single-cell"[TIAB] OR "single cell"[TIAB] OR scRNA-seq[TIAB])',
        "oscc_sc_recent": '(OSCC[TIAB] OR "oral squamous cell carcinoma"[TIAB] OR "oral cancer"[TIAB]) AND ("single-cell"[TIAB] OR scRNA-seq[TIAB]) AND 2023:2026[dp]',
        "hnsc_sc_recent": '("head and neck"[TIAB] OR HNSCC[TIAB]) AND ("single-cell"[TIAB] OR scRNA-seq[TIAB]) AND 2024:2026[dp]',
    },
    "C_niche": {
        "oscc_sc_splicing": '(OSCC[TIAB] OR "oral squamous"[TIAB] OR "oral cancer"[TIAB] OR "head and neck"[TIAB] OR HNSCC[TIAB]) AND (splicing[TIAB] OR isoform*[TIAB]) AND ("single-cell"[TIAB] OR scRNA-seq[TIAB] OR "single cell"[TIAB])',
        "caf_splicing": '("cancer-associated fibroblast"[TIAB] OR "cancer associated fibroblast"[TIAB] OR myofibroblast[TIAB]) AND (splicing[TIAB] OR isoform*[TIAB]) AND (cancer[TIAB] OR tumor[TIAB] OR carcinoma[TIAB])',
        "caf_sc_splicing": '("cancer-associated fibroblast"[TIAB] OR CAF[TIAB]) AND ("single-cell"[TIAB] OR scRNA-seq[TIAB]) AND (splicing[TIAB] OR isoform*[TIAB])',
        "mac_sc_splicing": '(macrophage[TIAB]) AND ("single-cell"[TIAB] OR scRNA-seq[TIAB]) AND (splicing[TIAB] OR isoform*[TIAB]) AND (cancer[TIAB] OR tumor[TIAB])',
        "opmd_sc": '("oral leukoplakia"[TIAB] OR "oral potentially malignant"[TIAB] OR "oral submucous fibrosis"[TIAB] OR erythroplakia[TIAB] OR "oral precancer"[TIAB]) AND ("single-cell"[TIAB] OR "single cell"[TIAB] OR scRNA-seq[TIAB])',
        "osmf": '"oral submucous fibrosis"[TIAB] AND ("single-cell"[TIAB] OR scRNA-seq[TIAB] OR fibroblast[TIAB])',
        "opmd_splicing": '("oral leukoplakia"[TIAB] OR "oral submucous fibrosis"[TIAB] OR "oral potentially malignant"[TIAB]) AND (splicing[TIAB] OR isoform*[TIAB])',
        "hnsc_immunotherapy_sc": '("head and neck"[TIAB] OR HNSCC[TIAB] OR "oral"[TIAB]) AND (immunotherapy[TIAB] OR "immune checkpoint"[TIAB] OR pembrolizumab[TIAB] OR nivolumab[TIAB]) AND ("single-cell"[TIAB] OR scRNA-seq[TIAB])',
    },
    "D_prior_art": {
        "sc_as_hnsc_or_oral": '(splicing[TIAB] OR isoform*[TIAB]) AND ("single-cell"[TIAB] OR scRNA-seq[TIAB]) AND ("head and neck"[TIAB] OR oral[TIAB] OR HNSCC[TIAB] OR OSCC[TIAB])',
        "isoform_switch_sc": '("isoform switch"[TIAB] OR "isoform switching"[TIAB] OR "isoform usage"[TIAB]) AND ("single-cell"[TIAB] OR scRNA-seq[TIAB])',
        "tme_splicing": '(splicing[TIAB] OR isoform*[TIAB]) AND ("tumor microenvironment"[TIAB] OR "tumour microenvironment"[TIAB] OR "immune microenvironment"[TIAB]) AND ("single-cell"[TIAB] OR scRNA-seq[TIAB])',
    },
}


def eutil(endpoint, params, tries=3):
    params = dict(params)
    params.setdefault("tool", "wisp-science")
    params.setdefault("email", "researcher@example.org")
    url = BASE + endpoint + "?" + urllib.parse.urlencode(params)
    err = None
    for _ in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(2)
    raise err


def search(term, retmax=40, sort="relevance"):
    r = json.loads(eutil("esearch.fcgi", {"db": "pubmed", "term": term, "retmax": retmax,
                                          "retmode": "json", "sort": sort}))
    return int(r["esearchresult"]["count"]), r["esearchresult"]["idlist"]


def fetch(pmids):
    if not pmids:
        return []
    xml = eutil("efetch.fcgi", {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"})
    root = ET.fromstring(xml)
    recs = []
    for art in root.findall(".//PubmedArticle"):
        doi = ""
        for aid in art.findall("./PubmedData/ArticleIdList/ArticleId"):
            if aid.get("IdType") == "doi":
                doi = (aid.text or "").strip().lower()
        ab = art.find(".//Abstract")
        recs.append({
            "pmid": art.findtext(".//MedlineCitation/PMID") or "",
            "doi": doi,
            "title": " ".join((art.findtext(".//ArticleTitle") or "").split()),
            "journal": art.findtext(".//Journal/Title") or "",
            "year": (art.findtext(".//JournalIssue/PubDate/Year")
                     or (art.findtext(".//JournalIssue/PubDate/MedlineDate") or "")[:4]),
            "abstract": " ".join("".join(ab.itertext()).split()) if ab is not None else "",
        })
    return recs


def main(retmax=40):
    counts, recs_all = {}, {}
    for group, qs in GROUPS.items():
        for key, term in qs.items():
            n, ids = search(term, retmax)
            recs = fetch(ids)
            for r in recs:
                r["query_key"] = key
                r["group"] = group
            counts[key] = {"group": group, "count": n, "query": term, "fetched": len(recs)}
            recs_all[key] = recs
            print(f"{group:12s} {key:28s} n={n:5d} fetched={len(recs)}")
            time.sleep(0.4)

    json.dump(counts, open(f"{OUT}/sc_competition_counts.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    flat = [r for v in recs_all.values() for r in v]
    json.dump(flat, open(f"{OUT}/sc_competition_records.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    with open("results/tables/sc_competition_sweep.csv", "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["group", "query_key", "pubmed_count", "query"])
        for k, v in counts.items():
            w.writerow([v["group"], k, v["count"], v["query"]])
    print("records:", len(flat))


if __name__ == "__main__":
    main()
