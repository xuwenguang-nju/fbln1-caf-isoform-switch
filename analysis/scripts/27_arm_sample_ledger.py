"""Stage 4 step 1: sample ledger for the three-arm oral-cavity CAF isoform study.

Answers one question with live API evidence: how many independent units
(patients / libraries / events) does each arm of the design actually have?

Queries (all live, no credentials):
  GDC API   - TCGA-HNSC splice-junction files + clinical fields -> subsite ledger
  GEO eutils - esummary for candidate scRNA / bulk resources   -> library counts,
               donor counts parsed from sample titles, platform chemistry

Outputs
  results/tables/tcga_hnsc_site_ledger.csv
  results/tables/tcga_hnsc_oral_power.csv
  results/tables/arm_sample_ledger.csv
  data/processed/arm_sample_ledger.json

Notes
  * Library/donor counts are re-derived on every run so the ledger cannot silently
    go stale; chemistry and PSI compatibility are curated annotations, each with an
    `evidence` string naming where the fact was read.
  * Every printed number in results/reports/arm_sample_ledger_and_framing.md comes
    from this script.
"""
import collections
import json
import re
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import norm

GDC = "https://api.gdc.cancer.gov/"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "Mozilla/5.0 (wisp-science; sample-ledger)"}

# TCGA-HNSC subsite grouping used throughout the study.
ORAL_STRICT = ["Other and unspecified parts of tongue", "Floor of mouth", "Gum",
               "Palate", "Lip", "Other and unspecified parts of mouth"]
ORAL_AMBIG = ["Other and ill-defined sites in lip, oral cavity and pharynx"]
OROPHARYNX = ["Tonsil", "Base of tongue", "Oropharynx"]
LARYNX = ["Larynx"]

# ----------------------------------------------------------------------------- GDC
def gdc(endpoint, payload=None):
    url = GDC + endpoint
    if payload is None:
        req = urllib.request.Request(url, headers=UA)
    else:
        req = urllib.request.Request(
            url, data=json.dumps(payload).encode(),
            headers={**UA, "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.loads(r.read().decode())


def fetch_tcga_hnsc():
    """One row per TCGA-HNSC splice-junction file (566 files, ~1.4 GB total)."""
    payload = {
        "filters": {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id",
                                     "value": ["TCGA-HNSC"]}},
            {"op": "in", "content": {"field": "data_type",
                                     "value": ["Splice Junction Quantification"]}}]},
        "fields": ("file_id,file_name,file_size,"
                   "cases.submitter_id,cases.primary_site,"
                   "cases.demographic.vital_status,cases.demographic.days_to_death,"
                   "cases.demographic.age_at_diagnosis,"
                   "cases.diagnoses.ajcc_pathologic_stage,"
                   "cases.diagnoses.days_to_last_follow_up,"
                   "cases.samples.sample_type"),
        "size": "2000", "format": "JSON"}
    hits = gdc("files", payload)["data"]["hits"]
    rows = []
    for h in hits:
        c = h["cases"][0]
        stype = {s["sample_type"] for s in c.get("samples", [])}
        dg = {}
        for d in (c.get("diagnoses") or []):
            for k, v in d.items():
                if v is not None and k not in dg:
                    dg[k] = v
        dem = c.get("demographic") or {}
        rows.append(dict(
            case=c["submitter_id"], site=c["primary_site"],
            sample_type="Primary Tumor" if "Primary Tumor" in stype
                        else ";".join(sorted(stype)),
            vital=dem.get("vital_status"), days_to_death=dem.get("days_to_death"),
            days_to_last_follow_up=dg.get("days_to_last_follow_up"),
            stage=dg.get("ajcc_pathologic_stage"),
            file_id=h["file_id"], file_mb=round(h["file_size"] / 1e6, 2)))
    return pd.DataFrame(rows)


def subsite_ledger(f, frames):
    """Per-subsite sample/event counts for primary-tumour junction files."""
    f = f.copy()
    f["os_time"] = f["days_to_death"].fillna(f["days_to_last_follow_up"])
    f["event"] = np.where(f["vital"] == "Dead", 1,
                          np.where(f["vital"] == "Alive", 0, np.nan))
    out = []
    for name, sites in frames:
        s = f[(f["sample_type"] == "Primary Tumor") & (f["site"].isin(sites))]
        ok = s["event"].notna() & s["os_time"].notna()
        out.append(dict(group=name, n_cases=s["case"].nunique(),
                        n_with_os=int(ok.sum()),
                        n_events=int(s.loc[ok, "event"].sum()),
                        median_os_or_fu_days=round(float(s.loc[ok, "os_time"].median()), 0)
                        if ok.any() else np.nan,
                        n_with_stage=int(s["stage"].notna().sum())))
    return pd.DataFrame(out)


def histology_check():
    """Two correctness checks a subsite analysis must declare.

    (1) Every TCGA-HNSC case carries a squamous diagnosis, but three carry a
        second, non-squamous diagnosis record (adenocarcinoma, lymphoma) -- an
        OSCC-only analysis should exclude them as a sensitivity step.
    (2) TCGA-HNSC has no finer resolution than primary_site: the tongue stratum
        ("Other and unspecified parts of tongue") cannot be split into oral
        tongue vs other tongue, so any "oral cavity" count derived from it is an
        upper bound, not a measurement.
    """
    payload = {"filters": {"op": "in", "content": {
                   "field": "cases.project.project_id", "value": ["TCGA-HNSC"]}},
               "fields": "submitter_id,primary_site,diagnoses.primary_diagnosis",
               "size": "1000", "format": "JSON"}
    hits = gdc("cases", payload)["data"]["hits"]
    rows = []
    for h in hits:
        dx = {d.get("primary_diagnosis") for d in h.get("diagnoses", [])
              if d.get("primary_diagnosis")}
        joined = ";".join(sorted(dx))
        rows.append(dict(case=h["submitter_id"], site=h["primary_site"],
                         has_scc="squamous" in joined.lower(),
                         mixed_histology=len(dx) > 1))
    d = pd.DataFrame(rows)
    oral_mask = d["site"].isin(ORAL_STRICT)
    return pd.DataFrame([
        dict(check="cases in TCGA-HNSC", n=len(d)),
        dict(check="cases with any squamous diagnosis", n=int(d.has_scc.sum())),
        dict(check="cases with >1 diagnosis record (mixed histology)",
             n=int(d.mixed_histology.sum())),
        dict(check="oral-cavity cases with mixed histology",
             n=int((d.mixed_histology & oral_mask).sum())),
        dict(check="cases whose subsite cannot be resolved finer than primary_site",
             n=int(oral_mask.sum())),
    ])


def site_breakdown(f):
    """Per-subsite counts for primary-tumour junction files, tagged by compartment."""
    def compartment(site):
        if site in ORAL_STRICT:
            return "oral cavity"
        if site in ORAL_AMBIG:
            return "lip/oral/pharynx, unspecified"
        if site in OROPHARYNX:
            return "oropharynx"
        if site in LARYNX:
            return "larynx"
        return "other"
    f = f.copy()
    f["os_time"] = f["days_to_death"].fillna(f["days_to_last_follow_up"])
    f["event"] = np.where(f["vital"] == "Dead", 1,
                          np.where(f["vital"] == "Alive", 0, np.nan))
    f = f[f["sample_type"] == "Primary Tumor"]
    out = []
    for site, s in f.groupby("site"):
        ok = s["event"].notna() & s["os_time"].notna()
        out.append(dict(site=site, compartment=compartment(site),
                        n_cases=s["case"].nunique(), n_with_os=int(ok.sum()),
                        n_deaths=int(s.loc[ok, "event"].sum()),
                        n_alive_followed=int((ok & (s["event"] == 0)).sum())))
    return (pd.DataFrame(out).sort_values("n_cases", ascending=False)
            .reset_index(drop=True))


def cox_power_table(n_cases, n_events):
    """Minimum detectable HR for a survival arm (Schoenfeld).

    Median split (p=0.5) and standardised continuous covariate (SD=1).
    A covariate with R^2=0.2 inflation is shown as the 'adjusted' column.
    """
    z = norm.ppf(0.975) + norm.ppf(0.80)
    rows = []
    for label, var in [("median split (binary, p1=p2=0.5)", 0.25),
                       ("continuous, per 1 SD", 1.0),
                       ("continuous, per 1 SD, adjusted R2=0.2", 1.0 * (1 - 0.2))]:
        if n_events > 0:
            ln_hr = np.sqrt(z ** 2 / (var * n_events))
            rows.append(dict(scenario=label, n_cases=n_cases, n_events=n_events,
                             min_detectable_HR=round(float(np.exp(ln_hr)), 2),
                             events_per_variable=round(n_events / 2, 1)))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- GEO
def esummary(acc):
    u = (EUTILS + "esearch.fcgi?db=gds&term=" + urllib.parse.quote(f"{acc}[Accession]")
         + "&retmode=json")
    with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=90) as r:
        ids = json.loads(r.read())["esearchresult"]["idlist"]
    if not ids:
        return None
    u2 = EUTILS + "esummary.fcgi?db=gds&id=" + ",".join(ids) + "&retmode=json"
    with urllib.request.urlopen(urllib.request.Request(u2, headers=UA), timeout=90) as r:
        res = json.loads(r.read())["result"]
    return res[res["uids"][0]]


def gsm_soft(gsm):
    u = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"
         "&targ=self&form=text&view=brief")
    return urllib.request.urlopen(urllib.request.Request(u, headers=UA),
                                 timeout=90).read().decode(errors="replace")


# ------------------------------------------------------------------ title parsing
# Each parser maps GEO sample titles to a list of (donor id, library class) pairs.
# They are written per accession because series number their patients differently;
# a single generic regex silently merges "C15" with "N15" and is not usable here.

def _parse_gse181919(titles):
    """C##/N##/LN## share one patient number (matched normal / metastasis); LP## are separate."""
    out = []
    for t in titles:
        m = re.search(r"\b(C|LN|N|LP)(\d+)\b", t)
        if not m:
            continue
        pre, num = m.groups()
        out.append((("LP" + num) if pre == "LP" else num,
                    {"C": "primary_cancer", "LN": "metastatic_ln",
                     "N": "normal", "LP": "leukoplakia"}[pre]))
    return out


def _parse_gse198315(titles):
    """P##-TC/TP/NT/mLN: tumour core / periphery / normal / metastatic LN."""
    out = []
    for t in titles:
        m = re.match(r"(P\d+)-(\w+)", t)
        if m:
            out.append((m.group(1), {"TC": "tumour_core", "TP": "tumour_periphery",
                                     "NT": "normal", "mLN": "met_ln"}.get(m.group(2), m.group(2))))
    return out


def _parse_gse234933(titles):
    out = []
    for t in titles:
        m = re.search(r"HN(\d+)", t)
        if m:
            out.append(("HN" + m.group(1),
                        "primary" if "Primary" in t else
                        "recurrence" if "Recurrence" in t else "metastasis"))
    return out


def _parse_gse215403(titles):
    out = []
    for t in titles:
        m = re.search(r"\(P(\d+)\)", t)
        if m:
            out.append(("P" + m.group(1), "oscc_buccal"))
    return out


def _parse_gse310797(titles):
    out = []
    for t in titles:
        m = re.match(r"OSCC_(\d+)", t)
        if m:
            out.append(("OSCC_" + m.group(1), "scrna"))
        elif t.startswith("SP_"):
            out.append((t, "spatial"))
    return out


def _parse_gse195832(titles):
    """Mixed series: single-cell, bulk RNA, sorted CAF bulk and ex vivo CAF.

    Only the single-cell libraries can serve Arm A, so donors are reported per
    library class -- pooling them would inflate the donor count.
    """
    out = []
    for t in titles:
        if t.startswith("ScRNA"):
            m = re.search(r"Pt(\d+)", t)
            if m:
                out.append(("Pt" + m.group(1), "single_cell"))
        elif t.startswith("Bulk RNA"):
            m = re.search(r"Pt(\d+)", t)
            if m:
                out.append(("Pt" + m.group(1), "bulk_rna"))
        elif t.startswith("Sorted"):
            out.append((t, "sorted_caf_bulk"))
        elif t.startswith("Ex Vivo"):
            out.append((t, "ex_vivo_caf"))
    return out


def _parse_gse301741(titles):
    """OCSCC##/LSCC## = patient, followed by Pre/Post and sort fraction."""
    out = []
    for t in titles:
        m = re.match(r"(OCSCC|LSCC)(\d+)", t)
        if m:
            out.append((m.group(1) + m.group(2),
                        f"{m.group(1)}_" + ("pre" if ", Pre," in t else "post")))
    return out


def _parse_gse303690(titles):
    out = []
    for t in titles:
        m = re.search(r"HT(\d+[A-B]?)", t)
        if m:
            out.append(("HT" + m.group(1),
                        "caf" if t.startswith("CAF") else "normal_fibroblast"))
    return out


def _parse_gse159067(titles):
    out = []
    for t in titles:
        m = re.search(r"patient number (\d+)", t)
        if m:
            out.append(("pt" + m.group(1), "patient"))
    return out


def _parse_gse270047(titles):
    """Specimen code = first comma field; the rest of the title is a technical replicate."""
    return [(t.split(",")[0].strip(), "library") for t in titles]


def _parse_gse306800(titles):
    out = []
    for t in titles:
        m = re.match(r"(TD|TJU)(\d+)", t)
        if m:
            out.append((m.group(1) + m.group(2),
                        "pre" if "pre-treatment" in t else "post"))
    return out


PARSERS = {"GSE181919": _parse_gse181919, "GSE198315": _parse_gse198315,
           "GSE234933": _parse_gse234933, "GSE215403": _parse_gse215403,
           "GSE310797": _parse_gse310797, "GSE195832": _parse_gse195832,
           "GSE301741": _parse_gse301741, "GSE303690": _parse_gse303690,
           "GSE159067": _parse_gse159067, "GSE270047": _parse_gse270047,
           "GSE306800": _parse_gse306800}

# Curated annotation per resource. `evidence` records where chemistry / PSI
# compatibility was read; library and donor counts are re-derived on every run.
RESOURCES = [
    dict(accession="GSE181919", arm="A", system="10x", clinical="OSCC + leukoplakia + normal",
         chemistry="10x Chromium 3' v2", psi="pseudobulk",
         evidence="PMC full text (choi 2023): 10x 3' v2; 54 239 cells, 12 336 fibroblasts",
         role="Arm A primary: CAF subtype x event-level PSI"),
    dict(accession="GSE198315", arm="A", system="10x", clinical="HPV-negative OSCC, multiregional",
         chemistry="10x Chromium 3'", psi="pseudobulk",
         evidence="GEO summary: 268 131 cells, 10 patients, core/periphery/adjacent/LN",
         role="Arm A: tumour core vs periphery CAF states"),
    dict(accession="GSE234933", arm="A", system="10x", clinical="HNSCC cohort (subsite not in titles)",
         chemistry="10x Chromium 3'", psi="pseudobulk",
         evidence="GEO summary: unsorted, no cell-type enrichment; 52 libraries, 1 per patient",
         role="Arm A donor pool (largest unsorted HNSCC series)"),
    dict(accession="GSE215403", arm="A", system="10x", clinical="OSCC gingivo-buccal",
         chemistry="10x Chromium 3'", psi="pseudobulk",
         evidence="GEO summary (pmid 37792582)",
         role="Arm A replication (OSCC-GB subsite)"),
    dict(accession="GSE195832", arm="A", system="10x", clinical="HNSCC, pre/post ICB + sorted CAF bulk",
         chemistry="10x Chromium 3' + bulk RNA-seq", psi="pseudobulk",
         evidence="GEO titles: mixed series (ScRNA PtN Pre/Post, Bulk RNA, Sorted CAF, Ex vivo)",
         role="Arm A/C bridge: CAF subtypes vs ICB"),
    dict(accession="GSE301741", arm="A", system="10x", clinical="HNSCC/Oral cavity, pre/post pembrolizumab",
         chemistry="10x 3' v2 + 5' v2, CD3/CD45 sorted", psi="pseudobulk (limited stromal)",
         evidence="GEO titles: OCSCC/LSCC + sort fraction; sorted libraries deplete fibroblasts",
         role="Arm A/C: treatment context, low CAF yield"),
    dict(accession="GSE310797", arm="A", system="10x+Visium", clinical="OSCC, different stages",
         chemistry="10x scRNA + Visium CytAssist FFPE", psi="pseudobulk",
         evidence="GEO summary; PMC full text earlier verified for chemistry",
         role="Arm A spatial anchor"),
    dict(accession="GSE303690", arm="B", system="bulk", clinical="Matched CAF vs normal fibroblast",
         chemistry="bulk RNA-seq, NextSeq 2000, 75 bp SE", psi="isoform (Salmon/StringTie)",
         evidence="GEO titles: 6 patients x CAF/NF x 3 biol reps = 36 libraries",
         role="Arm B only source of paired isoform switching"),
    dict(accession="TCGA-HNSC", arm="C", system="bulk", clinical="HNSCC, all subsites",
         chemistry="STAR - Counts (precomputed splice junctions)", psi="event + isoform",
         evidence="GDC files endpoint: 566 Splice Junction Quantification files (~1.4 GB total)",
         role="Arm C primary: PSI vs CAF abundance and survival"),
    dict(accession="GSE159067", arm="C", system="targeted", clinical="HNSCC anti-PD-1/PD-L1, RECIST + OS/PFS",
         chemistry="HTG EdgeSeq Oncology Biomarker Panel (probe capture, 50 bp SE)",
         psi="none (no junction reads)",
         evidence="GSM4818541 soft: HTG EdgeSeq + Oncology Biomarker Panel, 50-cycle SE",
         role="Arm C: gene-level only - CAF score vs ICB response"),
    dict(accession="GSE270047", arm="C", system="targeted", clinical="OncoPrism-HNSCC analytical validation",
         chemistry="targeted RNA-seq, FFPE, analytical replicates",
         psi="none (analytical reps, not patients)",
         evidence="GEO titles: 434 libraries from 28 specimen codes, rep/lot/operator/day",
         role="NOT a cohort: analytical performance series"),
    dict(accession="GSE306800", arm="C", system="bulk", clinical="HPV+ HNSCC, durvalumab+metformin, pre/post",
         chemistry="bulk RNA-seq, NEB Ultra II, NovaSeq 6000", psi="event (small n)",
         evidence="GSM soft: RNA-Seq, NEB Ultra II; GEO titles pre/post + pathologic response",
         role="Arm C exploratory: paired pre/post PSI (n=18 libraries)"),
]


def geo_ledger():
    rows = []
    for r in RESOURCES:
        if r["accession"] == "TCGA-HNSC":
            rows.append({**{k: r[k] for k in ("accession", "arm", "system", "clinical",
                                              "chemistry", "psi", "evidence", "role")},
                         "n_libraries": np.nan, "n_donors": np.nan,
                         "donors_by_class": "", "libraries_by_class": ""})
            continue
        d = esummary(r["accession"])
        if d is None:
            print("!! not found:", r["accession"])
            continue
        titles = [s.get("title", "") for s in d.get("samples", [])]
        pairs = PARSERS[r["accession"]](titles)
        if len(pairs) != int(d.get("n_samples", 0)):
            print(f"   WARN {r['accession']}: parsed {len(pairs)} of "
                  f"{d.get('n_samples')} libraries")
        lib_cls = collections.Counter(c for _, c in pairs)
        donor_cls = collections.defaultdict(set)
        for donor, cls in pairs:
            donor_cls[cls].add(donor)
        rows.append({**{k: r[k] for k in ("accession", "arm", "system", "clinical",
                                          "chemistry", "psi", "evidence", "role")},
                     "n_libraries": d.get("n_samples"),
                     "n_donors": len({dn for dn, _ in pairs}),
                     "donors_by_class": ";".join(f"{c}:{len(v)}"
                                                 for c, v in sorted(donor_cls.items())),
                     "libraries_by_class": ";".join(f"{c}:{n}"
                                                    for c, n in sorted(lib_cls.items()))})
        print(f"{r['accession']:11s} libs={d.get('n_samples'):>4} "
              f"donors={len({dn for dn, _ in pairs}):>4} "
              f"donors_by_class={{{', '.join(f'{c}:{len(v)}' for c, v in sorted(donor_cls.items()))}}}")
    return pd.DataFrame(rows)


def main():
    print("== TCGA-HNSC splice-junction files")
    f = fetch_tcga_hnsc()
    print(f"files={len(f)} primary_tumour={int((f.sample_type == 'Primary Tumor').sum())} "
          f"total_GB={f.file_mb.sum() / 1000:.2f}")
    frames = [("oral_strict", ORAL_STRICT),
              ("oral_ambig_lip_oral_pharynx", ORAL_AMBIG),
              ("oral_wide", ORAL_STRICT + ORAL_AMBIG),
              ("oropharynx_HPV_related", OROPHARYNX),
              ("larynx", LARYNX)]
    led = subsite_ledger(f, frames)
    led.to_csv("results/tables/tcga_hnsc_site_ledger.csv", index=False)
    print(led.to_string(index=False))

    sites = site_breakdown(f)
    sites.to_csv("results/tables/tcga_hnsc_site_breakdown.csv", index=False)

    oral = led[led.group == "oral_strict"].iloc[0]
    pow_tab = pd.concat([
        cox_power_table(int(oral.n_cases), int(oral.n_events)),
        cox_power_table(int(led[led.group == "oral_wide"].iloc[0].n_cases),
                        int(led[led.group == "oral_wide"].iloc[0].n_events))],
        ignore_index=True)
    pow_tab.to_csv("results/tables/tcga_hnsc_oral_power.csv", index=False)
    print(pow_tab.to_string(index=False))

    print("\n== GEO resources")
    geo = geo_ledger()
    geo.to_csv("results/tables/arm_sample_ledger.csv", index=False)

    hist = histology_check()
    print("\n== TCGA-HNSC histology")
    print(hist.to_string(index=False))

    json.dump(dict(tcga_files=len(f),
                   tcga_primary_tumours=int((f.sample_type == "Primary Tumor").sum()),
                   tcga_total_GB=round(float(f.file_mb.sum()) / 1000, 2),
                   subsite_ledger=led.to_dict("records"),
                   site_breakdown=sites.to_dict("records"),
                   histology=hist.to_dict("records"),
                   power=pow_tab.to_dict("records"),
                   resources=geo.to_dict("records")),
              open("data/processed/arm_sample_ledger.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    print("\nwritten: results/tables/{tcga_hnsc_site_ledger,tcga_hnsc_oral_power,"
          "arm_sample_ledger}.csv, data/processed/arm_sample_ledger.json")


if __name__ == "__main__":
    main()
