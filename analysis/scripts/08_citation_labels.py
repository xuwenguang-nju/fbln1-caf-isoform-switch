"""Step 8: fetch the exact first-author surname + year for every cited DOI, so that
inline citations in the report can be checked against PubMed records.

Inputs : results/tables/as_bioinformatics_exemplars.csv
Outputs: results/tables/as_citation_labels.csv   (doi, first_author, year, journal, title)
"""
import csv, json, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DOIS_TABLE = "results/tables/as_bioinformatics_exemplars.csv"
OUT = "results/tables/as_citation_labels.csv"


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


def main():
    dois = [r["doi"] for r in csv.DictReader(open(DOIS_TABLE, encoding="utf-8-sig"))]
    rows = []
    for doi in dois:
        res = json.loads(eutil("esearch.fcgi", {"db": "pubmed", "term": f'{doi}[doi]',
                                                "retmax": 3, "retmode": "json"}))
        ids = res["esearchresult"]["idlist"]
        if not ids:
            rows.append({"doi": doi, "first_author": "", "year": "", "journal": "", "title": ""})
            print("NO PM RECORD", doi)
            continue
        xml = eutil("efetch.fcgi", {"db": "pubmed", "id": ids[0], "retmode": "xml"})
        art = ET.fromstring(xml).find(".//PubmedArticle")
        surname = art.findtext(".//AuthorList/Author/LastName") or \
                  art.findtext(".//AuthorList/Author/CollectiveName") or ""
        initials = art.findtext(".//AuthorList/Author/Initials") or ""
        year = (art.findtext(".//JournalIssue/PubDate/Year")
                or (art.findtext(".//JournalIssue/PubDate/MedlineDate") or "")[:4])
        rows.append({"doi": doi, "first_author": surname, "initials": initials, "year": year,
                     "journal": art.findtext(".//Journal/Title") or "",
                     "title": " ".join((art.findtext(".//ArticleTitle") or "").split())})
        print(f"{surname} {year} | {doi}")
        time.sleep(0.35)
    with open(OUT, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["doi", "first_author", "initials", "year", "journal", "title"])
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
