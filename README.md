# VaxForge 🧬

**An agent-assisted, in silico reverse vaccinology pipeline with a web interface.**

🌐 English | [🇹🇷 Türkçe](README.tr.md)

Upload a pathogen file (FASTA/FASTQ) → the system **auto-detects** the input type →
**plans** what to do → mines virulence factors and vaccine targets → predicts B- and
T-cell epitopes across selectable hosts → ranks the strongest peptide candidates with a
**candidacy score** → assembles a multi-epitope **mRNA vaccine construct** → produces a
publication-style report with full tool citations.

---

## Pipeline

| Step | Tool / method | Status |
|------|---------------|--------|
| Input auto-detection | FASTA/FASTQ, nt/protein, genome/CDS/reads | built-in |
| **Discovery** (virulence factors) | DIAMOND + VFDB | ✅ real |
| **Antigen funnel** — localization | DeepLoc-2.1 | ✅ real |
| **Antigen funnel** — transmembrane | TMHMM-2.0 | ✅ real |
| **Antigen funnel** — signal peptide | SignalP-5.0 | ✅ real |
| **Antigen funnel** — antigenicity | IApred | ✅ real |
| **Antigen funnel** — host homology (safety) | DIAMOND vs human Swiss-Prot | ✅ real |
| **Epitopes** — B-cell | BepiPred-1.0 | ✅ real |
| **Epitopes** — MHC-I / MHC-II | NetMHCpan / NetMHCIIpan (local + IEDB) | ✅ real |
| Multi-host MHC panel | human, mouse, bovine, pig, chicken | ✅ real |
| **Survival** — toxicity | ToxinPred2 | ✅ real |
| **Survival** — allergenicity | FAO/WHO 6-mer + UniProt allergens | ✅ real |
| Candidacy scoring | weighted, configurable | ✅ real |
| mRNA construct | linkers + adjuvant + human codon-opt + GC/CAI/ProtParam | ✅ real |
| Peptide–MHC structure + MD | AlphaFold + docking/MD | ⏸️ deferred (needs GPU) |

The scientific core runs with **real, published tools** — each is cited in the report
(PDF, HTML, JSON) and in the UI.

## Key design decisions

- **Deterministic pipeline + LLM wrapper:** scientific steps are fixed/reproducible; the
  LLM only interprets and reports.
- **Host-selectable, multi-organism MHC:** every peptide is tested against the MHC alleles
  of one or more chosen hosts → a per-host presentation map ("which organism it works in").
- **No magic numbers:** all thresholds live in `config/thresholds.yaml` with organism
  presets; editable in the UI and written into every report for reproducibility.
- **Honest labeling:** where a real tool is unavailable, a fallback is used and clearly
  labeled as such.

## Installation

```bash
git clone git@github.com:aliarslan47/Vaxforge.git
cd Vaxforge
python3 -m venv --system-site-packages .venv
. .venv/bin/activate
pip install -r requirements.txt
```

### External tools (`tools/`, not included in the repo)

Open tools are auto-installable; some (DTU HealthTech) require an academic license and a
manual download. See `tools/README.md`. Summary:

- **Free / auto:** DIAMOND, VFDB, human Swiss-Prot, UniProt allergens, ToxinPred2 (`pip`),
  IApred (GitHub), DeepLoc-2.1 (`pip`).
- **Licensed / manual (DTU, Linux x86_64):** NetMHCpan, NetMHCIIpan, SignalP-5.0,
  TMHMM-2.0, BepiPred-1.0 — download and place under `tools/`.

If a tool is missing, the pipeline still runs using a clearly-labeled fallback.

## Usage

```bash
streamlit run app.py
```

Open the browser UI, drag in a FASTA/FASTQ file (or click a sample), pick the pathogen
profile and host(s), review/edit thresholds, and run. Outputs:

- **PDF** (publication-style, with References)
- **HTML** dashboard
- **CSV** ranked candidates, **FASTA** peptides, **GenBank** mRNA construct, **JSON** full run

## Project layout

```
config/            thresholds.yaml (per-organism presets), hosts.yaml (MHC alleles)
vaxforge/          pipeline modules (detect, discovery, funnel, epitope, survival,
                   scoring, mrna, report, report_pdf, citations, ...)
app.py             Streamlit interface
data/samples/      example FASTA/FASTQ
tools/             external tools + databases (gitignored)
```

## Citations

All tools used are cited in `vaxforge/citations.py` and reproduced in every report
(DIAMOND, VFDB, DeepLoc, TMHMM, SignalP, IApred, NetMHCpan/NetMHCIIpan, BepiPred,
ToxinPred2, FAO/WHO, and fallback methods) with DOIs.

## Status

Prototype. The full scientific core is real; GPU-dependent structural steps (AlphaFold
peptide–MHC, molecular dynamics) are deferred. Performance note: DeepLoc runs ESM on CPU,
so a full run takes a few minutes for a small proteome.
