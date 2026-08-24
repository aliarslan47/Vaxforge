# Vaxforge

An agent-assisted, in silico **reverse-vaccinology** pipeline with a web interface — from a pathogen sequence to a ranked multi-epitope mRNA vaccine construct and a fully cited report.

[![type](https://img.shields.io/badge/type-reverse%20vaccinology-0d6b8f)](https://github.com/aliarslan47/Vaxforge)
[![interface](https://img.shields.io/badge/interface-Streamlit%20web-2f8f5b)](https://github.com/aliarslan47/Vaxforge)
[![tools](https://img.shields.io/badge/tools-real%20%C2%B7%20cited-c07211)](https://github.com/aliarslan47/Vaxforge)

[Türkçe](README.tr.md) · **English**

## What is it?

Vaxforge is the vaccine-design member of the Forge family — a reverse-vaccinology pipeline that turns a pathogen file into ranked vaccine candidates. A deterministic scientific core does the biology; an LLM wrapper only plans, interprets and reports. It ships with a Streamlit web interface.

## What it does

Upload a pathogen file (FASTA/FASTQ) → the system auto-detects the input → mines virulence factors and vaccine targets → runs the antigen funnel → predicts B- and T-cell epitopes across selectable hosts → filters for survival (toxicity, allergenicity) → ranks candidates with a configurable candidacy score → assembles a multi-epitope mRNA construct → produces a publication-style report with full tool citations.

- **Real, published tools**, each cited in the report (PDF/HTML/JSON) and UI.
- **Host-selectable, multi-organism MHC** (human, mouse, bovine, pig, chicken) → per-host presentation map.
- **No magic numbers**: thresholds live in `config/thresholds.yaml` (organism presets), editable in the UI and written into every report.
- **Honest labeling**: where a real tool is unavailable, a clearly-labeled fallback is used.

## Installation

```bash
git clone https://github.com/aliarslan47/Vaxforge.git
cd Vaxforge
python3 -m venv --system-site-packages .venv
. .venv/bin/activate
pip install -r requirements.txt
```

External tools live under `tools/` (git-ignored). Free tools auto-install (DIAMOND, VFDB, Swiss-Prot, ToxinPred2, IApred, DeepLoc-2.1); licensed DTU HealthTech tools (NetMHCpan, NetMHCIIpan, SignalP-5.0, TMHMM-2.0, BepiPred-1.0) need a manual academic download — see `tools/README.md`. If a tool is missing, the pipeline still runs with a labeled fallback.

## Usage

```bash
streamlit run app.py
```

Drag in a FASTA/FASTQ file (or pick a sample), choose the pathogen profile and host(s), review thresholds, and run. Outputs: publication-style **PDF**, **HTML** dashboard, ranked **CSV**, peptide **FASTA**, **GenBank** mRNA construct, and full-run **JSON**.

## Modules

The scientific core runs with real, cited tools; GPU-dependent structural steps are deferred.

| Step | Tool / method | Status |
|---|---|---|
| Input auto-detection | FASTA/FASTQ · nt/protein · genome/CDS/reads | built-in |
| Discovery (virulence factors) | DIAMOND + VFDB | ✅ real |
| Antigen funnel — localization | DeepLoc-2.1 | ✅ real |
| Antigen funnel — transmembrane | TMHMM-2.0 | ✅ real |
| Antigen funnel — signal peptide | SignalP-5.0 | ✅ real |
| Antigen funnel — antigenicity | IApred | ✅ real |
| Antigen funnel — host homology (safety) | DIAMOND vs human Swiss-Prot | ✅ real |
| Epitopes — B-cell | BepiPred-1.0 | ✅ real |
| Epitopes — MHC-I / MHC-II | NetMHCpan / NetMHCIIpan (local + IEDB) | ✅ real |
| Survival — toxicity / allergenicity | ToxinPred2 · FAO/WHO 6-mer + UniProt | ✅ real |
| Candidacy scoring | weighted, configurable | ✅ real |
| mRNA construct | linkers + adjuvant + human codon-opt + GC/CAI | ✅ real |
| Peptide–MHC structure + MD | AlphaFold + docking/MD | ⏸️ deferred (GPU) |

Full citations (`vaxforge/citations.py`) and project layout live in the repo.

---

Forge family: [RNAForge](https://github.com/aliarslan47/RNAForge) (bulk RNA-seq) · [BacForge](https://github.com/aliarslan47/BacForge) (bacteria) · [VirusForge](https://github.com/aliarslan47/VirusForge) (virus/phage) · [MicrobiomeForge](https://github.com/aliarslan47/MicrobiomeForge) (microbiome) · **Vaxforge** (reverse vaccinology) · [ImmForge](https://github.com/aliarslan47/ImmForge) (immune simulation) · [PipelineForge](https://github.com/aliarslan47/PipelineForge) (DAG generator).
