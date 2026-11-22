# Vaxforge

An agent-assisted, in silico **reverse-vaccinology** pipeline with a web interface — from a pathogen sequence to a ranked multi-epitope mRNA vaccine construct and a fully cited report.

[![type](https://img.shields.io/badge/type-reverse%20vaccinology-0d6b8f)](https://github.com/aliarslan47/Vaxforge)
[![interface](https://img.shields.io/badge/interface-Streamlit%20%2B%20Next.js-2f8f5b)](https://github.com/aliarslan47/Vaxforge)
[![tools](https://img.shields.io/badge/tools-real%20%C2%B7%20cited-c07211)](https://github.com/aliarslan47/Vaxforge)

[Türkçe](README.tr.md) · **English**

## What is it?

Vaxforge is the vaccine-design member of the Forge family — a reverse-vaccinology pipeline that turns a pathogen file into ranked vaccine candidates. A deterministic scientific core does the biology; an LLM wrapper only plans, interprets and reports. It ships with two interfaces: a single-file Streamlit app for local analysis, and a FastAPI + Next.js stack for a multi-run web deployment.

## What it does

Upload a pathogen file (FASTA/FASTQ) → the system auto-detects the input → mines virulence factors and vaccine targets → runs the antigen funnel → predicts B- and T-cell epitopes across selectable hosts → filters for survival (toxicity, allergenicity) → ranks candidates with a configurable candidacy score → assembles a multi-epitope mRNA construct → produces a publication-style report with full tool citations.

- **Real, published tools**, each cited in the report (PDF/HTML/JSON) and UI.
- **Host-selectable, multi-organism MHC**: four hosts enabled in `config/hosts.yaml` (human, mouse, bovine, pig), every allele validated against the local NetMHCpan allele list. Pig MHC-II is deliberately disabled because NetMHCIIpan carries no SLA-DR models; a proxy would fabricate epitopes.
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

External tools live under `tools/` (git-ignored). Free tools auto-install (DIAMOND, VFDB, Swiss-Prot, ToxinPred2, IApred); licensed DTU HealthTech tools (NetMHCpan, NetMHCIIpan, SignalP-5.0, TMHMM-2.0, BepiPred-1.0) need a manual academic download — see `tools/README.md`. If a tool is missing, the pipeline still runs with a labeled fallback.

**After extracting any DTU HealthTech package, run this once:**

```bash
python3 scripts/configure_dtu_wrappers.py          # fix paths
python3 scripts/configure_dtu_wrappers.py --check  # verify only (exit 1 if not)
python3 scripts/setup_toxinpred.py                 # ToxinPred2: envfile, pandas API, RF model
python3 scripts/convert_ifn_models.py              # IFNepitope2 sklearn pickle bridge
```

**PSORTb (bacterial localization) runs in a rootless podman container.** If the
pull fails with `potentially insufficient UIDs or GIDs available in user namespace`,
the image needs a UID outside your `/etc/subuid` range. You do **not** need `sudo`
or `usermod`: pull it once with the overlay driver option that ignores out-of-range
chown failures, after which a plain `podman run` works.

```bash
podman --storage-opt ignore_chown_errors=true pull docker.io/brinkmanlab/psortb_commandline:1.0.2
podman run --rm docker.io/brinkmanlab/psortb_commandline:1.0.2 /usr/local/psortb/bin/psort --version
```

DTU packages ship with vendor **placeholder** paths that must be edited by hand
(`setenv NMHOME /tools/src/netMHCpan-4.2`, `setenv BP /usr/cbs/packages/...`).
Left unedited the tool never computes, and because `tools/` is git-ignored these
fixes are not kept in version control — re-extracting a package reintroduces the
placeholder. The pipeline now detects an unconfigured wrapper and reports the
tool as absent instead of silently falling back, and every run records a
**Tool status** section listing which external tools were really available.

## Usage

Local Streamlit app:

```bash
streamlit run app.py
```

Or the web stack, API on port 8000 and UI on port 3000:

```bash
cd web/backend && pip install -r requirements.txt && uvicorn main:app --port 8000
cd web/frontend && npm install && npm run dev
```

Drag in a FASTA/FASTQ file (or pick a sample), choose the pathogen profile and host(s), review thresholds, and run. Outputs: publication-style **PDF**, **HTML** dashboard, ranked **CSV**, peptide **FASTA**, **GenBank** mRNA construct, and full-run **JSON**.

## Benchmark

Head-to-head against two published reverse-vaccinology tools, every tool on default
parameters and identical input, no threshold tuning. Ground truth is licensed or
literature-confirmed protective antigens. Details, accessions and statistics live in
`benchmark/`.

Ten bacteria, full proteomes, 48 protective antigens recovered:

| Tool | Recall | Share |
|---|---|---|
| Vaxign-ML | 47/48 | 98% |
| VaxForge | 37/48 | 77% |
| NERVE 2.0 | 14/48 | 29% |

On the MenB exemplar, where the ground truth is the four licensed 4CMenB antigens,
VaxForge recovers 4/4 and NERVE 2.0 recovers 1/4.

The table reports the October 2026 run. It predates the Gram-positive localization
fix described below, so the VaxForge figure will change when the benchmark is re-run.
The competitor figures are unaffected.

**Second caveat, added November 2026.** That run also predates a set of tool-level
fixes found later: TMHMM, IApred and ToxinPred2 each reported themselves as ready
while never actually executing (a broken interpreter path, two missing Python
dependencies, and a stale model pickle). In this environment that meant transmembrane
topology and antigenicity were never measured and toxicity never filtered anything.
Whether the October run was affected **cannot be determined from the files**, because
the benchmark summaries of the time did not record which tools actually ran. They do
now, so this question is answerable for every future run. The VaxForge figure should
therefore be treated as provisional until the benchmark is re-run; the competitors run
in their own containers and are unaffected either way.

All eleven VaxForge misses are mechanistically attributed in
`benchmark/MISS_ATTRIBUTION.md`. Six are cytoplasmic moonlighting antigens removed by
the host-homology and localization gates, which is a deliberate trade-off that
suppresses cytoplasmic false positives. No scaling bug was found. AUC is reported only
as a secondary, caveated metric, because the positive counts are too small for the
differences to be significant.

### Known limitation, and the fix applied after that run

The bacterial allowed-localization list was written with Gram-negative anatomy in
mind. Gram-positive bacteria have no outer membrane, so the cytoplasmic membrane is
the outermost boundary and surface-anchored proteins are labelled `membrane` by
PSORTb. The single shared list therefore discarded legitimate Gram-positive surface
antigens. The allowed list is now selected by Gram type, so a Gram-positive run
accepts `extracellular`, `cell_wall` and `membrane`, while a Gram-negative run keeps
rejecting the shielded inner membrane. See the note in `config/thresholds.yaml`.

## Modules

The scientific core runs with real, cited tools. Construct-level structure prediction is live through the ESMAtlas API; only AlphaFold-grade peptide–MHC modelling and molecular dynamics remain deferred to a GPU machine.

| Step | Tool / method | Status |
|---|---|---|
| Input auto-detection | FASTA/FASTQ · nt/protein · genome/CDS/reads | built-in |
| Discovery (virulence factors) | DIAMOND + VFDB | ✅ real |
| Antigen funnel — localization | PSORTb 3.0 (bacteria) · Kyte-Doolittle + signal heuristic otherwise | ✅ real / labeled heuristic |
| Antigen funnel — transmembrane | TMHMM-2.0 | ✅ real |
| Antigen funnel — signal peptide | SignalP-5.0 | ✅ real |
| Antigen funnel — antigenicity | IApred | ✅ real |
| Antigen funnel — host homology (safety) | DIAMOND vs human Swiss-Prot | ✅ real |
| Epitopes — B-cell | BepiPred-1.0 | ✅ real |
| Epitopes — MHC-I / MHC-II | NetMHCpan / NetMHCIIpan (local + IEDB) | ✅ real |
| Survival — toxicity / allergenicity | ToxinPred2 · FAO/WHO 6-mer + UniProt | ✅ real |
| Candidacy scoring | weighted, configurable | ✅ real |
| mRNA construct | linkers + adjuvant + human codon-opt + GC/CAI | ✅ real |
| Antigen funnel — localization (bacteria) | PSORTb 3.0 (podman) | ✅ real |
| Antigen funnel — antigenicity cross-check | VaxiJen ACC method (independent 2nd evidence) | ✅ real |
| Safety — host mimicry | exact k-mer + BLOSUM62 conservative substitution | ✅ real |
| Safety — host DB auto-provisioning | UniProt Swiss-Prot + DIAMOND makedb | ✅ real |
| Conservation — multi-strain | NCBI Datasets auto-fetch + per-residue conservation | ✅ real |
| Epitopes — IFN-γ induction (HTL) | IFNepitope2 | ✅ real |
| Epitopes — population coverage | IEDB population coverage method | ✅ real |
| Glycosylation — glycan-masked epitopes | N-X-S/T sequon + O-site heuristic (soft flag) | ✅ real |
| Construct structure | ESMFold via ESMAtlas API (≤400 aa) → pLDDT | ✅ real |
| Construct surface | Shrake-Rupley RSA + conformational B patches | ✅ real |
| Manufacturability | E. coli codon optimization, CAI, GC%, restriction scan, pET vector | ✅ real |
| Immunization verdict | ImmForge engine, 9 species-calibrated | ✅ real |
| HLA–disease annotation | curated allele table, PMID/DOI per row (interpretation only) | ✅ real |
| Peptide–MHC structure + MD | AlphaFold + docking/MD | ⏸️ deferred (GPU) |

Full citations (`vaxforge/citations.py`) and project layout live in the repo.

---

Forge family: [RNAForge](https://github.com/aliarslan47/RNAForge) (bulk RNA-seq) · [BacForge](https://github.com/aliarslan47/BacForge) (bacteria) · [VirusForge](https://github.com/aliarslan47/VirusForge) (virus/phage) · [MicrobiomeForge](https://github.com/aliarslan47/MicrobiomeForge) (microbiome) · **Vaxforge** (reverse vaccinology) · [ImmForge](https://github.com/aliarslan47/ImmForge) (immune simulation) · [PipelineForge](https://github.com/aliarslan47/PipelineForge) (DAG generator).
