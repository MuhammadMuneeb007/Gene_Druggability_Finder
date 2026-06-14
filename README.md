# Gene Druggability Data and Feature Engineering Pipeline

This repository contains the data preparation and feature engineering code for a
human gene druggability modelling project. It builds:

1. A gene-level druggability label/evidence table from known target databases.
2. A collection of independent biological, genetic, structural, expression, and
   network feature tables.
3. Gene-level HGNC-merged outputs that can be combined into modelling datasets.

The most important design rule is the separation between **label evidence** and
**predictive features**:

- `Dataset0-DownloadDatabases.py` and `Dataset1-GenerateData.py` use known
  druggability evidence to create the target labels.
- `Feature0-DownloadDatabase.py` and `Feature1...Feature22` create predictors
  intended to avoid directly leaking those labels into the model.

## Repository Status

The current repository contains source scripts and one supplementary workbook.
Downloaded databases and generated feature directories are not currently stored
beside the scripts and must be downloaded or supplied before running the full
pipeline.

All Python files in the repository were syntax-checked successfully with
`python -m py_compile`.

`Supplementary Material 2.xlsx` currently contains:

- **Feature Statistics:** 2,285 feature-statistic records.
- **Dataset Statistics:** 9,680 modelling-dataset records.
- **Model Performance:** an empty placeholder sheet.

## Pipeline Overview

```text
Known druggability databases
    |
    +-- Dataset0-DownloadDatabases.py
    |
    +-- Dataset1-GenerateData.py
            |
            +-- Step0_Output/03_HumanGene_DruggabilityLabels.csv

Independent biological databases
    |
    +-- Feature0-DownloadDatabase.py
    |
    +-- Feature1 ... Feature22
            |
            +-- feature-specific gene tables
            +-- feature-specific HGNC-merged tables
```

Most feature scripts use this common HGNC input:

```text
databases/HGNC/hgnc_complete_set.txt
```

Unless `--all-hgnc-genes` or an equivalent option is supplied, most scripts
focus on protein-coding genes.

## Installation

Python 3.10 is recommended.

```bash
conda create -n druggability python=3.10 -y
conda activate druggability
pip install pandas numpy requests tqdm pyarrow networkx h5py scikit-learn openpyxl
```

Important optional system tools:

- `fpocket`: required only for real pocket scoring in Feature 10.
- SLURM: used by `Feature10-FPocket.sh` for HPC array execution.
- `curl` or `wget`: optional download fallbacks used by some scripts.

`Dataset1-GenerateData.py` reads Open Targets parquet data, so pandas needs a
parquet engine such as `pyarrow`.

## Recommended Execution Order

Run commands from the repository root.

### 1. Download Label/Evidence Databases

```bash
python Dataset0-DownloadDatabases.py --all --opentargets-release 26.03
```

This downloads HGNC, ChEMBL, Open Targets, DGIdb, and TCRD/Pharos data where
automatic download is possible.

Downloaded valid files are skipped by default. Use `--overwrite` only when a
deliberate refresh is required.

### 2. Build the Ground-Truth Label Table

```bash
python Dataset1-GenerateData.py --save-raw-evidence
python Dataset1.1-GenerateDataCheckFeatures.py
```

Primary output:

```text
Step0_Output/03_HumanGene_DruggabilityLabels.csv
```

This file is the target/evidence table. Do not use its ChEMBL, Open Targets,
DGIdb, or Pharos evidence columns as model input features.

### 3. Download Core No-Leakage Feature Databases

Recommended first run:

```bash
python Feature0-DownloadDatabase.py --recommended
```

This downloads the most useful core inputs:

- HGNC
- UniProt
- Ensembl
- STRING
- GTEx

Large optional downloads can be added later:

```bash
python Feature0-DownloadDatabase.py --recommended --include-alphafold
python Feature0-DownloadDatabase.py --sources interpro --include-interpro-huge --include-pfam
python Feature0-DownloadDatabase.py --full
```

The downloader safely skips existing valid files unless `--overwrite` is used.
AlphaFold archives are dynamically discovered from the latest EBI directory,
and the current Pfam bulk filenames are used.

### 4. Build Feature Tables

Run feature scripts individually. A sensible starting sequence is:

```bash
python Feature1_DeMap.py
python Feature2_String.py
python Feature3_Pathway.py
python Feature4_AlphaFold.py
python Feature5_InterProPfam.py
python Feature6_UniProt.py
python Feature7_GTEx.py
python Feature8_Ensmbl.py
python Feature9_Genetics_Contraint.py
```

Then run the remaining feature scripts as their source data becomes available.
Most scripts support a fast-test option such as `--limit-genes`, `--max-rows`,
or `--max-lines`.

There is currently no single orchestration script that executes all feature
builders or merges every final feature table into one master matrix.

## Expected Directory Structure

The pipeline creates or expects a structure similar to:

```text
databases/
    HGNC/
    ChEMBL/
    OpenTargets/
    DGIdb/
    TCRD_Pharos/

feature_databases/
    AlphaFold/
    BioGRID/
    CORUM/
    CTD/
    Ensembl/
    Ensembl_Paralogues/
    GeneOntology/
    GTEx/
    GWASCatalog/
    HPA/
    InterPro_Pfam/
    MGI/
    PhosphoSitePlus/
    ProteinEmbeddings/
    STRING/
    UniProt/
    gnomAD/

Step0_Output/
feature1.dmapp.database/
feature2.string.database/
feature3_pathway/
...
feature22_paralogues/
```

Most feature output directories contain:

- `downloads/` or `raw/`: cached source data.
- `processed/`: long-form, gene-level, and HGNC-merged tables.
- `*_summary.txt`: a human-readable run summary.
- `*_run_metadata.json`: run settings and provenance.

## Label and Dataset Scripts

| Script | Purpose | Main outputs |
|---|---|---|
| `Dataset0-DownloadDatabases.py` | Safely downloads label/evidence databases: HGNC, ChEMBL, Open Targets, DGIdb, and TCRD/Pharos. | `databases/`, download metadata, and manual-download notes where needed. |
| `Dataset1-GenerateData.py` | Builds the human gene-level druggability ground-truth/evidence table. Optionally queries the Pharos GraphQL API. | `Step0_Output/03_HumanGene_DruggabilityLabels.csv`, label summary, raw evidence, and API caches. |
| `Dataset1.1-GenerateDataCheckFeatures.py` | Prints diagnostics, label distributions, missingness, top genes, and sanity checks for the generated label table. | Console diagnostic report. |

Useful label-builder options:

```bash
python Dataset1-GenerateData.py --skip-pharos-api
python Dataset1-GenerateData.py --max-pharos-genes 100
python Dataset1-GenerateData.py --all-locus-types
```

## Feature Database Downloader

`Feature0-DownloadDatabase.py` downloads non-label feature sources. It
deliberately excludes ChEMBL, Open Targets tractability/knownDrugs, Pharos/TCRD
labels, DGIdb, DrugBank, and Guide to Pharmacology drug-target evidence.

Modes:

```bash
python Feature0-DownloadDatabase.py --minimal
python Feature0-DownloadDatabase.py --recommended
python Feature0-DownloadDatabase.py --full
python Feature0-DownloadDatabase.py --sources hgnc uniprot ensembl string gtex
```

Important behavior:

- Existing valid files are skipped.
- Partial `.part` downloads are resumed.
- `--overwrite` forces replacement.
- AlphaFold human archives are discovered dynamically.
- AlphaFold, full InterPro, and Pfam downloads are optional because they are
  very large.

## Feature Script Catalog

| Feature | Script | Biological signal | Primary output directory |
|---|---|---|---|
| 1 | `Feature1_DeMap.py` | DepMap CRISPR dependency, expression, mutation, copy-number, and optional functional-genomics summaries. | `feature1.dmapp.database/` |
| 1 download-only helper | `Feature1_DepMap.py` | Downloads selected DepMap release files through the Figshare API without feature generation. | User-selected output directory. |
| 2 | `Feature2_String.py` | STRING mapping, interaction burden, evidence channels, and network topology. | `feature2.string.database/` |
| 3 | `Feature3_Pathway.py` | Reactome pathway membership, hierarchy, and pathway diversity. | `feature3_pathway/` |
| 4 | `Feature4_AlphaFold.py` | Experimental PDB and AlphaFold structure availability, confidence, and geometry. | `feature4_structure/` |
| 5 | `Feature5_InterProPfam.py` | InterPro/Pfam domain, family, repeat, and architecture features. | `feature5_interpro_pfam/` |
| 6 | `Feature6_UniProt.py` | UniProt annotation, localization, sequence, motif, and functional text features. | `feature6_uniprot/` |
| 7 | `Feature7_GTEx.py` | GTEx tissue expression, specificity, entropy, and breadth. | `feature7_gtex/` |
| 8 | `Feature8_Ensmbl.py` | Ensembl genomic coordinates, transcripts, exons, CDS/UTR, and FASTA length summaries. | `feature8_ensembl/` |
| 9 | `Feature9_Genetics_Contraint.py` | Protein disorder plus gnomAD genetic constraint/intolerance. | `feature9_disorder_constraint/` |
| 10 | `Feature10-FPocket.py` | Per-chunk structure geometry and optional fpocket pocket scores. | `feature10_pocket_geometry/processed/chunks/` |
| 10 merge | `Feature10-FPocket-merge.py` | Merges Feature 10 array-job chunks into final gene-level tables. | `feature10_pocket_geometry/` |
| 11 | `Feature11-ProteinFeatures.py` | Protein sequence composition and physicochemical features from Ensembl peptide FASTA. | `feature11_protein_sequence/` |
| 12 | `Feature12-ProteinSequence.py` | UniProt ProtT5 embeddings and PCA-reduced gene-level embedding features. | `feature12_protein_embeddings/` |
| 13 | `Feature13_GWAS.py` | GWAS Catalog association burden and pleiotropy. | `feature13_gwas_catalog/` |
| 14 | `Feature14_GO.py` | Gene Ontology annotation burden, categories, evidence, diversity, and optional depth. | `feature14_gene_ontology/` |
| 15 | `Feature15_BioGRID.py` | Curated BioGRID physical/genetic interaction-network features. | `feature15_biogrid/` |
| 16 | `Feature16_HPA.py` | Human Protein Atlas localization, RNA, tissue, immune, single-cell, and protein expression. | `feature16_hpa/` |
| 17 | `Feature17_CTD.py` | CTD chemical-gene interaction burden. | `feature17_ctd/` |
| 18 | `Feature18_MGI.py` | Mouse orthology, knockout phenotype, lethality, and phenotype burden. | `feature18_mgi/` |
| 19 | `Feature19_gnomAD_full.py` | Expanded gnomAD gene-constraint metrics. | `feature19_gnomad_full/` |
| 20 | `Feature20_CORUM.py` | CORUM protein-complex membership and complex-level summaries. | `feature20_corum/` |
| 21 | `Feature21_PhosphoSitePlus.py` | PhosphoSitePlus PTM-site and kinase-substrate features. | `feature21_phosphositeplus/` |
| 22 | `Feature22_Paralogues.py` | Ensembl BioMart human paralogue count and sequence-identity features. | `feature22_paralogues/` |

## Feature-Specific Notes

### Features 1-3: Functional Genomics and Networks

- `Feature1_DeMap.py` is the complete DepMap feature builder. It defaults to
  `DepMap Public 26Q1`.
- `Feature1_DepMap.py` only downloads DepMap files; it does not create features.
- `Feature2_String.py` can resume from cached API data with `--no-api`.
  Expensive graph calculations can be disabled with `--skip-expensive-graph`.
- `Feature3_Pathway.py` downloads Reactome mapping files once and reuses them.

### Features 4 and 10: Structures and Pockets

`Feature4_AlphaFold.py` downloads per-protein AlphaFold/PDB structures and builds
gene-level structural features.

Feature 10 depends on Feature 4 structures and is designed for chunked HPC use:

```bash
sbatch Feature10-FPocket.sh
python Feature10-FPocket-merge.py
```

For local testing:

```bash
python Feature10-FPocket.py 1 --chunk-size 100
python Feature10-FPocket.py 1 --chunk-size 100 --run-fpocket
```

The shell script currently requests a 200-task SLURM array, 50 GB memory per
task, and a 24-hour time limit.

### Features 5-9 and 11-12: Protein Biology

- Feature 5 expects InterPro/Pfam files downloaded by Feature 0.
- Feature 6 uses local UniProt data first and can fall back to the UniProt REST
  API.
- Feature 7 expects a GTEx gene median TPM expression file.
- Feature 8 expects Ensembl GTF, cDNA FASTA, and peptide FASTA files.
- Feature 9 combines DisProt, MobiDB, gnomAD, and optional Feature 4 AlphaFold
  confidence proxies.
- Feature 11 computes direct sequence composition and physicochemical features.
- Feature 12 downloads a large UniProt ProtT5 HDF5 embedding file and requires
  `h5py` plus scikit-learn.

### Features 13-22: Extended Biological Evidence

Several scripts download their source data only when `--download` is supplied:

```bash
python Feature13_GWAS.py --download
python Feature14_GO.py --download
python Feature15_BioGRID.py --download
python Feature17_CTD.py --download
python Feature18_MGI.py --download
python Feature21_PhosphoSitePlus.py --download
```

Other behavior:

- Feature 16 automatically discovers/downloads many HPA files and supports
  manual file paths.
- Feature 19 reuses the gnomAD constraint file used by Feature 9.
- Feature 20 is local-file-first and expects CORUM files in
  `feature_databases/CORUM/`.
- Feature 22 queries Ensembl BioMart in batches and caches the raw paralogue
  table. Use `--no-download` to require the cache.

## Naming Notes

Some filenames differ from the script names shown inside their docstrings.
Always run the filenames that actually exist in this repository:

- Run `Feature11-ProteinFeatures.py`; it implements Feature 11 protein sequence
  composition.
- Run `Feature12-ProteinSequence.py`; it implements Feature 12 ProtT5 protein
  embeddings.
- Run `Feature8_Ensmbl.py`; the filename contains `Ensmbl`, while the internal
  feature name is Ensembl.
- Run `Feature9_Genetics_Contraint.py`; the filename contains `Contraint`, while
  the internal feature name is disorder/constraint.
- Run `Feature4_AlphaFold.py`; its internal documentation calls it
  `Feature4_Structure.py`.
- Run `Feature0-DownloadDatabase.py`; its internal documentation calls it
  `Step0B_Download_NoLeakage_Feature_Databases.py`.

## Data Leakage Policy

The ground-truth pipeline intentionally uses known drug-target evidence. Those
columns must remain on the label side of the modelling boundary.

Do not use the following as predictive input features:

- ChEMBL activity or target evidence.
- DGIdb drug-gene interactions.
- Open Targets tractability or known-drug evidence.
- Pharos/TCRD `Tclin`, `Tchem`, `Tbio`, or `Tdark` labels.
- Approved-drug counts or direct known-target flags.

Most feature scripts explicitly exclude these sources. However, GWAS, CTD,
pathway, disease, and annotation-derived features can still encode indirect
knowledge about well-studied genes. Their use should be documented and evaluated
with sensitivity analyses.

## Reproducibility and Caching

- Downloaders skip existing valid files unless explicitly forced.
- Most API-based feature scripts cache downloaded/raw responses.
- Summary text and JSON metadata files record run settings.
- Use the same HGNC file across all features to preserve a consistent gene
  universe.
- Keep release names and download dates with final modelling datasets.
- Avoid deleting raw or download directories until all processed outputs have
  been validated.

## Supplementary Workbook

`Supplementary Material 2.xlsx` is a summary workbook with these sheets:

| Sheet | Contents |
|---|---|
| `Feature Statistics` | Feature source, full/short feature name, gene coverage, missingness, and numeric distribution statistics. |
| `Dataset Statistics` | Dataset ID/name, tier, feature group, target definition, feature count, gene count, class balance, prevalence, and missingness. |
| `Model Performance` | Currently empty and reserved for future results. |

The workbook summarizes previously generated datasets; it is not consumed by
the current Python scripts.

## Troubleshooting

### No sources selected

This is expected when a downloader is run without a mode or source list:

```text
No sources selected. Use --minimal, --recommended, --full, or --sources ...
```

Use:

```bash
python Feature0-DownloadDatabase.py --recommended
```

### Existing files are skipped

Messages such as the following are expected and protect against unnecessary
downloads:

```text
[SKIP] Existing file: ...
```

Use `--overwrite`, `--force-download`, or `--force` only when the relevant
script documents that option and a refresh is genuinely required.

### Missing input database

Check the script's default `feature_databases/<source>/` directory, provide the
manual file option where supported, or rerun with `--download`.

### API or FTP failures

External endpoints can change or be temporarily unavailable. Keep existing
download caches, retry later, or place the required source file manually in the
documented database directory.

### Large runs

Start with a script's test option before a full run:

```bash
python Feature4_AlphaFold.py --limit-genes 100
python Feature13_GWAS.py --max-rows 100000
python Feature18_MGI.py --max-pheno-rows 200000
python Feature22_Paralogues.py --limit-genes 1000
```

## Source File Inventory

The repository contains:

- 29 Python scripts.
- 1 SLURM shell script.
- 1 supplementary Excel workbook.
- Generated Python bytecode under `__pycache__/`.
- A temporary hidden Excel lock file may appear while the workbook is open.

Generated bytecode and temporary lock files are not pipeline inputs.
