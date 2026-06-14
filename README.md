# Gene Druggability Finder

Gene Druggability Finder (GDF) is a leakage-aware, multi-evidence framework for
predicting and benchmarking human gene druggability.

The pipeline:

1. Builds eight druggability target labels for 19,295 HGNC-approved
   protein-coding genes.
2. Generates 22 independent biological feature blocks.
3. Merges 2,285 feature columns into a gene-level feature matrix.
4. Creates 440 target-feature datasets.
5. Benchmarks Random Forest and XGBoost using stratified five-fold
   cross-validation.
6. Summarises model performance and feature importance.

> **Interpretation:** GDF predicts database-constructed druggability labels. It
> is a target-prioritisation and evidence-integration framework, not definitive
> proof that a gene will become a clinically successful drug target.

## Key Results

The reported internal benchmark covers 440 target-feature datasets:

| Result | Value |
|---|---:|
| Human protein-coding genes | 19,295 |
| Druggability targets | 8 |
| Feature blocks | 22 |
| Feature columns before model sanitisation | 2,285 |
| Feature-group configurations | 55 |
| Train-ready datasets | 440 |
| Best mean cross-validated AUROC | 0.9809 |
| Median AUROC | 0.8488 |
| Mean AUROC | 0.8344 +/- 0.1005 |
| Median AUPRC | 0.7476 |
| Median Matthews correlation coefficient | 0.4411 |
| XGBoost selected as best model | 350 datasets (79.5%) |
| Random Forest selected as best model | 90 datasets (20.5%) |

The highest AUROC was achieved by `Dataset413_GRP_NoStructure_T5`, which
predicted the Biologic/Modality Target with:

- AUROC: `0.9809`
- AUPRC: `0.9817`
- Matthews correlation coefficient: `0.8786`
- Best model: XGBoost

Integrated feature sets performed better than most individual feature blocks:

| Feature tier | Datasets | Best AUROC | Median AUROC | Median AUPRC | Median MCC |
|---|---:|---:|---:|---:|---:|
| Cumulative | 176 | 0.9809 | 0.9346 | 0.9045 | 0.6497 |
| Thematic | 88 | 0.9809 | 0.8289 | 0.7620 | 0.4367 |
| Individual feature blocks | 176 | 0.9421 | 0.7475 | 0.6377 | 0.3093 |

These are internal cross-validation results. External and prospective
validation are still required.

## Pipeline Design

GDF separates target-label construction from predictive feature generation:

```text
LABEL SIDE
HGNC + ChEMBL + Open Targets + DGIdb + Pharos/TCRD
    -> eight gene-level druggability targets

FEATURE SIDE
22 independent biological and molecular evidence blocks
    -> merged feature matrix

ANALYSIS SIDE
feature matrix + target labels
    -> 440 datasets
    -> Random Forest and XGBoost benchmark
    -> performance and feature-importance summaries
```

Direct target-defining evidence from ChEMBL, Open Targets, DGIdb and
Pharos/TCRD is reserved for the label side and excluded from model inputs.

## Repository Workflow

Run commands from the repository root.

### Step 1: Install Dependencies

Python 3.10 is recommended.

```bash
conda create -n druggability python=3.10 -y
conda activate druggability

pip install pandas numpy requests tqdm pyarrow openpyxl scipy \
  scikit-learn xgboost networkx h5py biopython matplotlib
```

Optional system tools:

- [fpocket](https://github.com/Discngine/fpocket) for Feature 10 pocket scores.
- SLURM for HPC array execution.
- `curl` or `wget` as download fallbacks used by some scripts.

### Step 2: Download Target-Label Databases

```bash
python Dataset0-DownloadDatabases.py --all --opentargets-release 26.03
```

The downloader skips valid files that are already present. Use `--overwrite`
only when a deliberate refresh is required.

### Step 3: Build Druggability Labels

```bash
python Dataset1-GenerateData.py --save-raw-evidence
python Dataset1.1-GenerateDataCheckFeatures.py
```

Primary output:

```text
Step0_Output/03_HumanGene_DruggabilityLabels.csv
```

### Step 4: Download Core Feature Databases

Recommended first download:

```bash
python Feature0-DownloadDatabase.py --recommended
```

This downloads HGNC, UniProt, Ensembl, STRING and GTEx inputs.

Large optional downloads:

```bash
python Feature0-DownloadDatabase.py --recommended --include-alphafold
python Feature0-DownloadDatabase.py --sources interpro --include-interpro-huge --include-pfam
python Feature0-DownloadDatabase.py --full
```

The downloader skips existing valid files, resumes `.part` downloads,
dynamically discovers the current AlphaFold human archive and uses current Pfam
bulk filenames.

### Step 5: Generate the 22 Feature Blocks

Run the feature scripts after their source databases are available:

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
python Feature11-ProteinFeatures.py
python Feature12-ProteinSequence.py
python Feature13_GWAS.py --download
python Feature14_GO.py --download
python Feature15_BioGRID.py --download
python Feature16_HPA.py
python Feature17_CTD.py --download
python Feature18_MGI.py --download
python Feature19_gnomAD_full.py
python Feature20_CORUM.py
python Feature21_PhosphoSitePlus.py --download
python Feature22_Paralogues.py
```

Feature 10 is designed as a chunked structure/pocket workflow:

```bash
sbatch Feature10-FPocket.sh
python Feature10-FPocket-merge.py
```

For a small local Feature 10 test:

```bash
python Feature10-FPocket.py 1 --chunk-size 100
python Feature10-FPocket.py 1 --chunk-size 100 --run-fpocket
```

### Step 6: Merge All Feature Blocks

```bash
python Analysis1-MergeData.py --verbose
```

Output:

```text
Dataset/Features_All.csv
```

The merger:

- Locates the expected gene-level output from each feature block.
- Normalises gene symbols.
- Prefixes columns with `Feature1_` through `Feature22_`.
- Excludes identifiers, raw text, paths, URLs and other non-model columns.
- Outer-merges all available feature blocks.

### Step 7: Generate Feature Statistics

```bash
python Analysis2-ListFeaturesStatistics.py \
  --features Dataset/Features_All.csv \
  --out Feature_Statistics.xlsx

python Analysis2.1-ListFeaturesStatistics.py
```

The first script calculates gene coverage, missingness and numeric summary
statistics for every feature. The second prints feature counts per block.

Do not set `--out` to the existing `Supplementary Material 2.xlsx` unless you
intend to replace it. `Analysis2-ListFeaturesStatistics.py` creates a new
workbook and would overwrite the compiled performance and feature-importance
sheets.

### Step 8: Create Train-Ready Datasets

```bash
python Analysis3-MakeDatasets.py
```

This crosses eight target definitions with 55 feature-group configurations:

- 22 individual groups: `F01` to `F22`
- 22 cumulative groups: `CUM01` to `CUM22`
- 11 thematic groups

Thematic groups:

| Group | Included feature blocks |
|---|---|
| `GRP_Structure` | 4, 10 |
| `GRP_Network` | 2, 15 |
| `GRP_Expression` | 7, 16 |
| `GRP_Constraint` | 9, 19 |
| `GRP_Pathway` | 3, 14 |
| `GRP_Functional` | 5, 6, 11, 12 |
| `GRP_Omics` | 1, 18 |
| `GRP_Annotation` | 8, 13, 20, 21, 22 |
| `GRP_Literature` | 17 |
| `GRP_NoStructure` | All except 4 and 10 |
| `GRP_Full` | All 22 blocks |

Each generated dataset contains:

```text
Datasets/Dataset001_F01_T1/
    X_train.csv
    y_train.csv
    meta.json
    feature_columns.txt
    missingness_summary.csv
```

Missing values remain in `X_train.csv` and are imputed inside training folds.

### Step 9: Train Models

Train one dataset locally:

```bash
python Analysis4-TrainModels.py 1 --no-permutation
```

Train all 440 datasets on SLURM:

```bash
sbatch Analysis4-TrainModels.sh
```

The training script:

- Removes identifiers and suspicious leakage-related columns.
- Converts selected predictors to numeric values.
- Replaces infinite and unsafe values with missing values.
- Performs median imputation inside each cross-validation fold.
- Evaluates Random Forest and XGBoost.
- Uses stratified five-fold cross-validation by default.
- Selects the best model using mean cross-validated AUROC.
- Writes fold predictions, metrics, confusion matrices and feature importance.

Main training outputs inside each dataset's `Training/` directory:

```text
00_training_metadata.json
00_feature_cleaning_report.csv
01_feature_columns_used.csv
02_target_info.csv
03_cv_all_predictions.csv
04_cv_fold_metrics.csv
05_model_comparison_summary.csv
06_best_model_summary.csv
07_feature_importance_all_folds.csv
08_feature_importance_summary.csv
09_top_features_best_model.csv
10_confusion_matrix_per_fold.csv
11_confusion_matrix_aggregate.csv
12_permutation_test_results.csv
13_leakage_audit.csv
```

### Step 10: Audit Dataset and Training Outputs

```bash
python Analysis4.0-CheckDatasets.py
```

`Analysis4.0-CheckDatasets.py` currently contains an HPC-specific absolute path:

```text
/data/ascher02/uqmmune1/DrugableGeneFinder/Final
```

Update `BASE_DIR` in that script before running it elsewhere.

## Eight Druggability Targets

The target table contains eight operational definitions of druggability:

| Target | Definition | Positive genes | Prevalence |
|---|---|---:|---:|
| T1 Clinical Target | Approved or clinically established target; the strictest primary definition. | 1,052 | 5.45% |
| T2 Clinical Investigation Target | Phase 1-3 clinical-investigation evidence. | 504 | 2.61% |
| T3 Small-Molecule Target | Small-molecule tractability evidence. | 7,231 | 37.48% |
| T4 Chemical Tractability Target | Broader chemical tractability evidence. | 7,058 | 36.58% |
| T5 Biologic/Modality Target | Biologic, antibody, targeted-degradation or other modality evidence. | 9,699 | 50.27% |
| T6 Drug-Gene Interaction Target | DGIdb-supported drug-gene interaction evidence. | 4,575 | 23.71% |
| T7 Potentially Druggable Family Target | Potentially druggable family/category annotation. | 9,518 | 49.33% |
| T8 Broad Druggability Target | Positive for any preceding druggability definition. | 14,791 | 76.66% |

T1 is the primary conservative target. T8 is a broad evidence-union label and
should not be interpreted as a strict clinical gold standard.

## Label-Construction Databases

| Database | Purpose | Access |
|---|---|---|
| HGNC | Defines the human gene universe and stable identifier mappings. | [HGNC complete set](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt) |
| ChEMBL | Human targets, mechanisms, approved/phase 4 drugs and potent activity evidence. | [ChEMBL releases](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/) |
| Open Targets | Target identity, tractability and clinical-indication evidence. | [Open Targets platform downloads](https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/) |
| DGIdb | Drug-gene interactions and druggable-category evidence. | [DGIdb downloads](https://dgidb.org/downloads) |
| Pharos/TCRD | Target development levels: TCLIN, TCHEM, TBIO and TDARK. | [Pharos GraphQL API](https://pharos-api.ncats.io/graphql) |

## Feature Blocks and Database Links

The current feature-statistics workbook contains 2,285 feature columns:

| Feature | Source | Columns | Access |
|---:|---|---:|---|
| 1 | DepMap functional genomics | 90 | [DepMap downloads](https://depmap.org/portal/download/all/) / [file-index API](https://depmap.org/portal/api/download/files) |
| 2 | STRING v12.0 PPI network | 48 | [STRING API](https://string-db.org/api) |
| 3 | Reactome pathways | 30 | [Reactome downloads](https://reactome.org/download/current/) |
| 4 | AlphaFold DB + RCSB PDB structures | 33 | [AlphaFold DB](https://alphafold.ebi.ac.uk/) / [RCSB PDB](https://www.rcsb.org/) |
| 5 | InterPro + Pfam domains | 72 | [InterPro](https://www.ebi.ac.uk/interpro/) / [Pfam FTP](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) |
| 6 | UniProt protein annotation | 91 | [UniProt REST API](https://rest.uniprot.org/) |
| 7 | GTEx tissue expression | 108 | [GTEx portal](https://gtexportal.org/) / [v8 median TPM](https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz) |
| 8 | Ensembl gene structure | 52 | [Ensembl human GTF](https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/) |
| 9 | gnomAD + MobiDB + DisProt constraint/disorder | 74 | [gnomAD](https://gnomad.broadinstitute.org/downloads) / [MobiDB](https://mobidb.org/api) / [DisProt](https://disprot.org/download) |
| 10 | fpocket binding-pocket features | 186 | [fpocket](https://github.com/Discngine/fpocket) |
| 11 | Protein physicochemical/sequence features | 463 | Derived locally from protein FASTA |
| 12 | Protein language-model embeddings + PCA | 77 | [UniProt embedding archive](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000005640_9606/per-protein.h5) |
| 13 | NHGRI-EBI GWAS Catalog | 43 | [GWAS association download](https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0?split=false) |
| 14 | Gene Ontology / QuickGO | 82 | [GOA human GAF](https://current.geneontology.org/annotations/goa_human.gaf.gz) / [go-basic.obo](https://purl.obolibrary.org/obo/go/go-basic.obo) |
| 15 | BioGRID curated interactions | 50 | [BioGRID latest organism release](https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip) |
| 16 | Human Protein Atlas | 252 | [HPA downloads](https://www.proteinatlas.org/download/tsv) |
| 17 | CTD chemical-gene interaction burden | 72 | [CTD chemical-gene interactions](https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz) |
| 18 | Mouse Genome Informatics phenotypes | 56 | [MGI reports](https://www.informatics.jax.org/downloads/reports/) |
| 19 | Full gnomAD constraint | 109 | [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) |
| 20 | CORUM protein complexes | 122 | [CORUM](https://mips.helmholtz-muenchen.de/corum/) |
| 21 | PhosphoSitePlus PTM sites | 153 | [PhosphoSitePlus downloads](https://www.phosphosite.org/downloads) |
| 22 | Ensembl BioMart/Compara paralogues | 22 | [Ensembl BioMart](https://www.ensembl.org/biomart/martservice) |
| **Total** | **All 22 blocks** | **2,285** | |

## Feature Script Catalog

| Script | Purpose | Main output |
|---|---|---|
| `Feature0-DownloadDatabase.py` | Downloads core non-label feature databases. | `feature_databases/` |
| `Feature1_DeMap.py` | DepMap functional-genomics feature builder. | `feature1.dmapp.database/` |
| `Feature1_DepMap.py` | DepMap download-only helper. | Downloaded DepMap files |
| `Feature2_String.py` | STRING network topology and evidence-channel features. | `feature2.string.database/` |
| `Feature3_Pathway.py` | Reactome pathway membership and hierarchy. | `feature3_pathway/` |
| `Feature4_AlphaFold.py` | AlphaFold and PDB structure features. | `feature4_structure/` |
| `Feature5_InterProPfam.py` | InterPro/Pfam domain architecture. | `feature5_interpro_pfam/` |
| `Feature6_UniProt.py` | UniProt annotation and sequence features. | `feature6_uniprot/` |
| `Feature7_GTEx.py` | GTEx tissue-expression features. | `feature7_gtex/` |
| `Feature8_Ensmbl.py` | Ensembl genomic annotation and gene structure. | `feature8_ensembl/` |
| `Feature9_Genetics_Contraint.py` | Disorder and genetic constraint/intolerance. | `feature9_disorder_constraint/` |
| `Feature10-FPocket.py` | Chunked geometry and optional fpocket scoring. | `feature10_pocket_geometry/processed/chunks/` |
| `Feature10-FPocket-merge.py` | Merges Feature 10 chunks. | `feature10_pocket_geometry/` |
| `Feature11-ProteinFeatures.py` | Protein sequence composition and physicochemical features. | `feature11_protein_sequence/` |
| `Feature12-ProteinSequence.py` | Protein embedding and PCA features. | `feature12_protein_embeddings/` |
| `Feature13_GWAS.py` | GWAS association burden and pleiotropy. | `feature13_gwas_catalog/` |
| `Feature14_GO.py` | Gene Ontology features. | `feature14_gene_ontology/` |
| `Feature15_BioGRID.py` | Curated interaction-network features. | `feature15_biogrid/` |
| `Feature16_HPA.py` | HPA expression and localisation features. | `feature16_hpa/` |
| `Feature17_CTD.py` | Aggregate CTD chemical-gene burden. | `feature17_ctd/` |
| `Feature18_MGI.py` | Mouse orthology and phenotype features. | `feature18_mgi/` |
| `Feature19_gnomAD_full.py` | Expanded gnomAD constraint metrics. | `feature19_gnomad_full/` |
| `Feature20_CORUM.py` | Protein-complex membership features. | `feature20_corum/` |
| `Feature21_PhosphoSitePlus.py` | PTM and kinase-substrate features. | `feature21_phosphositeplus/` |
| `Feature22_Paralogues.py` | Human paralogue features. | `feature22_paralogues/` |

## Analysis and Modelling Scripts

| Script | Purpose |
|---|---|
| `Analysis1-MergeData.py` | Merges all available feature blocks into `Dataset/Features_All.csv`. |
| `Analysis2-ListFeaturesStatistics.py` | Produces per-feature coverage, missingness and numeric statistics. |
| `Analysis2.1-ListFeaturesStatistics.py` | Prints feature counts per block. |
| `Analysis3-MakeDatasets.py` | Builds the 440 train-ready target-feature datasets. |
| `Analysis4-TrainModels.py` | Trains Random Forest and XGBoost for one dataset. |
| `Analysis4-TrainModels.sh` | Runs all 440 datasets as a SLURM array. |
| `Analysis4.0-CheckDatasets.py` | Audits dataset and training output completeness. |

## Target-Level Benchmark Results

Results below use the intended 440-dataset benchmark:

| Target | Best AUROC | Median AUROC | Median AUPRC | Median MCC | Best feature group |
|---|---:|---:|---:|---:|---|
| T1 Clinical Target | 0.9738 | 0.9142 | 0.5582 | 0.4521 | CUM20 |
| T2 Clinical Investigation Target | 0.9302 | 0.8491 | 0.1966 | 0.1739 | CUM22 |
| T3 Small-Molecule Target | 0.9394 | 0.8505 | 0.7849 | 0.5198 | CUM18 |
| T4 Chemical Tractability Target | 0.9397 | 0.8488 | 0.7726 | 0.5147 | CUM22 |
| T5 Biologic/Modality Target | 0.9809 | 0.8496 | 0.8666 | 0.5325 | GRP_NoStructure |
| T6 Drug-Gene Interaction Target | 0.8914 | 0.8197 | 0.6402 | 0.4277 | CUM19 |
| T7 Potentially Druggable Family Target | 0.9142 | 0.8249 | 0.8312 | 0.4819 | CUM20 |
| T8 Broad Druggability Target | 0.9557 | 0.8345 | 0.9418 | 0.4339 | GRP_Full |

The strict Clinical Target has low prevalence, so AUPRC and AUPRC enrichment
are especially important alongside AUROC. The broad T8 target has high
prevalence and should be interpreted as an evidence-union target.

## Feature Importance

The supplied results contain:

- 516,672 model-derived feature-importance rows.
- 21,168 top-feature rows.
- Feature importance outputs for all trained datasets.

Highly ranked individual predictors include:

- PDB structure availability and PDB count.
- Reactome Ensembl/UniProt mapping indicators.
- STRING weighted degree, k-core and approximate closeness centrality.
- UniProt transmembrane, membrane and signal-peptide annotations.
- Gene Ontology author-statement/evidence-code counts.
- DepMap expression summaries.

The supplied manuscript reports the highest-ranked evidence blocks as STRING
network features, DepMap functional genomics, Reactome pathways, InterPro/Pfam
domains, UniProt annotations and AlphaFold/PDB structure features. This supports
the central conclusion that predictive signal is distributed across multiple
complementary evidence sources.

These values describe model contribution, not causal biological importance.
Correlated predictors, annotation density and tree-model splitting behavior can
affect importance rankings.

## Data Leakage Controls

Direct target-defining evidence must never be used as model input:

- ChEMBL target/activity evidence.
- Open Targets tractability and known-drug evidence.
- DGIdb direct drug-gene evidence.
- Pharos/TCRD target-development labels.
- Approved-drug counts.
- Clinical-target labels.

The modelling pipeline also removes identifiers, text, accessions, paths, URLs,
raw fields and columns matching suspicious leakage terms.

CTD is retained only as aggregate chemical-gene perturbation burden features.
Chemical names and direct drug-label fields are excluded.

## Supplementary Workbook

`Supplementary Material 2.xlsx` currently contains:

| Sheet | Current contents |
|---|---|
| `Feature Statistics` | 2,285 feature rows plus header. |
| `Dataset Statistics` | Dataset composition, target, feature-count and missingness records. |
| `Model Performance` | Best-model performance summaries for the benchmark datasets. |
| `Feature Importance` | Global and within-block feature-importance rankings. |

### Workbook Consistency Note

The current workbook contains 441 rows in `Model Performance` because Dataset
ID 97 appears twice. The intended experimental design and reported benchmark
contain 440 unique dataset IDs. Results in this README were calculated after
deduplicating by Dataset ID.

Treat the workbook as a compiled results artifact and preserve a backup before
running scripts that write Excel output.

## Current Naming and Consistency Notes

Use the filenames that actually exist in this repository:

- `Feature4_AlphaFold.py` is internally described as the structure feature.
- `Feature8_Ensmbl.py` contains the filename spelling `Ensmbl`.
- `Feature9_Genetics_Contraint.py` contains the filename spelling `Contraint`.
- `Feature11-ProteinFeatures.py` implements protein physicochemical/sequence
  features.
- `Feature12-ProteinSequence.py` implements protein embeddings.

The workbook and supplied manuscript label Feature 12 as **ESM-2 embeddings**,
while the current Feature 12 script documentation and download URL refer to the
UniProt **ProtT5** per-protein embedding archive. Confirm the intended embedding
source before rerunning Feature 12 or reporting the method.

## Expected Directory Structure

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
...
feature22_paralogues/

Dataset/
    Features_All.csv

Datasets/
    dataset_catalogue.csv
    Dataset001_F01_T1/
    ...
    Dataset440_GRP_Structure_T8/
```

## Caching and Reproducibility

- Downloaders skip existing valid files unless explicitly forced.
- API-based feature scripts cache downloaded/raw responses where possible.
- Most feature scripts write summary text and JSON run metadata.
- Use the same HGNC release across all feature blocks.
- Archive database versions and download dates with final outputs.
- Keep raw downloads until processed outputs and mappings are validated.
- Impute missing values only inside model-training folds.
- Treat model performance as internal validation until external validation is
  completed.

## Troubleshooting

### Downloader reports no sources selected

Supply a mode:

```bash
python Feature0-DownloadDatabase.py --recommended
```

### Existing files are skipped

This is expected:

```text
[SKIP] Existing file: ...
```

Use `--overwrite`, `--force-download` or `--force` only when the relevant
script documents the option and a refresh is required.

### A feature source is unavailable

External APIs and FTP layouts can change. Keep caches, retry later or manually
place the expected source file in the documented `feature_databases/`
subdirectory.

### A full run is too large

Use each script's testing options first, for example:

```bash
python Feature4_AlphaFold.py --limit-genes 100
python Feature13_GWAS.py --max-rows 100000
python Feature18_MGI.py --max-pheno-rows 200000
python Feature22_Paralogues.py --limit-genes 1000
python Analysis3-MakeDatasets.py --dry-run
```

## Data Availability

All source data used by the pipeline come from publicly available resources or
resources available for academic use. No private or patient-identifiable data
are used.

Repository:
[MuhammadMuneeb007/Gene_Druggability_Finder](https://github.com/MuhammadMuneeb007/Gene_Druggability_Finder)

## Citation and Reuse
MIT LICENSE
The repository currently does not contain a licence file. Verify reuse and
redistribution permissions before using the code outside the project.
