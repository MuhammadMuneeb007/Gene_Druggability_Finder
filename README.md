# Gene Druggability Finder (GDF)

**A comprehensive computational pipeline for assessing small molecule and biologic druggability of protein targets**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## 📋 Table of Contents
- [Overview](#overview)
- [Pipeline Workflow](#pipeline-workflow)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Detailed Step-by-Step Guide](#detailed-step-by-step-guide)
- [Composite Scoring System](#composite-scoring-system)
- [Performance Metrics](#performance-metrics)
- [Limitations](#limitations)
- [Future Improvements](#future-improvements)
- [Citation](#citation)

---

## 🎯 Overview

**Gene Druggability Finder (GDF)** is a 12-step computational pipeline that predicts protein druggability using dual-track scoring:
- **Small Molecule Track**: Assesses oral drug potential (pills, tablets)
- **Biologic Track**: Evaluates antibody/protein therapeutic potential

### Key Features
✅ Integrates **10+ authoritative databases** (UniProt, PDB, AlphaFold, GTEx, STRING, DGIdb, ChEMBL, Pharos, etc.)  
✅ **Dual-track scoring** distinguishes small molecule vs biologic druggability  
✅ **Validated performance**: 93% accuracy, 0.933 F1-score  
✅ **Comprehensive analysis**: Structure, expression, disorder, networks, and more  
✅ **Automated workflow**: Single command execution  

---

## 📊 Pipeline Workflow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            INPUT: Gene List                                  │
│                        (User-defined in main.py)                             │
│                        Example: ADRB2, JAK2, TNF, TP53                       │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 1: Input Preparation                                                  │
│  ├─ Query: MyGene.info, UniProt                                             │
│  ├─ Fields Extracted:                                                       │
│  │   • gene_symbol, ensembl_gene_id, entrez_gene_id                         │
│  │   • uniprot_swissprot, uniprot_trembl, protein_name                      │
│  │   • chromosome, start_pos, end_pos, strand                               │
│  │   • pdb_ids, hgnc_id, vega_id, refseq_ids                                │
│  │   • gene_description, gene_type, gene_aliases                            │
│  │   • 50+ database cross-references                                        │
│  └─ Output: Step1.csv (120 genes × 50+ columns)                             │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 2: Known Drug Interactions                                            │
│  ├─ Query: DGIdb, ChEMBL, Open Targets                                      │
│  ├─ Fields Extracted:                                                       │
│  │   • dgidb_num_approved, dgidb_approved_drugs (drug names)                │
│  │   • dgidb_num_clinical, dgidb_clinical_drugs                             │
│  │   • chembl_approved_drugs, chembl_clinical_drugs                         │
│  │   • chembl_bioactive_compounds (count)                                   │
│  │   • opentargets_approved, opentargets_clinical                           │
│  │   • has_approved_drugs (Boolean), druggability_evidence (HIGH/MED/LOW)   │
│  │   • interaction_types, source_databases                                  │
│  └─ Output: Step2.csv (adds 15+ drug interaction columns)                   │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 3: Protein Family Classification                                      │
│  ├─ Query: InterPro, Pfam, UniProt                                          │
│  ├─ Fields Extracted:                                                       │
│  │   • druggable_family_type (GPCR, Kinase, Ion Channel, etc.)              │
│  │   • is_druggable_family (Boolean)                                        │
│  │   • druggable_family_confidence (HIGH, MEDIUM, LOW)                      │
│  │   • pfam_domains (domain IDs), interpro_domains                          │
│  │   • protein_class (Enzyme, Receptor, Transporter, etc.)                  │
│  │   • go_molecular_function, go_biological_process                         │
│  │   • enzyme_classification (EC numbers)                                   │
│  └─ Output: Step3.csv (adds 12+ family/domain columns)                      │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 4: Protein Structure Retrieval                                        │
│  ├─ Query: PDB, AlphaFold                                                   │
│  ├─ Fields Extracted:                                                       │
│  │   • has_pdb_structure (Boolean), num_pdb_structures                      │
│  │   • best_pdb_id, best_pdb_resolution (Å), pdb_method                     │
│  │   • pdb_all_ids (list), pdb_release_dates                                │
│  │   • has_alphafold_structure (Boolean), alphafold_id                      │
│  │   • alphafold_avg_plddt (0-100 confidence), alphafold_url                │
│  │   • structure_source (PDB/AlphaFold), structure_file (path)              │
│  │   • sequence_length, molecular_weight                                    │
│  └─ Output: Step4.csv (adds 15+ structure columns) + .pdb files             │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 5: Binding Pocket Prediction                                          │
│  ├─ Tool: Fpocket (Voronoi tessellation-based pocket detection)             │
│  ├─ Fields Extracted:                                                       │
│  │   • pockets_found (count), has_druggable_pocket (Boolean)                │
│  │   • best_pocket_volume (Å³), best_pocket_druggability (0-1 score)        │
│  │   • best_pocket_hydrophobicity, best_pocket_polarity_score               │
│  │   • best_pocket_charge, best_pocket_alpha_sphere_density                 │
│  │   • best_pocket_mean_alpha_sphere_radius                                 │
│  │   • best_pocket_category (Excellent/Good/Moderate/Poor)                  │
│  │   • favorable_features (list), concerns (list)                           │
│  │   • all_pocket_scores (top 5 pockets)                                    │
│  └─ Output: Step5.csv (adds 18+ pocket columns) + pocket files              │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 6: Subcellular Location Assessment                                    │
│  ├─ Query: UniProt (subcellular_location, topology fields)                  │
│  ├─ Fields Extracted:                                                       │
│  │   • Step6_Subcellular_Location (Membrane, Cytoplasm, Nucleus, etc.)      │
│  │   • Step6_Topology (Extracellular, Transmembrane, Intracellular)         │
│  │   • Step6_Accessibility_Score (0-4, higher = more accessible)            │
│  │   • Step6_Druggability_Tier (Tier 1-4)                                   │
│  │   • Step6_Location_Category (Secreted, Membrane, Intracellular)          │
│  │   • Step6_Drug_Modality_Options (Small Molecule, Antibody, Both)         │
│  │   • Step6_Targeting_Difficulty (Easy, Moderate, Difficult, Very Diff.)   │
│  │   • Step6_Transmembrane_Domains (count), signal_peptide (Boolean)        │
│  └─ Output: Step6.csv (adds 10+ location columns)                           │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 7: Tissue Expression Analysis                                         │
│  ├─ Query: GTEx v8 (54 tissues)                                             │
│  ├─ Fields Extracted:                                                       │
│  │   • Step7_Max_TPM (maximum expression in any tissue)                     │
│  │   • Step7_Mean_TPM, Step7_Median_TPM, Step7_Std_TPM                      │
│  │   • Step7_Tissues_With_Expression (count ≥0.1 TPM)                       │
│  │   • Step7_Tissue_Specificity (Highly/Moderately/Broadly/Ubiquitous)      │
│  │   • Step7_Top_3_Tissues (highest expressing tissues)                     │
│  │   • Step7_Brain_Expression (Boolean), Step7_Blood_Expression             │
│  │   • Step7_Expression_In_[54_tissues] (individual TPM values)             │
│  │   • Step7_Overall_Expression_Score (0-10)                                │
│  └─ Output: Step7.csv (adds 60+ expression columns)                         │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 8: Protein Disorder Analysis                                          │
│  ├─ Query: MobiDB, AlphaFold pLDDT, Composition-based prediction            │
│  ├─ Fields Extracted:                                                       │
│  │   • Step8_Disorder_Percentage (% of residues disordered)                 │
│  │   • Step8_Passes_50_Percent_Rule (Boolean, <50% = PASS)                  │
│  │   • Step8_Disorder_Category (Highly Ordered/Ordered/Disordered/Highly)   │
│  │   • Step8_Disorder_Method (MobiDB/AlphaFold/Composition/Unavailable)     │
│  │   • Step8_Druggability_Impact (Excellent/Good/Acceptable/Challenging)    │
│  │   • Step8_Num_Disordered_Regions (count), disorder_region_coords         │
│  │   • Step8_Longest_Ordered_Region (residue count)                         │
│  │   • Step8_Confidence (High/Medium/Low based on method)                   │
│  └─ Output: Step8.csv (adds 12+ disorder columns)                           │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 9: Druggable Genome Classification                                    │
│  ├─ Query: Pharos/TCRD GraphQL API                                          │
│  ├─ Fields Extracted:                                                       │
│  │   • Step9_Target_Development_Level (Tclin/Tchem/Tbio/Tdark)              │
│  │   • Step9_TDL_Description (explanation of category)                      │
│  │   • Step9_Protein_Family (Pharos classification)                         │
│  │   • Step9_Novelty_Score (0-1, where 0=well-studied, 1=novel)             │
│  │   • Step9_PubMed_Count (number of publications)                          │
│  │   • Step9_Patent_Count, Step9_Grant_Count                                │
│  │   • Step9_Development_Status (Clinical/Pre-clinical/Research/Unknown)    │
│  │   • Step9_Overall_Pharos_Score (0-10 composite)                          │
│  └─ Output: Step9.csv (adds 10+ Pharos columns)                             │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 10: Network Analysis                                                  │
│  ├─ Query: STRING (protein-protein interaction database)                    │
│  ├─ Fields Extracted:                                                       │
│  │   • Step10_Total_Interactions (count of all partners)                    │
│  │   • Step10_Drug_Target_Interactions (partners that are drug targets)     │
│  │   • Step10_Network_Hub (Boolean, >15 interactions)                       │
│  │   • Step10_Network_Connectivity (Highly/Well/Moderately/Poorly)          │
│  │   • Step10_Top_Interactors (list of top 10 partners)                     │
│  │   • Step10_Drug_Target_Families (families of interacting drug targets)   │
│  │   • Step10_Average_Interaction_Score (0-1000, STRING confidence)         │
│  │   • Step10_Pathways (KEGG, Reactome pathways)                            │
│  │   • Step10_GO_Processes (shared biological processes)                    │
│  └─ Output: Step10.csv (adds 12+ network columns)                           │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 11: Composite Druggability Score ★★★                                  │
│  ├─ Integration: Weighted scoring from ALL prior steps                      │
│  ├─ Fields Extracted:                                                       │
│  │   • Step11_SmallMolecule_Score (0-100, weighted 29-point system)         │
│  │   • Step11_Biologic_Score (0-100, weighted 20-point system)              │
│  │   • Step11_Composite_Score (max of SM or Bio score)                      │
│  │   • Step11_Best_Modality (Small Molecule/Biologic/Both/Neither)          │
│  │   • Step11_Classification (High/Medium/Low Druggability)                 │
│  │   • Step11_Priority (High/Medium/Low Priority for development)           │
│  │   • Step11_Is_Biologic_Only (Boolean, cytokines/secreted only)           │
│  │   • Step11_Score_Breakdown (individual component scores)                 │
│  │   • Step11_Confidence (High/Medium/Low based on data completeness)       │
│  │   • Step11_Rationale (text explanation of score)                         │
│  └─ Output: Step11.csv (adds 15+ scoring columns)                           │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 12: Validation & Performance Metrics                                  │
│  ├─ Query: PubMed, ClinicalTrials.gov                                       │
│  ├─ Fields Extracted:                                                       │
│  │   • Step12_Literature_Support (PubMed articles mentioning druggability)  │
│  │   • Step12_Clinical_Trials (count of active trials)                      │
│  │   • Step12_Trial_Phases (Phase I/II/III/IV breakdown)                    │
│  │   • Step12_Validation_Status (Validated/Partially/Not Validated)         │
│  │   • Step12_Known_Drugs_Literature (drugs mentioned in papers)            │
│  │   • Step12_Predicted_vs_Actual (comparison with ground truth)            │
│  │   • Step12_Confidence_Score (validation confidence)                      │
│  └─ Output: Step12.csv (adds 10+ validation columns)                        │
│  └─ Also Creates:                                                           │
│      • Step12_final_predictions.csv (high-priority targets)                 │
│      • Step12_performance_metrics.txt (AUC, Accuracy, F1, etc.)             │
│      • Step12_confusion_matrix.txt (TP/FP/TN/FN breakdown)                  │
└──────────────────────────────┬──────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           FINAL OUTPUT FILES                                 │
│                                                                              │
│  📊 Step12.csv - Complete dataset (120 genes × 250+ columns)                │
│      ├─ All original gene identifiers (Step 1)                              │
│      ├─ Drug interaction data (Step 2)                                      │
│      ├─ Protein family classifications (Step 3)                             │
│      ├─ Structure quality metrics (Step 4)                                  │
│      ├─ Binding pocket properties (Step 5)                                  │
│      ├─ Subcellular accessibility (Step 6)                                  │
│      ├─ Tissue expression profiles (Step 7)                                 │
│      ├─ Disorder analysis (Step 8)                                          │
│      ├─ Pharos TDL classifications (Step 9)                                 │
│      ├─ Network connectivity (Step 10)                                      │
│      ├─ Druggability scores (Step 11) ★                                     │
│      └─ Validation metrics (Step 12)                                        │
│                                                                              │
│  🎯 Step12_final_predictions.csv - Top candidates sorted by score           │
│                                                                              │
│  📈 Step12_performance_metrics.txt - Model performance                      │
│      ├─ Accuracy: 93% (93/100 correct)                                      │
│      ├─ Sensitivity: 96% (48/50 druggable caught)                           │
│      ├─ Specificity: 90% (45/50 non-druggable excluded)                     │
│      ├─ F1-Score: 0.933                                                     │
│      ├─ AUC-ROC: [calculated from scores]                                   │
│      └─ Confusion matrix: TP=48, FP=5, FN=2, TN=45                          │
│                                                                              │
│  📁 Additional Outputs:                                                     │
│      ├─ protein_structures/ - Downloaded PDB/AlphaFold files                │
│      ├─ binding_pockets/ - Fpocket output for each protein                  │
│      └─ gene_list.json - Input gene list                                    │
└─────────────────────────────────────────────────────────────────────────────┘

TOTAL DATA GENERATED: ~120 genes × 250+ features = 30,000+ data points
```

---

## 🚀 Installation

### Prerequisites
- **Python 3.8+**
- **fpocket** (for binding pocket analysis)
- Internet connection (for API queries)

### Step-by-Step Installation

```bash
# 1. Clone or download the repository
cd DrugFinalCode

# 2. Install Python dependencies
pip install -r requirements.txt

# 3. Install fpocket (Linux/Mac)
# Download from: https://github.com/Discngine/fpocket
# Follow fpocket installation instructions (see README.md)
# Verify installation:
fpocket -h

# 4. Verify Python packages
python -c "import pandas, numpy, scipy, sklearn, requests; print('✓ All packages installed')"
```

---

## ⚡ Quick Start

### 1. Define Your Gene List

Edit `main.py` and modify the `GENE_LIST`:

```python
GENE_LIST = [
    'ADRB2',   # Beta-2 adrenergic receptor
    'JAK2',    # JAK2 kinase
    'TNF',     # TNF cytokine
    'TP53',    # p53 tumor suppressor
    # Add your genes here
]
```

### 2. Run the Complete Pipeline

```bash
python main.py
```

### 3. View Results

- **Step12.csv** - Complete dataset with all features
- **Step12_final_predictions.csv** - Prioritized drug targets
- **Step12_performance_metrics.txt** - Validation metrics

---

## 📖 Detailed Step-by-Step Guide

### STEP 1: Input Preparation

**Python File:** `step01_input_preparation.py`

**Execution:**
```bash
# Automatically executed by main.py
# Or run standalone (requires gene_list.json from main.py):
python step01_input_preparation.py
```

**What It Does:**
- Queries **MyGene.info** for comprehensive gene identifiers
- Queries **UniProt** for protein-level annotations
- Retrieves 50+ identifier types per gene

**Output:** `Step1.csv`

**Key Columns:**
| Column | Description | Why Important |
|--------|-------------|---------------|
| `gene_symbol` | Official gene symbol | Primary identifier |
| `ensembl_gene_id` | Ensembl ID | For GTEx queries |
| `uniprot_swissprot` | UniProt reviewed ID | Highest quality annotation |
| `entrez_gene_id` | NCBI Gene ID | Cross-database linking |
| `pdb_ids` | PDB structure IDs | Indicates structural data |
| `gene_description` | Functional summary | Context understanding |

**Why Important:**
- Comprehensive IDs enable querying diverse databases
- Prevents lookup failures due to ID mismatches
- Reviewed UniProt IDs ensure data quality

**How to Improve:**
- Add isoform-specific analysis
- Include protein complex information
- Add disease-gene associations

**Limitations:**
- Uses canonical sequences only (no splice variants)
- Limited to human proteins (organism: 9606)
- API rate limits may slow large gene lists

**Sample Output (Step1.csv - First 3 genes):**

```
gene_symbol  ensembl_gene_id  uniprot_swissprot  entrez_gene_id  chromosome  protein_length  pdb_ids
ADRB2        ENSG00000169252  P07550             154             5           413             1GQ4|2R4R|2R4S|2RH1|...
ADRB1        ENSG00000043591  P08588             153             10          477             2LSQ|7BTS|7BU6|7BU7|...
ADRA1A       ENSG00000171873  P25100             146             20          572             (None)

gene_description                                                      subcellular_location
Beta-2-adrenergic receptor, GPCR superfamily, calcium channel...     Cell membrane; Early endosome; Golgi
Alpha-1-adrenergic receptor, Gq/11 G-protein coupled...              Cell membrane; Early endosome  
Alpha-1D-adrenergic receptor, activates mitogenic responses...       Cell membrane
```

---

### STEP 2: Known Drug Interactions

**Python File:** `step02_known_drug_interactions.py`

**Execution:**
```bash
python step02_known_drug_interactions.py
```

**What It Does:**
- Queries **DGIdb** for drug-gene interactions
- Queries **ChEMBL** for bioactive compounds
- Queries **Open Targets** for clinical drugs
- Identifies approved drugs, clinical candidates, and tool compounds

**Output:** `Step2.csv` (all Step1 columns + drug data)

**Key Columns:**
| Column | Description | Why Important |
|--------|-------------|---------------|
| `dgidb_num_approved` | Number of approved drugs | Proven druggability |
| `dgidb_approved_drugs` | Drug names | Identifies existing therapies |
| `chembl_approved_drugs` | ChEMBL approved count | Independent validation |
| `opentargets_approved` | Open Targets approved | Clinical precedent |
| `has_approved_drugs` | Boolean flag | Quick filter |
| `chembl_bioactive_compounds` | Bioactive molecule count | Tool compound availability |
| `druggability_evidence` | Evidence level | Overall assessment |

**Why Important:**
- Existing drugs prove target is druggable
- Clinical compounds indicate feasibility
- Bioactive compounds enable target validation

**How to Improve:**
- Add drug-target binding affinity data
- Include failed clinical trials (safety flags)
- Add patent information

**Limitations:**
- Drug databases may have annotation delays
- Small molecule bias (biologics underrepresented in DGIdb)
- Does not distinguish direct vs indirect effects

**Sample Output (Step2.csv - First 3 genes, drug columns):**

```
gene_symbol  dgidb_num_approved  dgidb_approved_drugs                                                            has_approved_drugs
ADRB2        85                  CARVEDILOL|SODIUM CHLORIDE|PENBUTOLOL SULFATE|DESIPRAMINE|METIPRANOLOL|...      True
ADRB1        69                  BRETYLIUM|PINDOLOL|BISOPROLOL|FLECAINIDE ACETATE|METIPRANOLOL|...               True
ADRA1A       81                  PROPIOMAZINE|METHOTRIMEPRAZINE|METHOXAMINE|EPHEDRINE|ALFUZOSIN|...              True
```

*Note: These GPCRs have 69-85 approved drugs each, demonstrating proven druggability.*

---

### STEP 3: Protein Family Classification

**Python File:** `step03_protein_family_classification.py`

**Execution:**
```bash
python step03_protein_family_classification.py
```

**What It Does:**
- Queries **InterPro** for domain architecture
- Queries **Pfam** for protein families
- Classifies into 30+ druggable families
- Assigns confidence levels (HIGH, MEDIUM, LOW)

**Output:** `Step3.csv`

**Key Columns:**
| Column | Description | Why Important |
|--------|-------------|---------------|
| `druggable_family_type` | Family classification | Predicts druggability |
| `is_druggable_family` | Boolean flag | Quick filter |
| `druggable_family_confidence` | Confidence level | Reliability assessment |
| `pfam_domains` | Pfam domain IDs | Functional modules |
| `interpro_domains` | InterPro annotations | Comprehensive classification |
| `protein_class` | Broad protein class | High-level categorization |

**Druggable Families Identified:**
- **Classic:** GPCRs, Kinases, Ion Channels, Nuclear Receptors, Proteases
- **Emerging:** E3 Ligases, Bromodomains, HDACs, Methyltransferases
- **Structural:** Tubulin, Actin (for specific applications)

**Why Important:**
- Protein family is strongest druggability predictor
- GPCRs and kinases account for 50%+ of drugs
- Guides drug discovery strategy

**How to Improve:**
- Add allosteric site prediction
- Include post-translational modifications
- Add protein-protein interaction interface analysis

**Limitations:**
- Domain-based classification may miss atypical targets
- Novel families (not in databases) scored as non-druggable
- Does not distinguish isoforms

**Sample Output (Step3.csv - First 3 genes, family columns):**

```
gene_symbol  druggable_family_type  is_druggable_family  druggable_family_confidence  protein_class
ADRB2        GPCR                   True                 HIGH                         GPCR
ADRB1        GPCR                   True                 HIGH                         GPCR
ADRA1A       GPCR                   True                 HIGH                         GPCR
```

*Note: All three are classified as GPCRs with HIGH confidence - the gold standard for druggability.*

---

### STEP 4: Protein Structure Retrieval

**Python File:** `step04_protein_structure.py`

**Execution:**
```bash
python step04_protein_structure.py
```

**What It Does:**
- Searches **PDB** for experimental structures
- Queries **AlphaFold** for predicted structures
- Downloads structure files (.pdb format)
- Prioritizes experimental > predicted structures

**Output:** `Step4.csv` + structure files in `protein_structures/`

**Key Columns:**
| Column | Description | Why Important |
|--------|-------------|---------------|
| `has_pdb_structure` | Experimental structure exists | Highest confidence |
| `num_pdb_structures` | Count of PDB entries | Multiple structures available |
| `best_pdb_resolution` | Resolution (Å) | Quality metric |
| `has_alphafold_structure` | Predicted structure exists | Backup for no PDB |
| `alphafold_avg_plddt` | Confidence score (0-100) | Prediction quality |
| `structure_source` | PDB or AlphaFold | Data provenance |
| `structure_file` | File path | For pocket analysis |

**Why Important:**
- Structure enables binding site identification
- Resolution <2.5Å is ideal for drug design
- AlphaFold provides structures for 98%+ of proteome

**How to Improve:**
- Add cryo-EM structures
- Include structure quality validation (Ramachandran, etc.)
- Add conformational ensemble analysis

**Limitations:**
- AlphaFold struggles with disordered regions
- Membrane protein structures often incomplete
- Predicted structures lack ligand-bound conformations

**Sample Output (Step4.csv - First 3 genes, structure columns):**

```
gene_symbol  has_pdb_structure  best_pdb_resolution  has_alphafold_structure  alphafold_avg_plddt  structure_source
ADRB2        True               1.9                  True                     NaN                  PDB
ADRB1        True               NaN                  True                     NaN                  PDB
ADRA1A       False              NaN                  True                     NaN                  AlphaFold
```

*Note: ADRB2 has excellent experimental structure at 1.9Å resolution. ADRA1A uses AlphaFold prediction.*

---

### STEP 5: Binding Pocket Prediction

**Python File:** `step05_find_binding_pockets.py`

**Execution:**
```bash
python step05_find_binding_pockets.py
```

**What It Does:**
- Runs **Fpocket** on all protein structures
- Identifies potential binding pockets using alpha spheres
- Calculates pocket properties (volume, druggability, hydrophobicity)
- Ranks pockets by druggability score

**Output:** `Step5.csv` + pocket files in `binding_pockets/`

**Key Columns:**
| Column | Description | Why Important |
|--------|-------------|---------------|
| `pockets_found` | Number of pockets detected | Binding site availability |
| `has_druggable_pocket` | At least one good pocket | Druggability indicator |
| `best_pocket_volume` | Volume (Å³) | Size for small molecules |
| `best_pocket_druggability` | Fpocket score (0-1) | Binding likelihood |
| `best_pocket_hydrophobicity` | Hydrophobic score | Compound binding |
| `best_pocket_category` | Excellent/Good/Moderate | Overall quality |
| `favorable_features` | Positive features | Strengths |
| `concerns` | Negative features | Weaknesses |

**Druggability Criteria:**
- **Volume:** ≥500 Å³ optimal, ≥200 Å³ minimum
- **Druggability Score:** ≥0.7 good, ≥0.5 acceptable
- **Hydrophobicity:** >0.5 favorable

**Why Important:**
- Binding pockets are essential for small molecules
- Fpocket druggability score correlates with success
- Pocket volume determines compound size

**How to Improve:**
- Add pocket flexibility analysis (molecular dynamics)
- Include allosteric pocket detection
- Add fragment screening predictions

**Limitations:**
- Static structures (no dynamics)
- May miss cryptic pockets
- Fpocket optimized for small molecules (not biologics)

**Sample Output (Step5.csv - First 3 genes, pocket columns):**

```
gene_symbol  pockets_found  has_druggable_pocket  best_pocket_volume  best_pocket_druggability  best_pocket_category
ADRB2        7              True                  389.616             0.095                     Moderate
ADRB1        1              True                  388.968             0.157                     Moderate
ADRA1A       0              False                 NaN                 NaN                       NaN
```

*Note: ADRB2 has 7 pockets detected. Volume ~390Å³ is adequate for small molecules. ADRA1A has no structure available for pocket analysis.*

---

### STEP 6: Subcellular Location Assessment

**Python File:** `step06_subcellular_location.py`

**Execution:**
```bash
python step06_subcellular_location.py
```

**What It Does:**
- Queries **UniProt** for subcellular localization
- Analyzes topology (extracellular, transmembrane, intracellular)
- Assigns accessibility tiers (1-4)
- Determines drug modality feasibility

**Output:** `Step6.csv`

**Key Columns:**
| Column | Description | Why Important |
|--------|-------------|---------------|
| `Step6_Subcellular_Location` | Cellular compartment | Determines accessibility |
| `Step6_Topology` | Membrane topology | Drug access route |
| `Step6_Accessibility_Score` | Score 0-4 | Higher = easier access |
| `Step6_Druggability_Tier` | Tier 1-4 classification | Quick assessment |
| `Step6_Location_Category` | Categorized location | Grouped analysis |
| `Step6_Drug_Modality_Options` | Small molecule/Antibody | Guides strategy |
| `Step6_Targeting_Difficulty` | Easy/Moderate/Difficult | Development challenge |

**Tier System:**
- **Tier 1 (Score 3-4):** Extracellular, secreted, membrane with ECD → Small molecules + Antibodies
- **Tier 2 (Score 3):** Membrane (intramembrane) → Small molecules
- **Tier 3 (Score 2):** Cytoplasmic → Small molecules (cell-permeable)
- **Tier 4 (Score 1):** Nuclear, organellar → Challenging

**Why Important:**
- Location determines drug delivery feasibility
- Antibodies cannot access intracellular targets
- Membrane proteins are 60% of drug targets

**How to Improve:**
- Add blood-brain barrier (BBB) permeability prediction
- Include lysosomal targeting assessment
- Add mitochondrial accessibility analysis

**Limitations:**
- Does not predict BBB penetration for CNS targets
- Dynamic relocalization not captured
- Post-translational modifications affect location

**Sample Output (Step6.csv - First 3 genes, location columns):**

```
gene_symbol  Step6_Subcellular_Location                      Step6_Topology                                       Step6_Accessibility_Score  Step6_Druggability_Tier  Step6_Drug_Modality_Options
ADRB2        Cell membrane; Early endosome; Golgi apparatus  Multi-pass membrane protein; Extracellular; Cytoplasmic  3                     Tier 1 - High            Small molecule; Antibody (for extracellular epitopes)
ADRB1        Cell membrane; Early endosome                   Multi-pass membrane protein; Extracellular; Cytoplasmic  3                     Tier 1 - High            Small molecule; Antibody (for extracellular epitopes)
ADRA1A       Cell membrane                                   Multi-pass membrane protein; Extracellular; Cytoplasmic  3                     Tier 1 - High            Small molecule; Antibody (for extracellular epitopes)
```

*Note: All three are Tier 1 (membrane proteins with extracellular domains) - accessible to both small molecules and antibodies.*

---

### STEP 7: Tissue Expression Analysis

**Python File:** `step07_tissue_expression.py`

**Execution:**
```bash
python step07_tissue_expression.py
```

**What It Does:**
- Queries **GTEx v8** for expression across 54 tissues
- Calculates tissue specificity (ubiquitous vs specific)
- Identifies top expressed tissues
- Assesses accessibility (blood, brain, systemic)

**Output:** `Step7.csv`

**Key Columns:**
| Column | Description | Why Important |
|--------|-------------|---------------|
| `Step7_Max_TPM` | Maximum expression (TPM) | Expression level |
| `Step7_Tissues_With_Expression` | Number of tissues (≥0.1 TPM) | Breadth |
| `Step7_Tissue_Specificity` | Highly/Moderately/Broadly/Ubiquitous | Selectivity |
| `Step7_Top_3_Tissues` | Highest expressing tissues | Target location |
| `Step7_Brain_Expression` | Expressed in brain? | CNS target |
| `Step7_Blood_Expression` | Expressed in blood? | Systemic accessible |
| `Step7_Overall_Expression_Score` | Score 0-10 | Combined metric |

**Expression Categories:**
- **Highly Specific:** ≤5 tissues (good for selectivity)
- **Moderately Specific:** 6-15 tissues
- **Broadly Expressed:** 16-30 tissues
- **Ubiquitous:** >30 tissues (off-target risk)

**Why Important:**
- Confirms target is expressed where disease occurs
- High expression enables robust validation
- Tissue specificity predicts off-target effects

**How to Improve:**
- Add disease-specific expression (TCGA, etc.)
- Include developmental stage expression
- Add single-cell resolution data

**Limitations:**
- GTEx uses healthy tissues (no disease state)
- Bulk RNA-seq (no cell-type specificity)
- Does not capture dynamic regulation

**Sample Output (Step7.csv - First 3 genes, expression columns):**

```
gene_symbol  Step7_Max_TPM  Step7_Tissues_With_Expression  Step7_Tissue_Specificity  Step7_Top_3_Tissues                                                      Step7_Brain_Expression
ADRB2        31.851         54                             Ubiquitously Expressed    Skin_Not_Sun_Exposed_Suprapubic:31.9; Skin_Sun_Exposed_Lower_leg:27.5; Breast_Mammary_Tissue:27.4   Yes
ADRB1        24.293         54                             Ubiquitously Expressed    Lung:24.3; Heart_Atrial_Appendage:14.6; Heart_Left_Ventricle:14.1                                   Yes
ADRA1A       14.138         54                             Ubiquitously Expressed    Liver:14.1; Lung:3.7; Spleen:3.3                                                                    Yes
```

*Note: All expressed in all 54 tissues (ubiquitous). Max TPM 14-32 indicates moderate-to-high expression levels.*

---

### STEP 8: Protein Disorder Analysis

**Python File:** `step08_protein_disorder.py`

**Execution:**
```bash
python step08_protein_disorder.py
```

**What It Does:**
- Uses **MobiDB** for consensus disorder predictions
- Analyzes **AlphaFold pLDDT** scores (<50 = disordered)
- Applies composition-based prediction as fallback
- Applies **50% disorder rule** for druggability

**Output:** `Step8.csv`

**Key Columns:**
| Column | Description | Why Important |
|--------|-------------|---------------|
| `Step8_Disorder_Percentage` | % disordered residues | Druggability predictor |
| `Step8_Passes_50_Percent_Rule` | <50% disorder | Critical threshold |
| `Step8_Disorder_Category` | Highly Ordered/Disordered | Classification |
| `Step8_Disorder_Method` | MobiDB/AlphaFold/Composition | Data source |
| `Step8_Druggability_Impact` | Excellent/Good/Challenging | Impact assessment |
| `Step8_Num_Disordered_Regions` | Count of regions | Fragmentation |

**50% Disorder Rule:**
- **<10% disorder:** Excellent (highly structured)
- **10-30% disorder:** Good (mostly structured)
- **30-50% disorder:** Acceptable (moderately disordered)
- **>50% disorder:** Poor (intrinsically disordered proteins)

**Why Important:**
- Disorder >50% indicates lack of stable binding sites
- Intrinsically disordered proteins (IDPs) resist traditional drugs
- Most successful drugs target ordered proteins

**How to Improve:**
- Add disorder-to-order transition prediction
- Include phase separation propensity
- Add short linear motif (SLiM) analysis

**Limitations:**
- Does not predict disorder-to-order transitions upon binding
- Disordered proteins can still be targeted (e.g., molecular glues)
- Composition-based prediction is approximate

**Sample Output (Step8.csv - First 3 genes, disorder columns):**

```
gene_symbol  Step8_Disorder_Percentage  Step8_Passes_50_Percent_Rule  Step8_Disorder_Category  Step8_Disorder_Method  Step8_Druggability_Impact
ADRB2        19.73                      True                          Mostly Ordered           AlphaFold pLDDT        Excellent - Well-structured
ADRB1        23.36                      True                          Mostly Ordered           AlphaFold pLDDT        Good - Druggable (<50% disorder)
ADRA1A       38.74                      True                          Moderately Disordered    AlphaFold pLDDT        Good - Druggable (<50% disorder)
```

*Note: All three pass the 50% disorder rule (19-39% disorder). GPCRs have structured transmembrane domains with flexible loops.*

---

### STEP 9: Druggable Genome Classification

**Python File:** `step09_druggable_genome.py`

**Execution:**
```bash
python step09_druggable_genome.py
```

**What It Does:**
- Queries **Pharos/TCRD** via GraphQL API
- Retrieves Target Development Level (TDL)
- Extracts protein family and novelty score
- Classifies development stage

**Output:** `Step9.csv`

**Key Columns:**
| Column | Description | Why Important |
|--------|-------------|---------------|
| `Step9_Target_Development_Level` | Tclin/Tchem/Tbio/Tdark | Development stage |
| `Step9_TDL_Description` | Category description | Interpretation |
| `Step9_Protein_Family` | Pharos family | Cross-validation |
| `Step9_Novelty_Score` | 0-1 (0=well-studied) | Literature coverage |
| `Step9_Development_Status` | Status description | Current stage |
| `Step9_Overall_Pharos_Score` | Score 0-10 | Combined metric |

**TDL Levels:**
- **Tclin:** Targets with approved drugs or clinical candidates (BEST)
- **Tchem:** Targets with bioactive compounds (IC50/Ki <30nM)
- **Tbio:** Targets with disease association but no compounds
- **Tdark:** Minimal knowledge, no compounds (WORST)

**Why Important:**
- TDL reflects community validation
- Tclin/Tchem targets have higher success rates
- Novelty score identifies understudied targets

**How to Improve:**
- Add patent landscape analysis
- Include target safety data (adverse effects)
- Add target-disease association strength

**Limitations:**
- TDL biased toward well-studied targets
- Does not capture failed drug programs
- Novelty score penalizes new targets

**Sample Output (Step9.csv - First 3 genes, Pharos columns):**

```
gene_symbol  Step9_Target_Development_Level  Step9_TDL_Description                                          Step9_Protein_Family  Step9_Novelty_Score  Step9_Development_Status
ADRB2        Tclin                           Clinical - Targets with approved drugs or clinical candidates  GPCR                  5                    Clinical stage or approved
ADRB1        Tclin                           Clinical - Targets with approved drugs or clinical candidates  GPCR                  5                    Clinical stage or approved
ADRA1A       Tclin                           Clinical - Targets with approved drugs or clinical candidates  GPCR                  5                    Clinical stage or approved
```

*Note: All three are Tclin (clinical targets) - the highest development level. Novelty score of 5 indicates moderate literature coverage.*

---

### STEP 10: Network Analysis

**Python File:** `step10_network_analysis.py`

**Execution:**
```bash
python step10_network_analysis.py
```

**What It Does:**
- Queries **STRING** for protein-protein interactions
- Identifies drug target connectivity
- Determines network hubs
- Retrieves pathway annotations

**Output:** `Step10.csv`

**Key Columns:**
| Column | Description | Why Important |
|--------|-------------|---------------|
| `Step10_Total_Interactions` | Interaction count | Connectivity |
| `Step10_Drug_Target_Interactions` | Connected to known targets | Proximity to drugs |
| `Step10_Network_Hub` | High connectivity (>15) | Hub status |
| `Step10_Network_Connectivity` | Highly/Well/Moderately/Poorly | Category |
| `Step10_Drug_Target_Families` | Connected target families | Combination potential |
| `Step10_Pathways` | KEGG/Reactome pathways | Biological context |

**Why Important:**
- Network hubs affect multiple pathways (broad impact)
- Connectivity to drug targets enables polypharmacology
- Pathway context guides combination therapy

**How to Improve:**
- Add tissue-specific interaction networks
- Include dynamic network analysis
- Add synthetic lethality prediction

**Limitations:**
- STRING includes computational predictions (not all validated)
- Does not capture context-dependent interactions
- Hub proteins may have off-target effects

**Sample Output (Step10.csv - First 3 genes, network columns):**

```
gene_symbol  Step10_Total_Interactions  Step10_Drug_Target_Interactions  Step10_Network_Hub  Step10_Network_Connectivity
ADRB2        20                         1                                True                Highly Connected (Hub)
ADRB1        20                         0                                True                Highly Connected (Hub)
ADRA1A       20                         0                                True                Highly Connected (Hub)
```

*Note: All three are network hubs with 20 interactions each. ADRB2 has 1 drug target interaction. High connectivity suggests broad pathway involvement.*

---

## 🎯 STEP 11: Composite Druggability Score ★

**Python File:** `step11_composite_score.py`

**Execution:**
```bash
python step11_composite_score.py
```

**What It Does:**
- **Integrates all 10 prior steps** into final scores
- **Dual-track scoring:** Small Molecule AND Biologic
- **Weighted multi-criteria** decision making
- **Identifies best modality** for each target

**Output:** `Step11.csv`

### Key Output Columns

| Column | Description | Range |
|--------|-------------|-------|
| `Step11_Composite_Score` | Best of SM or Bio score | 0-100 |
| `Step11_SmallMolecule_Score` | Small molecule druggability | 0-100 |
| `Step11_Biologic_Score` | Biologic druggability | 0-100 |
| `Step11_Best_Modality` | Small Molecule/Biologic/Both | Categorical |
| `Step11_Classification` | High/Medium/Low Druggability | Categorical |
| `Step11_Priority` | High/Medium/Low Priority | Categorical |
| `Step11_Is_Biologic_Only` | Cannot use small molecules | Boolean |

---

### Small Molecule Scoring System

**Total Weights: 29 points (normalized to 0-100)**

| Component | Weight | Max Points | % of Total | Description |
|-----------|--------|------------|------------|-------------|
| **Known Drugs** | 10 | 10 | 34.5% | Approved small molecule drugs |
| **Druggable Family** | 6 | 10 | 20.7% | GPCR/Kinase/IC/NR/Enzyme |
| **Binding Pocket** | 5 | 10 | 17.2% | Fpocket druggability score |
| **Pharos TDL** | 2 | 10 | 6.9% | Tclin/Tchem/Tbio/Tdark |
| **Accessibility** | 2 | 10 | 6.9% | Subcellular location tier |
| **Low Disorder** | 2 | 10 | 6.9% | <50% disorder rule |
| **Expression** | 1 | 10 | 3.4% | GTEx tissue expression |
| **Network** | 1 | 10 | 3.4% | STRING connectivity |

**Scoring Logic:**

#### 1. Known Drugs (10/29 weight)
- **10 points:** ≥10 approved small molecule drugs
- **9 points:** 5-9 approved drugs
- **8 points:** 3-4 approved drugs
- **7 points:** 1-2 approved drugs
- **4 points:** Clinical candidates
- **1-3 points:** Bioactive compounds (10-100+)
- **0 points:** No compounds

*Why weighted highest:* Proven druggability is strongest predictor

#### 2. Druggable Family (6/29 weight)
- **10 points:** GPCR or Kinase (gold standard)
- **9 points:** Ion Channel or Nuclear Receptor
- **8 points:** Protease (many drugs)
- **7 points:** Transporter or Enzyme (HIGH confidence)
- **5 points:** Enzyme (any confidence)
- **2-4 points:** Other druggable family
- **0 points:** Non-druggable family

*Why important:* Family membership predicts success rate

#### 3. Binding Pocket (5/29 weight)
- **10 points:** Druggability ≥0.8, Volume ≥500 Å³
- **8 points:** Druggability ≥0.6, Volume ≥400 Å³
- **6 points:** Druggability ≥0.5, Volume ≥300 Å³
- **3 points:** Druggability ≥0.3, Volume ≥200 Å³
- **0 points:** No pockets or poor quality

*Why important:* Binding site is necessary for small molecules

#### 4. Pharos TDL (2/29 weight)
- **10 points:** Tclin (clinical targets)
- **7 points:** Tchem (chemical targets)
- **3 points:** Tbio (biological targets)
- **0 points:** Tdark (understudied)

#### 5. Accessibility (2/29 weight)
- **10 points:** Tier 1 (Extracellular/Membrane+ECD)
- **7 points:** Tier 2 (Intramembrane)
- **4 points:** Tier 3 (Cytoplasmic)
- **1 point:** Tier 4 (Nuclear/Organellar)

#### 6. Low Disorder (2/29 weight)
- **10 points:** <10% disorder
- **8 points:** 10-20% disorder
- **6 points:** 20-30% disorder
- **4 points:** 30-40% disorder
- **2 points:** 40-50% disorder
- **0 points:** >50% disorder (fails 50% rule)

#### 7. Expression (1/29 weight)
- **10 points:** Expressed in ≥40 tissues, TPM ≥10
- **7 points:** Expressed in ≥30 tissues
- **5 points:** Expressed in ≥20 tissues
- **3 points:** Expressed in ≥10 tissues
- **0 points:** Poor expression

#### 8. Network (1/29 weight)
- **10 points:** ≥10 drug target interactions or hub with ≥5
- **7 points:** 5-9 drug target interactions
- **4 points:** 2-4 drug target interactions
- **2 points:** 1 drug target interaction
- **0 points:** No drug target connectivity

---

### Biologic Scoring System

**Total Weights: 20 points (normalized to 0-100)**

| Component | Weight | % of Total | Description |
|-----------|--------|------------|-------------|
| **Known Biologics** | 8 | 40% | Approved antibodies/proteins |
| **Secreted/Accessible** | 6 | 30% | Extracellular location |
| **Pharos TDL** | 3 | 15% | Development level |
| **Expression** | 2 | 10% | Tissue expression |
| **Network** | 1 | 5% | Connectivity |

**Scoring Logic:**

#### 1. Known Biologics (8/20 weight)
- **10 points:** ≥5 approved biologics
- **9 points:** 3-4 approved biologics
- **8 points:** 1-2 approved biologics
- **5 points:** Clinical biologics
- **0 points:** No biologics

#### 2. Secreted/Accessible (6/20 weight)
- **10 points:** Secreted protein
- **9 points:** Membrane + extracellular domain
- **6 points:** Membrane (any)
- **2 points:** Collagen (insoluble, not antibody-accessible)
- **1 point:** Intracellular

*Why important:* Antibodies cannot cross cell membranes

---

### Classification Thresholds

| Score Range | Classification | Priority | Typical Outcome |
|-------------|----------------|----------|-----------------|
| ≥80 | High Druggability | High Priority | Strong candidates |
| 50-79 | Medium Druggability | Medium Priority | Feasible targets |
| <50 | Low Druggability | Low Priority | Challenging targets |

---

### Example Scored Targets

#### High Small Molecule Score
```
HTR1B (Serotonin 1B Receptor)
├─ SM Score: 98.62/100
├─ Known Drugs: 10/10 (27 approved drugs)
├─ Family: 10/10 (GPCR)
├─ Pocket: 10/10 (excellent binding site)
├─ Classification: HIGH DRUGGABILITY
└─ Best Modality: Small Molecule
```

#### High Biologic Score
```
TNF (Tumor Necrosis Factor)
├─ Bio Score: 88.5/100
├─ SM Score: 6.9/100 (Biologic-only)
├─ Known Biologics: 10/10 (Humira, Enbrel, etc.)
├─ Secreted: 10/10 (fully accessible)
├─ Classification: HIGH DRUGGABILITY
└─ Best Modality: Biologic
```

#### Low Score (Undruggable)
```
TP53 (p53 Tumor Suppressor)
├─ SM Score: 72.07/100 → FALSE POSITIVE
├─ Family: 0/10 (Transcription factor)
├─ Pocket: 1/10 (poor pockets)
├─ Classification: LOW DRUGGABILITY
└─ Reality: Undruggable (40+ years of failure)
```

---

### Biologic-Only Flagging

**Genes Flagged as Biologic-Only:**
- Cytokines (TNF, IL1B, IL6, IL10, IL17A, IFNG)
- Growth factors (BDNF, NGF, VEGFA, TSLP)
- Collagens (COL1A1, COL1A2, COL3A1, COL4A1, COL5A1)
- Structural proteins (VIM, KRT1)

**Criteria:**
1. Gene symbol starts with IL/TNF/COL/BDNF
2. Protein family contains "cytokine", "growth factor", "collagen"
3. Secreted + no enzyme/protease activity
4. No small molecule drugs despite being Tier 1 accessible

**Why Important:**
- Prevents over-scoring cytokines for small molecules
- Reflects biological reality (antibodies work, pills don't)
- Guides drug discovery strategy

**Sample Output (Step11.csv - First 3 genes, scoring columns):**

```
gene_symbol  Step11_SmallMolecule_Score  Step11_Biologic_Score  Step11_Composite_Score  Step11_Best_Modality  Step11_Classification  Step11_Priority
ADRB2        80.34                       53.0                   80.34                   Small Molecule        High Druggability      High Priority
ADRB1        78.28                       52.0                   78.28                   Small Molecule        Medium Druggability    Medium Priority
ADRA1A       75.17                       52.0                   75.17                   Small Molecule        Medium Druggability    Medium Priority
```

*Note: ADRB2 scores 80.34/100 for small molecules (HIGH druggability). All three score ~52/100 for biologics (membrane accessibility). Small molecule is the best modality.*

---

## 📊 STEP 12: Validation & Performance

**Python File:** `step12_validate_candidates.py`

**Execution:**
```bash
python step12_validate_candidates.py
```

**What It Does:**
- Queries **PubMed** for literature validation
- Queries **ClinicalTrials.gov** for active trials
- Calculates **AUC, Accuracy, Precision, Recall**
- Generates performance metrics report

**Output:** 
- `Step12.csv` - Complete dataset with validation
- `Step12_final_predictions.csv` - Prioritized targets
- `Step12_performance_metrics.txt` - Detailed metrics

**Sample Output (Step12.csv - First 3 genes, final columns):**

```
gene_symbol  Step11_SmallMolecule_Score  Step11_Composite_Score  Step11_Classification  Step11_Best_Modality
ADRB2        80.34                       80.34                   High Druggability      Small Molecule
ADRB1        78.28                       78.28                   Medium Druggability    Small Molecule
ADRA1A       75.17                       75.17                   Medium Druggability    Small Molecule
```

*Note: Step12.csv contains all 250+ columns from previous steps. This shows the final druggability scores that are validated against literature and clinical data.*

---

## 🏆 Performance Metrics

### Validation Dataset
- **100 genes** (50 druggable, 50 non-druggable)
- **Druggable:** GPCRs, Kinases, Ion Channels, Nuclear Receptors, Enzymes
- **Non-druggable:** Cytokines, TFs, Structural, Adapters, GTPases, Ribosomal

### Classification Threshold
- **Druggable:** Small Molecule Score ≥70
- **Non-druggable:** Small Molecule Score <70

---

### Overall Performance

| Metric | Value | Interpretation |
|--------|-------|----------------|
| **Accuracy** | 93.0% (93/100) | Excellent overall performance |
| **Sensitivity** | 96.0% (48/50) | Catches 96% of druggable targets |
| **Specificity** | 90.0% (45/50) | Correctly excludes 90% of non-druggable |
| **Precision** | 90.6% (48/53) | 91% of predictions are correct |
| **F1-Score** | 0.933 | Excellent balance of precision/recall |

---

### Confusion Matrix

```
                           PREDICTED
                    Druggable (≥70)    Non-druggable (<70)
         ┌─────────────────────────────────────────────────┐
ACTUAL   │                                                 │
Druggable│     48 (TRUE POSITIVE)      2 (FALSE NEGATIVE) │ 50
         │                                                 │
Non-Drug │      5 (FALSE POSITIVE)    45 (TRUE NEGATIVE)  │ 50
         │                                                 │
         └─────────────────────────────────────────────────┘
              53 predicted              47 predicted
```

---

### Errors Analysis

#### ❌ 2 False Negatives (Druggable scored <70)

**Should be druggable but scored too low:**

1. **MET (c-MET Kinase)** - Score: 46.9 (Expected: >80)
   - **Reality:** FDA-approved drugs (crizotinib, capmatinib)
   - **Why missed:** Limited PDB structures, underweighted in DGIdb
   - **Impact:** CRITICAL miss (approved drug target)

2. **NR1I2 (PXR Nuclear Receptor)** - Score: 53.79 (Expected: >75)
   - **Reality:** Nuclear receptor, multiple binders
   - **Why missed:** Fewer approved drugs, nuclear location
   - **Impact:** Moderate (known druggable family)

#### ❌ 5 False Positives (Non-druggable scored ≥70)

**Should be undruggable but scored too high:**

1. **TP53 (p53)** - Score: 72.07 (Expected: <30)
   - **Reality:** Transcription factor, 40+ years of failure
   - **Why wrong:** Has clinical compounds (not FDA-approved)
   - **Impact:** MAJOR (classic undruggable target)

2. **JUN (c-Jun)** - Score: 70.34 (Expected: <30)
   - **Reality:** Transcription factor, no drugs
   - **Why wrong:** Some bioactive compounds reported
   - **Impact:** Moderate (on threshold)

3. **RIPK1** - Score: 73.10 (Expected: <30)
   - **Reality:** Difficult kinase target
   - **Why wrong:** Kinase family bias
   - **Impact:** Moderate (kinase but challenging)

4. **RAC1 (Small GTPase)** - Score: 76.55 (Expected: <30)
   - **Reality:** "Undruggable" small GTPase
   - **Why wrong:** Tool compounds exist
   - **Impact:** Moderate (emerging target area)

5. **RHOA (Small GTPase)** - Score: 78.62 (Expected: <30)
   - **Reality:** "Undruggable" small GTPase
   - **Why wrong:** Some bioactive compounds
   - **Impact:** Moderate (historically undruggable)

---

### Performance by Target Class

| Target Class | Accuracy | Notes |
|--------------|----------|-------|
| GPCRs | 100% (15/15) | Perfect classification |
| Kinases | 93% (14/15) | Missed MET |
| Ion Channels | 100% (8/8) | Perfect classification |
| Nuclear Receptors | 80% (4/5) | Missed NR1I2 |
| Enzymes | 100% (7/7) | Perfect classification |
| Cytokines | 100% (10/10) | Perfect biologic-only |
| Transcription Factors | 80% (8/10) | TP53, JUN false positives |
| Structural Proteins | 100% (10/10) | Perfect exclusion |
| GTPases | 50% (3/6) | RAC1, RHOA false positives |

---

### Key Insights

✅ **Strengths:**
- Excellent at identifying classic druggable families (GPCRs, Kinases, Ion Channels)
- Perfect biologic-only classification (cytokines)
- High sensitivity (catches 96% of druggable targets)

⚠️ **Weaknesses:**
- Small GTPases over-scored (emerging target class confusion)
- Transcription factors with tool compounds over-scored
- Some approved targets under-scored due to limited PDB data

---

## ⚠️ Limitations

### 1. **Known-Target Bias**
- **Issue:** 34.5% weight on existing drugs favors well-studied targets
- **Impact:** May under-score novel therapeutic targets
- **Mitigation:** Consider separate "novelty opportunity score"

### 2. **Isoform Specificity Not Assessed**
- **Issue:** Analyzes canonical sequences only
- **Impact:** Misses isoform-specific druggability (e.g., kinase splice variants)
- **Example:** PI3K has multiple isoforms with different druggability
- **Mitigation:** Add isoform-level analysis using Ensembl

### 3. **Post-Translational Modifications (PTMs) Ignored**
- **Issue:** Glycosylation, phosphorylation not considered
- **Impact:** Affects antibody epitope accessibility
- **Example:** Heavily glycosylated proteins harder to target
- **Mitigation:** Extract PTM features from UniProt

### 4. **Allosteric vs Orthosteric Sites Not Distinguished**
- **Issue:** Fpocket finds pockets but doesn't classify type
- **Impact:** May prioritize wrong binding site
- **Example:** Allosteric sites often better for selectivity
- **Mitigation:** Add allosteric site prediction tools

### 5. **Blood-Brain Barrier Not Considered**
- **Issue:** Brain expression ≠ CNS drug accessibility
- **Impact:** Over-scores CNS targets without BBB assessment
- **Example:** Brain target may need specialized delivery
- **Mitigation:** Add BBB permeability predictor

### 6. **Ubiquitous Expression Penalized**
- **Issue:** Housekeeping proteins scored lower on expression
- **Impact:** May under-score valid targets (e.g., COX, HMG-CoA reductase)
- **Reality:** Aspirin targets COX (ubiquitous), statins target HMGCR
- **Mitigation:** Reframe as "off-target risk" not "undruggability"

### 7. **Network Hub Bias**
- **Issue:** Highly connected proteins score higher
- **Impact:** Hubs may have more off-target effects
- **Example:** Hub disruption can cause toxicity
- **Mitigation:** Add selectivity challenge flag for hubs

### 8. **Disease Context Absent**
- **Issue:** Pipeline is disease-agnostic
- **Impact:** Doesn't assess disease-specific pathway importance
- **Example:** Target may be druggable but irrelevant for disease
- **Mitigation:** Add disease-gene association module

### 9. **Static Structures Only**
- **Issue:** No protein dynamics or flexibility
- **Impact:** Misses cryptic pockets that open upon binding
- **Example:** Some pockets only visible with molecular dynamics
- **Mitigation:** Add MD simulation or ensemble docking

### 10. **Small GTPase Classification Issues**
- **Issue:** Emerging target class confuses scoring
- **Impact:** KRAS was "undruggable" until 2021 (now druggable)
- **Evidence:** RAC1, RHOA scored 76-78 (false positives)
- **Mitigation:** Add special handling for emerging target classes

---

## 🚀 Future Improvements

### Short-Term (Weeks)

1. **Adjust Small GTPase Handling**
   - Add historical context flag
   - Distinguish KRAS (now druggable) from others
   - **Impact:** Fix 2 false positives

2. **Refine Transcription Factor Scoring**
   - Add TF-specific penalty
   - Distinguish DNA-binding (hard) from protein-binding (easier)
   - **Impact:** Fix TP53, JUN false positives

3. **Add MET Special Case**
   - Boost kinases with multiple approved drugs
   - **Impact:** Fix critical false negative

4. **Reweight Expression Component**
   - Reduce penalty for ubiquitous expression
   - Reframe as "off-target risk" metric
   - **Impact:** Better reflects biology

### Mid-Term (Months)

5. **Isoform-Level Analysis**
   - Query Ensembl for splice variants
   - Score each isoform separately
   - Flag isoform-specific drugs
   - **Databases:** Ensembl, APPRIS
   - **Impact:** Captures isoform biology

6. **Post-Translational Modifications**
   - Extract glycosylation sites from UniProt
   - Predict phosphorylation (NetPhos, etc.)
   - Flag heavily modified proteins
   - **Impact:** Better biologic predictions

7. **Blood-Brain Barrier Prediction**
   - Integrate BBBPredict or similar
   - Flag CNS targets requiring BBB penetration
   - Adjust scoring for brain-expressed targets
   - **Impact:** Realistic CNS druggability

8. **Allosteric Site Detection**
   - Add AlloPred or AlloSigMA
   - Distinguish orthosteric vs allosteric pockets
   - Prioritize allosteric for selectivity
   - **Impact:** Better pocket quality assessment

9. **Disease-Gene Association**
   - Integrate DisGeNET or OpenTargets disease links
   - Score pathway relevance
   - Prioritize disease-relevant targets
   - **Impact:** Therapeutic relevance

### Long-Term (Months+)

10. **Protein Dynamics**
    - Add molecular dynamics (GROMACS, AMBER)
    - Identify cryptic pockets
    - Assess binding site flexibility
    - **Impact:** Discover hidden sites

11. **Machine Learning Optimization**
    - Train on FDA-approved targets
    - Optimize weights via cross-validation
    - Add deep learning features
    - **Impact:** Improved accuracy

12. **Drug-Target Binding Affinity**
    - Add docking simulations (AutoDock, Glide)
    - Predict IC50 values
    - Rank compounds by affinity
    - **Impact:** Lead optimization guidance

13. **Safety Assessment**
    - Integrate adverse effect data (SIDER, FAERS)
    - Add toxicity predictions (ProTox, etc.)
    - Flag high-risk targets
    - **Impact:** Clinical safety

14. **Patent Landscape**
    - Search patent databases (USPTO, EPO)
    - Identify freedom-to-operate issues
    - Flag crowded target space
    - **Impact:** Commercial viability

---

## �‍💼 Author Information

- **Name**: Muhammad Muneeb
- **Affiliation**: The University of Queensland, Australia
- **Email**: [m.muneeb@uq.edu.au](mailto:m.muneeb@uq.edu.au)
- **Gmail**: [muneebsiddique007@gmail.com](mailto:muneebsiddique007@gmail.com)
- **GitHub**: [MuhammadMuneeb007](https://github.com/MuhammadMuneeb007/)
- **Google Scholar**: [Muhammad Muneeb](https://scholar.google.com/citations?hl=en&user=X0xdltIAAAAJ&view_op=list_works&sortby=pubdate)
- **ResearchGate**: [Muhammad Muneeb](https://www.researchgate.net/profile/Muhammad-Muneeb-5)
- **Supervisor**: [Prof. David Ascher](https://scmb.uq.edu.au/profile/8654/david-ascher)
- **Lab**: [BioSig Lab](https://biosig.lab.uq.edu.au/)

---

## �📚 Citation

If you use this pipeline in your research, please cite:

```
Muhammad Muneeb, David B. Ascher (2025). Gene Druggability Finder (GDF): A Comprehensive 
Dual-Track Computational Pipeline for Small Molecule and Biologic Target Assessment.
[Journal Name], [Volume], [Pages].
```

---

## 📞 Support

For questions, issues, or contributions:
- **Email**: [m.muneeb@uq.edu.au](mailto:m.muneeb@uq.edu.au) or [muneebsiddique007@gmail.com](mailto:muneebsiddique007@gmail.com)
- **GitHub**: [MuhammadMuneeb007](https://github.com/MuhammadMuneeb007/)
- **Lab Website**: [BioSig Lab](https://biosig.lab.uq.edu.edu/)

---

## 📄 License

This project is licensed under the MIT License - see LICENSE file for details.

---

## 🙏 Acknowledgments

### Databases Used
- **MyGene.info** - Gene annotation
- **UniProt** - Protein annotation
- **PDB** - Experimental structures
- **AlphaFold** - Predicted structures
- **DGIdb** - Drug-gene interactions
- **ChEMBL** - Bioactive compounds
- **Open Targets** - Clinical drugs
- **GTEx** - Tissue expression
- **STRING** - Protein interactions
- **Pharos/TCRD** - Target development levels
- **MobiDB** - Protein disorder

### Tools Used
- **Fpocket** - Binding pocket prediction
- **Python** - Data processing
- **pandas, numpy, scipy** - Numerical computing
- **scikit-learn** - Performance metrics

---

**Last Updated:** January 2025  
**Version:** 1.0  
**Pipeline Steps:** 12  
**Default Gene Set:** 120 genes  
**Estimated Runtime:** ~25-30 minutes for 120 genes
