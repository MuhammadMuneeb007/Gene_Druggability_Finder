"""
MAIN ORCHESTRATOR - Druggability Assessment Pipeline
=====================================================
Execute all 12 steps sequentially with progress tracking
User defines gene list here - pipeline processes automatically
"""

import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
import json


# ============================================================================
# USER-DEFINED GENE LIST - MODIFY THIS SECTION
# ============================================================================

GENE_LIST = [
    # GPCRs (15 genes) - Gold Standard
    'ADRB2',      # Beta-2 adrenergic receptor - 85+ drugs
    'ADRB1',      # Beta-1 adrenergic receptor - 50+ drugs
    'ADRA1A',     # Alpha-1A adrenergic - prazosin
    'DRD2',       # Dopamine D2 receptor - 131+ drugs
    'DRD3',       # Dopamine D3 receptor - antipsychotics
    'HTR1A',      # Serotonin 1A receptor - buspirone
    'HTR1B',      # Serotonin 1B receptor - 27 drugs
    'HTR2A',      # Serotonin 2A receptor - 69+ drugs
    'HTR2C',      # Serotonin 2C receptor - lorcaserin
    'CHRM1',      # Muscarinic M1
    'CHRM3',      # Muscarinic M3 - tiotropium
    'OPRM1',      # Mu opioid receptor - 76+ drugs
    'OPRD1',      # Delta opioid receptor
    'CNR1',       # Cannabinoid receptor 1
    'HRH1',       # Histamine H1 receptor - 100+ drugs
    
    # Kinases (15 genes) - Highly Druggable
    'ABL1',       # BCR-ABL kinase - imatinib
    'EGFR',       # EGFR kinase - gefitinib
    'BRAF',       # BRAF kinase - vemurafenib
    'ERBB2',      # HER2 kinase - lapatinib
    'ALK',        # ALK kinase - crizotinib
    'JAK2',       # JAK2 kinase - ruxolitinib
    'JAK1',       # JAK1 kinase - tofacitinib
    'BTK',        # Bruton's kinase - ibrutinib
    'SRC',        # SRC kinase - dasatinib
    'KIT',        # c-KIT kinase - imatinib
    'MET',        # MET kinase - crizotinib
    'RET',        # RET kinase - selpercatinib
    'FLT3',       # FLT3 kinase - midostaurin
    'CDK4',       # CDK4 kinase - palbociclib
    'MAPK1',      # ERK2 kinase - 30+ drugs
    
    # Ion Channels (8 genes)
    'SCN5A',      # Cardiac sodium channel - lidocaine
    'SCN9A',      # Nav1.7 sodium channel - 54+ drugs
    'KCNH2',      # hERG potassium channel
    'CACNA1C',    # L-type calcium channel - amlodipine
    'TRPV1',      # TRPV1 ion channel - capsaicin
    'KCNQ2',      # KCNQ2 potassium channel - retigabine
    'GABRA1',     # GABA-A receptor - benzodiazepines
    'GRIN1',      # NMDA receptor - ketamine
    
    # Nuclear Receptors (5 genes)
    'ESR1',       # Estrogen receptor - 109+ drugs
    'NR3C1',      # Glucocorticoid receptor - 50+ steroids
    'PPARG',      # PPAR-gamma - pioglitazone
    'THRB',       # Thyroid hormone receptor
    'NR1I2',      # PXR receptor
    
    # Enzymes (7 genes)
    'ACE',        # ACE enzyme - 42+ inhibitors
    'HMGCR',      # HMG-CoA reductase - statins
    'CYP2D6',     # Cytochrome P450 - 328+ drugs
    'PTGS2',      # COX-2 enzyme - celecoxib
    'PDE5A',      # PDE5 enzyme - sildenafil
    'ACHE',       # Acetylcholinesterase - donepezil
    'DPP4',       # DPP-4 enzyme - sitagliptin

    # Cytokines & Growth Factors (10 genes) - Biologic Only
    'TNF',        # TNF - ONLY antibodies (Humira)
    'IL1B',       # IL-1 beta - ONLY antibodies
    'IL6',        # IL-6 - ONLY antibodies (Tocilizumab)
    'IL10',       # IL-10 - antibodies only
    'IL17A',      # IL-17A - ONLY antibodies (Cosentyx)
    'TSLP',       # TSLP - antibody (Tezepelumab)
    'IFNG',       # Interferon gamma - biologic only
    'BDNF',       # BDNF - no SM drugs
    'NGF',        # NGF - antibody (tanezumab)
    'VEGFA',      # VEGF - antibodies (Avastin)
    
    # Transcription Factors (10 genes) - Undruggable
    'TP53',       # p53 - undruggable despite 40 years
    'MYC',        # c-Myc - undruggable TF
    'JUN',        # c-Jun - no direct inhibitors
    'FOS',        # c-Fos - no drugs
    'STAT3',      # STAT3 - no approved drugs
    'NFKB1',      # NF-kappa-B - indirect only
    'SOX2',       # SOX2 - stem cell TF
    'NANOG',      # Nanog - no drugs
    'POU5F1',     # OCT4 - undruggable
    'HIF1A',      # HIF-1 alpha - no direct drugs
    
    # Structural Proteins (10 genes) - Not Targetable
    'COL1A1',     # Collagen I alpha 1 - structural
    'COL1A2',     # Collagen I alpha 2 - structural
    'COL3A1',     # Collagen III alpha 1 - structural
    'COL4A1',     # Collagen IV alpha 1 - structural
    'COL5A1',     # Collagen V alpha 1 - structural
    'ACTB',       # Beta-actin - cytoskeleton
    'TUBB',       # Beta-tubulin - cytoskeleton
    'VIM',        # Vimentin - intermediate filament
    'LMNA',       # Lamin A/C - nuclear envelope
    'KRT1',       # Keratin 1 - structural
    
    # Adapter/Scaffold Proteins (8 genes) - No Binding Sites
    'GRB2',       # GRB2 - adapter protein
    'SHC1',       # SHC1 - adapter protein
    'GAB1',       # GAB1 - scaffold protein
    'IRS1',       # IRS1 - adapter protein
    'TRAF6',      # TRAF6 - scaffold protein
    'MYD88',      # MYD88 - adapter protein
    'TIRAP',      # TIRAP - adapter protein
    'RIPK1',      # RIPK1 - difficult target
    
    # Small GTPases (6 genes) - "Undruggable"
    'KRAS',       # KRAS - undruggable until 2021
    'NRAS',       # NRAS - no approved drugs
    'HRAS',       # HRAS - no approved drugs
    'RAC1',       # RAC1 - no drugs
    'CDC42',      # CDC42 - no drugs
    'RHOA',       # RhoA - no drugs
    
    # Ribosomal Proteins (6 genes) - Essential, Not Targetable
    'RPS6',       # Ribosomal protein S6
    'RPL5',       # Ribosomal protein L5
    'RPS3',       # Ribosomal protein S3
    'RPL10',      # Ribosomal protein L10
    'RPS19',      # Ribosomal protein S19
    'RPL11'       # Ribosomal protein L11
]


# ============================================================================
# PIPELINE CONFIGURATION
# ============================================================================

PIPELINE_STEPS = [
    {
        'step': 1,
        'name': 'Input Preparation',
        'script': 'step01_input_preparation.py',
        'description': 'Convert gene symbols to comprehensive IDs'
    },
    {
        'step': 2,
        'name': 'Known Drug Interactions',
        'script': 'step02_known_drug_interactions.py',
        'description': 'Query drug databases (DGIdb, ChEMBL, Open Targets)'
    },
    {
        'step': 3,
        'name': 'Protein Family Classification',
        'script': 'step03_protein_family_classification.py',
        'description': 'Classify into druggable families (GPCRs, Kinases, etc.)'
    },
    {
        'step': 4,
        'name': 'Protein Structure Retrieval',
        'script': 'step04_protein_structure.py',
        'description': 'Download structures from PDB and AlphaFold'
    },
    {
        'step': 5,
        'name': 'Binding Pocket Prediction',
        'script': 'step05_find_binding_pockets.py',
        'description': 'Run Fpocket to identify druggable binding sites'
    },
    {
        'step': 6,
        'name': 'Subcellular Location Analysis',
        'script': 'step06_subcellular_location.py',
        'description': 'Assess accessibility based on cellular location'
    },
    {
        'step': 7,
        'name': 'Tissue Expression Analysis',
        'script': 'step07_tissue_expression.py',
        'description': 'Query GTEx for tissue-specific expression'
    },
    {
        'step': 8,
        'name': 'Protein Disorder Analysis',
        'script': 'step08_protein_disorder.py',
        'description': 'Assess structural disorder (50% rule)'
    },
    {
        'step': 9,
        'name': 'Druggable Genome Classification',
        'script': 'step09_druggable_genome.py',
        'description': 'Query Pharos/TCRD for TDL classification'
    },
    {
        'step': 10,
        'name': 'Network Analysis',
        'script': 'step10_network_analysis.py',
        'description': 'Analyze protein-protein interactions (STRING)'
    },
    {
        'step': 11,
        'name': 'Composite Druggability Score',
        'script': 'step11_composite_score.py',
        'description': 'Calculate dual-track scores (SM + Biologic)'
    },
    {
        'step': 12,
        'name': 'Validate Candidates',
        'script': 'step12_validate_candidates.py',
        'description': 'Literature validation and performance metrics'
    }
]


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def print_banner():
    """Print pipeline banner"""
    print("\n" + "="*80)
    print(" DRUGGABILITY ASSESSMENT PIPELINE - COMPLETE WORKFLOW ".center(80))
    print("="*80)
    print(f" Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ".center(80))
    print(f" Total Genes: {len(GENE_LIST)} ".center(80))
    print(f" Total Steps: {len(PIPELINE_STEPS)} ".center(80))
    print("="*80 + "\n")


def print_step_header(step_info):
    """Print step header"""
    print("\n" + "─"*80)
    print(f"STEP {step_info['step']}: {step_info['name']}")
    print("─"*80)
    print(f"Description: {step_info['description']}")
    print(f"Script: {step_info['script']}")
    print(f"Time: {datetime.now().strftime('%H:%M:%S')}")
    print("─"*80 + "\n")


def print_step_complete(step_info, duration):
    """Print step completion"""
    print("\n" + "─"*80)
    print(f"✓ STEP {step_info['step']} COMPLETE - {step_info['name']}")
    print(f"Duration: {duration:.1f} seconds")
    print("─"*80 + "\n")


def save_gene_list(gene_list, output_file='gene_list.json'):
    """Save gene list to JSON for step01 to read"""
    with open(output_file, 'w') as f:
        json.dump({'genes': gene_list}, f, indent=2)
    print(f"✓ Gene list saved to {output_file}")


def run_step(step_info):
    """Execute a pipeline step"""
    print_step_header(step_info)
    
    start_time = time.time()
    
    try:
        # Run the Python script
        result = subprocess.run(
            [sys.executable, step_info['script']],
            capture_output=False,
            text=True,
            check=True
        )
        
        duration = time.time() - start_time
        print_step_complete(step_info, duration)
        
        return True, duration
        
    except subprocess.CalledProcessError as e:
        duration = time.time() - start_time
        print(f"\n✗ ERROR in Step {step_info['step']}: {step_info['name']}")
        print(f"Duration before error: {duration:.1f} seconds")
        print(f"Error details: {e}")
        return False, duration
    
    except FileNotFoundError:
        print(f"\n✗ ERROR: Script not found: {step_info['script']}")
        print(f"Make sure all step scripts are in the same directory as main.py")
        return False, 0
    
    except Exception as e:
        duration = time.time() - start_time
        print(f"\n✗ UNEXPECTED ERROR in Step {step_info['step']}")
        print(f"Error: {e}")
        return False, duration


def print_final_summary(total_duration, completed_steps):
    """Print final summary"""
    print("\n" + "="*80)
    print(" PIPELINE EXECUTION COMPLETE ".center(80))
    print("="*80)
    print(f"\nTotal Duration: {total_duration/60:.1f} minutes")
    print(f"Completed Steps: {completed_steps}/{len(PIPELINE_STEPS)}")
    print(f"Genes Processed: {len(GENE_LIST)}")
    print(f"\nOutput Files Generated:")
    
    for i in range(1, completed_steps + 1):
        output_file = f"Step{i}.csv"
        if Path(output_file).exists():
            file_size = Path(output_file).stat().st_size / 1024  # KB
            print(f"  ✓ {output_file} ({file_size:.1f} KB)")
    
    print("\n" + "="*80)
    print(" RESULTS READY FOR ANALYSIS ".center(80))
    print("="*80 + "\n")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main pipeline orchestrator"""
    
    print_banner()
    
    # Validate gene list
    if not GENE_LIST or len(GENE_LIST) == 0:
        print("✗ ERROR: GENE_LIST is empty!")
        print("Please add genes to GENE_LIST at the top of main.py")
        return
    
    print(f"📋 Gene List ({len(GENE_LIST)} genes):")
    for i, gene in enumerate(GENE_LIST, 1):
        print(f"  {i:3d}. {gene}")
    print()
    
    # Save gene list for step01
    save_gene_list(GENE_LIST)
    
    # Confirm execution
    print("\n" + "="*80)
    response = input("Proceed with pipeline execution? (yes/no): ").strip().lower()
    if response not in ['yes', 'y']:
        print("Pipeline execution cancelled.")
        return
    
    print("\n🚀 Starting pipeline execution...\n")
    
    # Track progress
    pipeline_start_time = time.time()
    step_durations = []
    completed_steps = 0
    
    # Execute each step
    for step_info in PIPELINE_STEPS:
        success, duration = run_step(step_info)
        step_durations.append(duration)
        
        if success:
            completed_steps += 1
        else:
            print(f"\n✗ Pipeline stopped at Step {step_info['step']}")
            print(f"Please fix errors and re-run from this step")
            break
        
        # Brief pause between steps
        time.sleep(1)
    
    # Calculate total duration
    total_duration = time.time() - pipeline_start_time
    
    # Print summary
    print_final_summary(total_duration, completed_steps)
    
    # Print step timing breakdown
    if completed_steps > 0:
        print("\n📊 Step Timing Breakdown:")
        print("─"*80)
        for i, (step_info, duration) in enumerate(zip(PIPELINE_STEPS[:completed_steps], step_durations)):
            pct = (duration / total_duration) * 100
            print(f"  Step {step_info['step']:2d}: {duration:6.1f}s ({pct:5.1f}%) - {step_info['name']}")
        print("─"*80)
    
    # Success message
    if completed_steps == len(PIPELINE_STEPS):
        print("\n✅ ALL STEPS COMPLETED SUCCESSFULLY!")
        print(f"\n📁 Final results: Step12.csv")
        print(f"📁 Final predictions: Step12_final_predictions.csv")
        print(f"📁 Performance metrics: Step12_performance_metrics.txt")
    else:
        print(f"\n⚠️  Pipeline incomplete: {completed_steps}/{len(PIPELINE_STEPS)} steps completed")


if __name__ == "__main__":
    main()
