"""
Step 12: Validate Top Candidates & Clinical Intelligence with Performance Metrics
Enhanced with AUC, Accuracy, and Correlation analysis between predicted scores and literature evidence
"""

import pandas as pd
import numpy as np
import requests
import time
from typing import Dict, List, Optional, Set, Tuple
from pathlib import Path
import logging
from datetime import datetime
from scipy import stats
from sklearn.metrics import roc_auc_score, accuracy_score, precision_recall_fscore_support, confusion_matrix
import warnings
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('step12_validation.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class CandidateValidator:
    """Validate top candidates with clinical trials, literature, and performance metrics"""
    
    def __init__(self):
        self.pubmed_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.clinicaltrials_base_url = "https://clinicaltrials.gov/api/v2"
        self.cache = {}
        
    def load_step11_data(self, filepath: str = 'Step11.csv') -> pd.DataFrame:
        """Load Step 11 data"""
        logger.info(f"Loading Step 11 data from {filepath}")
        df = pd.read_csv(filepath)
        logger.info(f"Successfully loaded {len(df)} entries")
        return df
    
    def _find_gene_column(self, df: pd.DataFrame) -> str:
        """Find gene identifier column"""
        for col in ['gene_symbol', 'Gene', 'gene', 'Gene_Symbol']:
            if col in df.columns:
                return col
        raise ValueError(f"Could not find gene column")
    
    def search_pubmed(self, gene_symbol: str, disease_term: str = None) -> Dict:
        """
        Search PubMed for gene publications (disease-agnostic)
        """
        
        logger.info(f"  Searching PubMed for {gene_symbol}...")
        
        try:
            # Disease-agnostic query - just search for the gene
            if disease_term:
                query = f"({gene_symbol}[Title/Abstract]) AND ({disease_term}[Title/Abstract])"
            else:
                # General search for gene
                query = f"{gene_symbol}[Gene Name] OR {gene_symbol}[Title/Abstract]"
            
            # Get count
            search_url = f"{self.pubmed_base_url}/esearch.fcgi"
            params = {
                'db': 'pubmed',
                'term': query,
                'retmax': 10,
                'sort': 'relevance',
                'retmode': 'json'
            }
            
            response = requests.get(search_url, params=params, timeout=15)
            
            if response.status_code == 200:
                data = response.json()
                
                if 'esearchresult' in data:
                    result = data['esearchresult']
                    count = int(result.get('count', 0))
                    pmids = result.get('idlist', [])
                    
                    logger.info(f"    Found {count} publications")
                    
                    # Get recent paper details
                    recent_papers = []
                    if pmids and len(pmids) > 0:
                        time.sleep(0.5)
                        
                        summary_url = f"{self.pubmed_base_url}/esummary.fcgi"
                        summary_params = {
                            'db': 'pubmed',
                            'id': ','.join(pmids[:5]),
                            'retmode': 'json'
                        }
                        
                        summary_response = requests.get(summary_url, params=summary_params, timeout=15)
                        
                        if summary_response.status_code == 200:
                            summary_data = summary_response.json()
                            
                            if 'result' in summary_data:
                                for pmid in pmids[:5]:
                                    if pmid in summary_data['result']:
                                        paper = summary_data['result'][pmid]
                                        recent_papers.append({
                                            'pmid': pmid,
                                            'title': paper.get('title', 'N/A'),
                                            'year': paper.get('pubdate', 'N/A').split()[0]
                                        })
                    
                    return {
                        'total_publications': count,
                        'recent_papers': recent_papers,
                        'has_evidence': count > 0
                    }
            
            logger.warning(f"    PubMed search failed")
            return {
                'total_publications': 0,
                'recent_papers': [],
                'has_evidence': False
            }
            
        except Exception as e:
            logger.warning(f"    Error: {e}")
            return {
                'total_publications': 0,
                'recent_papers': [],
                'has_evidence': False
            }
    
    def search_clinical_trials(self, gene_symbol: str, protein_name: str = None) -> Dict:
        """
        Search ClinicalTrials.gov for ongoing trials
        """
        
        logger.info(f"  Searching ClinicalTrials.gov for {gene_symbol}...")
        
        try:
            # Build search query
            search_terms = [gene_symbol]
            if protein_name and protein_name not in ['Unknown', 'N/A']:
                search_terms.append(protein_name.split()[0])
            
            query = ' OR '.join(search_terms)
            
            url = f"{self.clinicaltrials_base_url}/studies"
            params = {
                'query.term': query,
                'pageSize': 10,
                'format': 'json'
            }
            
            response = requests.get(url, params=params, timeout=15)
            
            if response.status_code == 200:
                data = response.json()
                
                studies = data.get('studies', [])
                total_count = len(studies)
                
                if total_count > 0:
                    logger.info(f"    Found {total_count} clinical trials")
                    
                    trials = []
                    active_trials = 0
                    recruiting = 0
                    
                    for study in studies[:5]:
                        protocol = study.get('protocolSection', {})
                        identification = protocol.get('identificationModule', {})
                        status = protocol.get('statusModule', {})
                        
                        trial_info = {
                            'nct_id': identification.get('nctId', 'N/A'),
                            'title': identification.get('briefTitle', 'N/A'),
                            'status': status.get('overallStatus', 'Unknown'),
                            'phase': status.get('phase', 'N/A')
                        }
                        
                        trials.append(trial_info)
                        
                        if 'ACTIVE' in trial_info['status'].upper() or 'RECRUITING' in trial_info['status'].upper():
                            active_trials += 1
                        
                        if 'RECRUITING' in trial_info['status'].upper():
                            recruiting += 1
                    
                    return {
                        'total_trials': total_count,
                        'active_trials': active_trials,
                        'recruiting_trials': recruiting,
                        'trials': trials,
                        'has_trials': True
                    }
                else:
                    logger.info(f"    No trials found")
                    return {
                        'total_trials': 0,
                        'active_trials': 0,
                        'recruiting_trials': 0,
                        'trials': [],
                        'has_trials': False
                    }
            else:
                logger.warning(f"    ClinicalTrials API returned {response.status_code}")
                return {
                    'total_trials': 0,
                    'active_trials': 0,
                    'recruiting_trials': 0,
                    'trials': [],
                    'has_trials': False
                }
                
        except Exception as e:
            logger.warning(f"    Error: {e}")
            return {
                'total_trials': 0,
                'active_trials': 0,
                'recruiting_trials': 0,
                'trials': [],
                'has_trials': False
            }
    
    def identify_tool_compounds(self, row: pd.Series) -> Dict:
        """Identify available tool compounds and chemical probes"""
        
        tool_info = {
            'has_tool_compounds': False,
            'approved_drugs': 0,
            'bioactive_compounds': 0,
            'chembl_compounds': 0,
            'tool_compound_summary': '',
            'probe_quality': 'Unknown'
        }
        
        # Check approved drugs
        has_drugs = row.get('has_approved_drugs', False)
        dgidb_approved = row.get('dgidb_num_approved', 0)
        if pd.isna(dgidb_approved):
            dgidb_approved = 0
        
        tool_info['approved_drugs'] = int(dgidb_approved)
        
        # Check bioactive compounds
        has_bioactive = row.get('has_bioactive_compounds', False)
        chembl_bioactive = row.get('chembl_bioactive_compounds', 0)
        if pd.isna(chembl_bioactive):
            chembl_bioactive = 0
        
        tool_info['bioactive_compounds'] = int(chembl_bioactive)
        
        # ChEMBL compounds
        chembl_total = row.get('chembl_total_compounds', 0)
        if pd.isna(chembl_total):
            chembl_total = 0
        
        tool_info['chembl_compounds'] = int(chembl_total)
        
        # Determine if tool compounds exist
        if tool_info['approved_drugs'] > 0:
            tool_info['has_tool_compounds'] = True
            tool_info['probe_quality'] = 'Excellent - Approved drugs available'
            tool_info['tool_compound_summary'] = f"{tool_info['approved_drugs']} approved drug(s) as chemical probes"
        elif tool_info['bioactive_compounds'] >= 10:
            tool_info['has_tool_compounds'] = True
            tool_info['probe_quality'] = 'Good - Multiple bioactive compounds'
            tool_info['tool_compound_summary'] = f"{tool_info['bioactive_compounds']} bioactive compounds available"
        elif tool_info['bioactive_compounds'] > 0:
            tool_info['has_tool_compounds'] = True
            tool_info['probe_quality'] = 'Moderate - Limited compounds'
            tool_info['tool_compound_summary'] = f"{tool_info['bioactive_compounds']} bioactive compound(s)"
        elif tool_info['chembl_compounds'] > 0:
            tool_info['has_tool_compounds'] = False
            tool_info['probe_quality'] = 'Poor - Only screening compounds'
            tool_info['tool_compound_summary'] = f"{tool_info['chembl_compounds']} screening compounds (no validated probes)"
        else:
            tool_info['probe_quality'] = 'None - No compounds available'
            tool_info['tool_compound_summary'] = 'No chemical tools available - de novo discovery required'
        
        return tool_info
    
    def identify_combination_approaches(self, row: pd.Series) -> Dict:
        """Identify potential combination therapy approaches"""
        
        combination_info = {
            'has_combination_potential': False,
            'connected_drug_targets': 0,
            'target_families': [],
            'combination_strategies': [],
            'synthetic_lethality_candidates': [],
            'pathway_based_combinations': []
        }
        
        # Network-based combinations
        drug_target_interactions = row.get('Step10_Drug_Target_Interactions', 0)
        if not pd.isna(drug_target_interactions) and drug_target_interactions > 0:
            combination_info['has_combination_potential'] = True
            combination_info['connected_drug_targets'] = int(drug_target_interactions)
            
            partners = row.get('Step10_Drug_Target_Partners', '')
            if partners and partners != '':
                partner_list = [p.strip() for p in str(partners).split(';')]
                combination_info['combination_strategies'].append(
                    f"Polypharmacology: Target alongside {', '.join(partner_list[:3])}"
                )
            
            families = row.get('Step10_Drug_Target_Families', '')
            if families and families != '':
                family_list = [f.strip() for f in str(families).split(';')]
                combination_info['target_families'] = family_list
                combination_info['combination_strategies'].append(
                    f"Multi-target approach within {', '.join(family_list[:2])} pathways"
                )
        
        # Pathway-based combinations
        pathways = row.get('Step10_Pathways', '')
        if pathways and pathways != '':
            pathway_list = [p.strip() for p in str(pathways).split(';')]
            
            for pathway in pathway_list[:3]:
                if 'inflammation' in pathway.lower():
                    combination_info['pathway_based_combinations'].append(
                        "Anti-inflammatory combination (e.g., with NSAIDs or biologics)"
                    )
                elif 'serotonin' in pathway.lower():
                    combination_info['pathway_based_combinations'].append(
                        "Serotonergic modulation combination (e.g., with triptans or SSRIs)"
                    )
                elif 'calcium' in pathway.lower() or 'ion' in pathway.lower():
                    combination_info['pathway_based_combinations'].append(
                        "Ion channel modulation (e.g., with calcium channel blockers)"
                    )
        
        # Synthetic lethality
        family = row.get('druggable_family_type', 'Unknown')
        if pd.isna(family):
            family = 'Unknown'
        else:
            family = str(family)
        
        if 'GPCR' in family:
            combination_info['synthetic_lethality_candidates'].append(
                "GPCR desensitization pathway components (GRKs, arrestins)"
            )
        elif 'Kinase' in family:
            combination_info['synthetic_lethality_candidates'].append(
                "Synthetic lethality with phosphatase inhibitors"
            )
        elif 'Enzyme' in family:
            combination_info['synthetic_lethality_candidates'].append(
                "Metabolic compensation pathway targets"
            )
        
        if combination_info['has_combination_potential']:
            combination_info['combination_strategies'].append(
                "Sequential dosing to minimize resistance"
            )
        
        return combination_info
    
    def generate_literature_verification_protocol(self, row: pd.Series, gene_symbol: str) -> Dict:
        """
        Generate detailed literature verification protocol
        """
        
        # Get target characteristics
        best_modality = row.get('Step11_Best_Modality', 'Unknown')
        sm_score = row.get('Step11_SmallMolecule_Score', 0)
        bio_score = row.get('Step11_Biologic_Score', 0)
        is_biologic_only = row.get('Step11_Is_Biologic_Only', False)
        
        family = str(row.get('Step9_Protein_Family', 'Unknown'))
        has_pocket = row.get('has_druggable_pocket', False)
        has_drugs = row.get('has_approved_drugs', False)
        location = str(row.get('subcellular_location', 'Unknown'))
        
        verification_protocol = {
            'gene': gene_symbol,
            'predicted_modality': best_modality,
            'literature_search_queries': [],
            'key_evidence_to_find': [],
            'validation_experiments': [],
            'red_flags_to_check': [],
            'expected_druggability_markers': [],
            'literature_validation_checklist': []
        }
        
        # Basic queries
        verification_protocol['literature_search_queries'].append(
            f'"{gene_symbol}" AND (drug OR inhibitor OR agonist OR antagonist)'
        )
        
        if best_modality == 'Small Molecule' or sm_score > bio_score:
            verification_protocol['literature_search_queries'].extend([
                f'"{gene_symbol}" AND ("small molecule" OR "chemical probe" OR "tool compound")',
                f'"{gene_symbol}" AND (IC50 OR Ki OR Kd OR "binding affinity")',
                f'"{gene_symbol}" AND ("structure activity relationship" OR SAR OR "medicinal chemistry")',
                f'"{gene_symbol}" AND ("crystal structure" OR "binding site" OR "active site")'
            ])
        
        if best_modality == 'Biologic' or bio_score > sm_score or is_biologic_only:
            verification_protocol['literature_search_queries'].extend([
                f'"{gene_symbol}" AND (antibody OR biologics OR "monoclonal antibody")',
                f'"{gene_symbol}" AND (secreted OR extracellular OR membrane)',
                f'"{gene_symbol}" AND (neutralizing OR blocking OR "antibody therapy")'
            ])
        
        # Family-specific queries
        if family == 'GPCR':
            verification_protocol['literature_search_queries'].extend([
                f'"{gene_symbol}" AND (radioligand OR "binding assay" OR "functional assay")',
                f'"{gene_symbol}" AND ("G protein" OR "beta-arrestin" OR signaling)'
            ])
        elif family == 'Kinase':
            verification_protocol['literature_search_queries'].extend([
                f'"{gene_symbol}" AND (kinase AND (inhibitor OR selectivity))',
                f'"{gene_symbol}" AND ("kinome" OR "kinase panel" OR "off-target")'
            ])
        elif family == 'IC':
            verification_protocol['literature_search_queries'].append(
                f'"{gene_symbol}" AND ("ion channel" AND (blocker OR modulator))'
            )
        
        # Key evidence
        if has_drugs:
            verification_protocol['key_evidence_to_find'].extend([
                '1. Approved drug mechanisms',
                '2. Clinical efficacy data',
                '3. Drug binding sites'
            ])
        else:
            verification_protocol['key_evidence_to_find'].extend([
                '1. Tool compound papers',
                '2. Proof-of-concept studies',
                '3. Target validation'
            ])
        
        # Validation experiments
        if best_modality == 'Small Molecule':
            verification_protocol['validation_experiments'] = [
                'Biochemical assays', 'Cellular assays', 'Structural biology', 'In vivo validation'
            ]
        else:
            verification_protocol['validation_experiments'] = [
                'Antibody generation', 'Functional assays', 'Characterization', 'In vivo validation'
            ]
        
        # Red flags
        verification_protocol['red_flags_to_check'] = [
            'Target validity issues', 'Druggability concerns', 'Safety concerns'
        ]
        
        # Expected markers
        if best_modality == 'Small Molecule':
            verification_protocol['expected_druggability_markers'] = [
                'Tool compounds with IC50 < 1 micro M',
                'Co-crystal structures',
                'Multiple chemical scaffolds',
                'Patent literature'
            ]
        else:
            verification_protocol['expected_druggability_markers'] = [
                'Secreted or extracellular',
                'Neutralizing antibodies',
                'Epitope mapping',
                'Antibody therapies'
            ]
        
        # Checklist
        verification_protocol['literature_validation_checklist'] = [
            'Tool compounds found',
            'Binding site characterized',
            'Target validated',
            'Clinical precedent',
            'Safety assessed'
        ]
        
        return verification_protocol
    
    def calculate_literature_druggability_ground_truth(self, row: pd.Series) -> Dict:
        """
        Calculate ground truth druggability from literature evidence
        This serves as the "true label" for performance metrics
        """
        
        # Get literature evidence
        has_approved = row.get('has_approved_drugs', False)
        num_approved = row.get('dgidb_num_approved', 0)
        if pd.isna(num_approved):
            num_approved = 0
        
        has_clinical = row.get('has_clinical_drugs', False)
        has_bioactive = row.get('has_bioactive_compounds', False)
        num_bioactive = row.get('chembl_bioactive_compounds', 0)
        if pd.isna(num_bioactive):
            num_bioactive = 0
        
        # Get drug types
        drug_types = str(row.get('opentargets_drug_types', ''))
        
        # Get TDL
        tdl = str(row.get('Step9_Target_Development_Level', 'Unknown'))
        
        # Ground truth for Small Molecule
        sm_ground_truth = 0
        
        if 'Small molecule' in drug_types and num_approved >= 10:
            sm_ground_truth = 100  # Highly druggable
        elif 'Small molecule' in drug_types and num_approved >= 5:
            sm_ground_truth = 90
        elif 'Small molecule' in drug_types and num_approved >= 3:
            sm_ground_truth = 80
        elif 'Small molecule' in drug_types and num_approved >= 1:
            sm_ground_truth = 70
        elif has_clinical and 'Small molecule' in drug_types:
            sm_ground_truth = 60
        elif num_bioactive >= 100:
            sm_ground_truth = 50
        elif num_bioactive >= 50:
            sm_ground_truth = 40
        elif num_bioactive >= 10:
            sm_ground_truth = 30
        elif tdl == 'Tclin':
            sm_ground_truth = 50
        elif tdl == 'Tchem':
            sm_ground_truth = 40
        else:
            sm_ground_truth = 20  # Low/unknown
        
        # Ground truth for Biologic
        bio_ground_truth = 0
        
        if any(x in drug_types for x in ['Antibody', 'Protein']) and num_approved >= 5:
            bio_ground_truth = 100
        elif any(x in drug_types for x in ['Antibody', 'Protein']) and num_approved >= 1:
            bio_ground_truth = 90
        elif has_clinical and any(x in drug_types for x in ['Antibody', 'Protein']):
            bio_ground_truth = 70
        elif 'Secreted' in str(row.get('subcellular_location', '')):
            bio_ground_truth = 60
        elif 'Cell membrane' in str(row.get('subcellular_location', '')):
            bio_ground_truth = 50
        else:
            bio_ground_truth = 30
        
        # Binary labels (for classification metrics)
        sm_binary = 1 if sm_ground_truth >= 70 else 0
        bio_binary = 1 if bio_ground_truth >= 70 else 0
        
        return {
            'SM_Literature_Score': sm_ground_truth,
            'Bio_Literature_Score': bio_ground_truth,
            'SM_Is_Druggable': sm_binary,
            'Bio_Is_Druggable': bio_binary
        }
    
    def validate_candidate(self, row: pd.Series, gene_symbol: str, disease: str = "migraine") -> Dict:
        """Complete validation for a high-priority candidate"""
        
        cache_key = f"{gene_symbol}_{disease}"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        logger.info(f"[Validating] {gene_symbol}")
        
        validation_data = {
            'gene': gene_symbol,
            'composite_score': row.get('Step11_Composite_Score', 0),
            'priority': row.get('Step11_Priority', 'Unknown')
        }
        
        # Literature evidence
        protein_name = row.get('uniprot_protein_name', None)
        pubmed_data = self.search_pubmed(gene_symbol, disease)
        validation_data.update({
            'literature_publications': pubmed_data['total_publications'],
            'has_literature_evidence': pubmed_data['has_evidence'],
            'recent_papers': pubmed_data['recent_papers']
        })
        
        time.sleep(0.5)
        
        # Clinical trials
        trials_data = self.search_clinical_trials(gene_symbol, protein_name)
        validation_data.update({
            'total_trials': trials_data['total_trials'],
            'active_trials': trials_data['active_trials'],
            'recruiting_trials': trials_data['recruiting_trials'],
            'has_clinical_trials': trials_data['has_trials'],
            'trial_details': trials_data['trials']
        })
        
        time.sleep(0.5)
        
        # Tool compounds
        tool_data = self.identify_tool_compounds(row)
        validation_data.update({
            'has_tool_compounds': tool_data['has_tool_compounds'],
            'tool_compound_count': tool_data['approved_drugs'] + tool_data['bioactive_compounds'],
            'probe_quality': tool_data['probe_quality'],
            'tool_compound_summary': tool_data['tool_compound_summary']
        })
        
        # Combination approaches
        combo_data = self.identify_combination_approaches(row)
        validation_data.update({
            'has_combination_potential': combo_data['has_combination_potential'],
            'combination_strategies': combo_data['combination_strategies'],
            'synthetic_lethality_candidates': combo_data['synthetic_lethality_candidates'],
            'pathway_combinations': combo_data['pathway_based_combinations']
        })
        
        # Literature verification protocol
        verification_protocol = self.generate_literature_verification_protocol(row, gene_symbol)
        validation_data['verification_protocol'] = verification_protocol
        
        # Ground truth from literature
        ground_truth = self.calculate_literature_druggability_ground_truth(row)
        validation_data.update(ground_truth)
        
        # Overall validation score
        validation_score = 0
        if validation_data['has_literature_evidence']:
            validation_score += 3
        if validation_data['has_clinical_trials']:
            validation_score += 3
        if validation_data['has_tool_compounds']:
            validation_score += 2
        if validation_data['has_combination_potential']:
            validation_score += 2
        
        validation_data['validation_score'] = validation_score
        validation_data['validation_level'] = self._classify_validation(validation_score)
        
        self.cache[cache_key] = validation_data
        
        logger.info(f"  Validation: {validation_data['validation_level']} " +
                   f"(Lit: {validation_data['literature_publications']}, " +
                   f"Trials: {validation_data['total_trials']})")
        
        return validation_data
    
    def _classify_validation(self, score: int) -> str:
        """Classify validation level based on score"""
        if score >= 8:
            return 'Excellent Validation'
        elif score >= 6:
            return 'Strong Validation'
        elif score >= 4:
            return 'Moderate Validation'
        elif score >= 2:
            return 'Preliminary Validation'
        else:
            return 'Limited Validation'
    
    def calculate_performance_metrics(self, df: pd.DataFrame) -> Dict:
        """
        Calculate AUC, Accuracy, and Correlation for Small Molecule and Biologic predictions
        """
        
        logger.info("\n" + "="*80)
        logger.info("CALCULATING PERFORMANCE METRICS")
        logger.info("="*80)
        
        # Filter for genes with ground truth
        valid_genes = df[
            (df['SM_Literature_Score'].notna()) &
            (df['Bio_Literature_Score'].notna()) &
            (df['Step11_SmallMolecule_Score'].notna()) &
            (df['Step11_Biologic_Score'].notna())
        ].copy()
        
        if len(valid_genes) == 0:
            logger.warning("No genes with both predictions and ground truth")
            return {}
        
        logger.info(f"Analyzing {len(valid_genes)} genes with complete data")
        
        metrics = {}
        
        # ==================== SMALL MOLECULE METRICS ====================
        
        logger.info("\n--- Small Molecule Druggability ---")
        
        sm_predicted = valid_genes['Step11_SmallMolecule_Score'].values
        sm_truth_continuous = valid_genes['SM_Literature_Score'].values
        sm_truth_binary = valid_genes['SM_Is_Druggable'].values
        
        # Normalize predicted scores to 0-100 range (they should already be)
        sm_predicted_normalized = np.clip(sm_predicted, 0, 100)
        
        # Binarize predictions (threshold at 70)
        sm_predicted_binary = (sm_predicted_normalized >= 70).astype(int)
        
        # Calculate metrics
        try:
            # AUC (using continuous scores)
            sm_auc = roc_auc_score(sm_truth_binary, sm_predicted_normalized / 100.0)
            logger.info(f"  Small Molecule AUC: {sm_auc:.3f}")
            
            # Accuracy
            sm_accuracy = accuracy_score(sm_truth_binary, sm_predicted_binary)
            logger.info(f"  Small Molecule Accuracy: {sm_accuracy:.3f} ({sm_accuracy*100:.1f}%)")
            
            # Precision, Recall, F1
            sm_precision, sm_recall, sm_f1, _ = precision_recall_fscore_support(
                sm_truth_binary, sm_predicted_binary, average='binary', zero_division=0
            )
            logger.info(f"  Small Molecule Precision: {sm_precision:.3f}")
            logger.info(f"  Small Molecule Recall: {sm_recall:.3f}")
            logger.info(f"  Small Molecule F1-Score: {sm_f1:.3f}")
            
            # Confusion Matrix
            sm_cm = confusion_matrix(sm_truth_binary, sm_predicted_binary)
            logger.info(f"  Confusion Matrix:")
            logger.info(f"    TN={sm_cm[0,0]}, FP={sm_cm[0,1]}")
            logger.info(f"    FN={sm_cm[1,0]}, TP={sm_cm[1,1]}")
            
            # Correlation (continuous scores)
            sm_pearson_r, sm_pearson_p = stats.pearsonr(sm_predicted_normalized, sm_truth_continuous)
            sm_spearman_r, sm_spearman_p = stats.spearmanr(sm_predicted_normalized, sm_truth_continuous)
            
            logger.info(f"  Pearson Correlation: r={sm_pearson_r:.3f}, p={sm_pearson_p:.4f}")
            logger.info(f"  Spearman Correlation: ?={sm_spearman_r:.3f}, p={sm_spearman_p:.4f}")
            
            metrics['Small_Molecule'] = {
                'AUC': sm_auc,
                'Accuracy': sm_accuracy,
                'Precision': sm_precision,
                'Recall': sm_recall,
                'F1_Score': sm_f1,
                'Pearson_r': sm_pearson_r,
                'Pearson_p': sm_pearson_p,
                'Spearman_r': sm_spearman_r,
                'Spearman_p': sm_spearman_p,
                'Confusion_Matrix': sm_cm.tolist(),
                'N_Genes': len(valid_genes)
            }
            
        except Exception as e:
            logger.error(f"Error calculating SM metrics: {e}")
            metrics['Small_Molecule'] = {'Error': str(e)}
        
        # ==================== BIOLOGIC METRICS ====================
        
        logger.info("\n--- Biologic Druggability ---")
        
        bio_predicted = valid_genes['Step11_Biologic_Score'].values
        bio_truth_continuous = valid_genes['Bio_Literature_Score'].values
        bio_truth_binary = valid_genes['Bio_Is_Druggable'].values
        
        # Normalize predicted scores
        bio_predicted_normalized = np.clip(bio_predicted, 0, 100)
        
        # Binarize predictions (threshold at 70)
        bio_predicted_binary = (bio_predicted_normalized >= 70).astype(int)
        
        try:
            # AUC
            bio_auc = roc_auc_score(bio_truth_binary, bio_predicted_normalized / 100.0)
            logger.info(f"  Biologic AUC: {bio_auc:.3f}")
            
            # Accuracy
            bio_accuracy = accuracy_score(bio_truth_binary, bio_predicted_binary)
            logger.info(f"  Biologic Accuracy: {bio_accuracy:.3f} ({bio_accuracy*100:.1f}%)")
            
            # Precision, Recall, F1
            bio_precision, bio_recall, bio_f1, _ = precision_recall_fscore_support(
                bio_truth_binary, bio_predicted_binary, average='binary', zero_division=0
            )
            logger.info(f"  Biologic Precision: {bio_precision:.3f}")
            logger.info(f"  Biologic Recall: {bio_recall:.3f}")
            logger.info(f"  Biologic F1-Score: {bio_f1:.3f}")
            
            # Confusion Matrix
            bio_cm = confusion_matrix(bio_truth_binary, bio_predicted_binary)
            logger.info(f"  Confusion Matrix:")
            logger.info(f"    TN={bio_cm[0,0]}, FP={bio_cm[0,1]}")
            logger.info(f"    FN={bio_cm[1,0]}, TP={bio_cm[1,1]}")
            
            # Correlation
            bio_pearson_r, bio_pearson_p = stats.pearsonr(bio_predicted_normalized, bio_truth_continuous)
            bio_spearman_r, bio_spearman_p = stats.spearmanr(bio_predicted_normalized, bio_truth_continuous)
            
            logger.info(f"  Pearson Correlation: r={bio_pearson_r:.3f}, p={bio_pearson_p:.4f}")
            logger.info(f"  Spearman Correlation: ?={bio_spearman_r:.3f}, p={bio_spearman_p:.4f}")
            
            metrics['Biologic'] = {
                'AUC': bio_auc,
                'Accuracy': bio_accuracy,
                'Precision': bio_precision,
                'Recall': bio_recall,
                'F1_Score': bio_f1,
                'Pearson_r': bio_pearson_r,
                'Pearson_p': bio_pearson_p,
                'Spearman_r': bio_spearman_r,
                'Spearman_p': bio_spearman_p,
                'Confusion_Matrix': bio_cm.tolist(),
                'N_Genes': len(valid_genes)
            }
            
        except Exception as e:
            logger.error(f"Error calculating Biologic metrics: {e}")
            metrics['Biologic'] = {'Error': str(e)}
        
        return metrics
    
    def generate_metrics_report(self, metrics: Dict, output_file: str = 'Step12_performance_metrics.txt'):
        """Generate detailed metrics report"""
        
        if not metrics:
            logger.warning("No metrics to report")
            return
        
        with open(output_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("DRUGGABILITY PREDICTION PERFORMANCE METRICS\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("="*80 + "\n\n")
            
            f.write("This report compares pipeline predictions against literature-derived ground truth\n")
            f.write("to assess the accuracy of small molecule and biologic druggability scoring.\n\n")
            
            # Small Molecule Metrics
            if 'Small_Molecule' in metrics and 'Error' not in metrics['Small_Molecule']:
                sm = metrics['Small_Molecule']
                
                f.write("="*80 + "\n")
                f.write("SMALL MOLECULE DRUGGABILITY METRICS\n")
                f.write("="*80 + "\n\n")
                
                f.write(f"Number of Genes Analyzed: {sm['N_Genes']}\n\n")
                
                f.write("CLASSIFICATION METRICS:\n")
                f.write(f"  AUC (Area Under ROC Curve):    {sm['AUC']:.3f}\n")
                f.write(f"  Accuracy:                       {sm['Accuracy']:.3f} ({sm['Accuracy']*100:.1f}%)\n")
                f.write(f"  Precision:                      {sm['Precision']:.3f}\n")
                f.write(f"  Recall (Sensitivity):           {sm['Recall']:.3f}\n")
                f.write(f"  F1-Score:                       {sm['F1_Score']:.3f}\n\n")
                
                cm = sm['Confusion_Matrix']
                f.write("CONFUSION MATRIX:\n")
                f.write(f"                    Predicted Non-Druggable    Predicted Druggable\n")
                f.write(f"  Actual Non-Drug   {cm[0][0]:^23} {cm[0][1]:^23}\n")
                f.write(f"  Actual Druggable  {cm[1][0]:^23} {cm[1][1]:^23}\n\n")
                
                f.write("CORRELATION METRICS:\n")
                f.write(f"  Pearson Correlation:   r = {sm['Pearson_r']:>6.3f}  (p = {sm['Pearson_p']:.4f})\n")
                f.write(f"  Spearman Correlation:  ? = {sm['Spearman_r']:>6.3f}  (p = {sm['Spearman_p']:.4f})\n\n")
                
                f.write("INTERPRETATION:\n")
                if sm['AUC'] >= 0.9:
                    f.write("  AUC = 0.90: Excellent discrimination\n")
                elif sm['AUC'] >= 0.8:
                    f.write("  AUC = 0.80: Good discrimination\n")
                elif sm['AUC'] >= 0.7:
                    f.write("  AUC = 0.70: Acceptable discrimination\n")
                else:
                    f.write("  AUC < 0.70: Poor discrimination\n")
                
                if abs(sm['Pearson_r']) >= 0.7:
                    f.write("  |r| = 0.70: Strong correlation\n")
                elif abs(sm['Pearson_r']) >= 0.5:
                    f.write("  |r| = 0.50: Moderate correlation\n")
                else:
                    f.write("  |r| < 0.50: Weak correlation\n")
                
                f.write("\n")
            
            # Biologic Metrics
            if 'Biologic' in metrics and 'Error' not in metrics['Biologic']:
                bio = metrics['Biologic']
                
                f.write("="*80 + "\n")
                f.write("BIOLOGIC DRUGGABILITY METRICS\n")
                f.write("="*80 + "\n\n")
                
                f.write(f"Number of Genes Analyzed: {bio['N_Genes']}\n\n")
                
                f.write("CLASSIFICATION METRICS:\n")
                f.write(f"  AUC (Area Under ROC Curve):    {bio['AUC']:.3f}\n")
                f.write(f"  Accuracy:                       {bio['Accuracy']:.3f} ({bio['Accuracy']*100:.1f}%)\n")
                f.write(f"  Precision:                      {bio['Precision']:.3f}\n")
                f.write(f"  Recall (Sensitivity):           {bio['Recall']:.3f}\n")
                f.write(f"  F1-Score:                       {bio['F1_Score']:.3f}\n\n")
                
                cm = bio['Confusion_Matrix']
                f.write("CONFUSION MATRIX:\n")
                f.write(f"                    Predicted Non-Druggable    Predicted Druggable\n")
                f.write(f"  Actual Non-Drug   {cm[0][0]:^23} {cm[0][1]:^23}\n")
                f.write(f"  Actual Druggable  {cm[1][0]:^23} {cm[1][1]:^23}\n\n")
                
                f.write("CORRELATION METRICS:\n")
                f.write(f"  Pearson Correlation:   r = {bio['Pearson_r']:>6.3f}  (p = {bio['Pearson_p']:.4f})\n")
                f.write(f"  Spearman Correlation:  ? = {bio['Spearman_r']:>6.3f}  (p = {bio['Spearman_p']:.4f})\n\n")
                
                f.write("INTERPRETATION:\n")
                if bio['AUC'] >= 0.9:
                    f.write("  AUC = 0.90: Excellent discrimination\n")
                elif bio['AUC'] >= 0.8:
                    f.write("  AUC = 0.80: Good discrimination\n")
                elif bio['AUC'] >= 0.7:
                    f.write("  AUC = 0.70: Acceptable discrimination\n")
                else:
                    f.write("  AUC < 0.70: Poor discrimination\n")
                
                if abs(bio['Pearson_r']) >= 0.7:
                    f.write("  |r| = 0.70: Strong correlation\n")
                elif abs(bio['Pearson_r']) >= 0.5:
                    f.write("  |r| = 0.50: Moderate correlation\n")
                else:
                    f.write("  |r| < 0.50: Weak correlation\n")
                
                f.write("\n")
            
            f.write("="*80 + "\n")
            f.write("NOTES\n")
            f.write("="*80 + "\n\n")
            f.write("Ground Truth Definition:\n")
            f.write("  - Based on approved drugs, clinical trials, and bioactive compounds\n")
            f.write("  - Druggable threshold: =70 (on 0-100 scale)\n")
            f.write("  - Separate scoring for small molecule vs biologic modalities\n\n")
            f.write("Prediction Threshold:\n")
            f.write("  - Druggable: Score =70\n")
            f.write("  - Non-druggable: Score <70\n\n")
        
        logger.info(f"Performance metrics report saved: {output_file}")
    
    def generate_final_prediction(self, row: pd.Series, validation_data: Dict) -> Dict:
        """Generate final druggability prediction combining pipeline and literature"""
        
        gene = validation_data['gene']
        
        # Pipeline prediction
        composite_score = row.get('Step11_Composite_Score', 0)
        if pd.isna(composite_score):
            composite_score = 0
        
        sm_score = row.get('Step11_SmallMolecule_Score', 0)
        bio_score = row.get('Step11_Biologic_Score', 0)
        best_modality = row.get('Step11_Best_Modality', 'Unknown')
        
        classification = row.get('Step11_Classification', 'Unknown')
        
        # Literature validation
        has_lit = validation_data.get('has_literature_evidence', False)
        lit_count = validation_data.get('literature_publications', 0)
        has_trials = validation_data.get('has_clinical_trials', False)
        has_tools = validation_data.get('has_tool_compounds', False)
        
        # Ground truth scores
        sm_literature = validation_data.get('SM_Literature_Score', 0)
        bio_literature = validation_data.get('Bio_Literature_Score', 0)
        
        # Generate pipeline-based prediction
        if composite_score >= 80:
            pipeline_prediction = "Highly Druggable"
            pipeline_confidence = "High"
        elif composite_score >= 50:
            pipeline_prediction = "Moderately Druggable"
            pipeline_confidence = "Medium"
        else:
            pipeline_prediction = "Challenging Target"
            pipeline_confidence = "Low"
        
        # Pipeline reasoning
        reasons = []
        if row.get('has_approved_drugs', False):
            reasons.append("Approved drugs exist")
        if row.get('has_druggable_pocket', False):
            reasons.append("Druggable binding pocket")
        if row.get('is_druggable_family', False):
            reasons.append("Druggable protein family")
        if row.get('Step8_Passes_50_Percent_Rule', False):
            reasons.append("Low disorder (<50%)")
        
        tdl = row.get('Step9_Target_Development_Level', 'Unknown')
        if tdl == 'Tclin':
            reasons.append("Clinical target (Tclin)")
        elif tdl == 'Tchem':
            reasons.append("Chemical target (Tchem)")
        
        pipeline_reason = "; ".join(reasons) if reasons else "Limited evidence"
        
        # Literature-based validation
        if lit_count >= 10 and has_trials:
            lit_prediction = "Strong Clinical Evidence"
            lit_confidence = "High"
        elif lit_count >= 5 or has_trials:
            lit_prediction = "Moderate Clinical Evidence"
            lit_confidence = "Medium"
        elif lit_count > 0:
            lit_prediction = "Preliminary Evidence"
            lit_confidence = "Low"
        else:
            lit_prediction = "No Direct Evidence"
            lit_confidence = "Very Low"
        
        # Literature reasoning
        lit_reasons = []
        if lit_count > 0:
            lit_reasons.append(f"{lit_count} publications")
        if has_trials:
            trial_count = validation_data.get('total_trials', 0)
            lit_reasons.append(f"{trial_count} clinical trial(s)")
        if has_tools:
            lit_reasons.append("Chemical tools available")
        
        lit_reason = "; ".join(lit_reasons) if lit_reasons else "No literature evidence found"
        
        # Final integrated prediction
        if composite_score >= 80 and (lit_count >= 5 or has_trials):
            final_prediction = "PRIORITY TARGET - High druggability with clinical validation"
            final_confidence = "Very High"
            recommendation = "PROCEED: Initiate drug discovery program"
        elif composite_score >= 80:
            final_prediction = "HIGH POTENTIAL - Excellent druggability, needs validation"
            final_confidence = "High"
            recommendation = "VALIDATE: Confirm disease association before proceeding"
        elif composite_score >= 50 and lit_count >= 10:
            final_prediction = "VIABLE TARGET - Moderate druggability with strong evidence"
            final_confidence = "Medium-High"
            recommendation = "CONSIDER: Address druggability challenges"
        elif composite_score >= 50:
            final_prediction = "POSSIBLE TARGET - Requires optimization"
            final_confidence = "Medium"
            recommendation = "EXPLORE: Fragment-based or alternative approaches"
        else:
            final_prediction = "CHALLENGING - Alternative strategies needed"
            final_confidence = "Low"
            recommendation = "CAUTION: Consider biologics or indirect targeting"
        
        return {
            'Gene': gene,
            'Pipeline_Druggability_Score': round(composite_score, 1),
            'Small_Molecule_Score': round(sm_score, 1),
            'Biologic_Score': round(bio_score, 1),
            'SM_Literature_Score': sm_literature,
            'Bio_Literature_Score': bio_literature,
            'Best_Modality': best_modality,
            'Pipeline_Prediction': pipeline_prediction,
            'Pipeline_Confidence': pipeline_confidence,
            'Pipeline_Reasoning': pipeline_reason,
            'Literature_Publications': lit_count,
            'Clinical_Trials': validation_data.get('total_trials', 0),
            'Literature_Prediction': lit_prediction,
            'Literature_Confidence': lit_confidence,
            'Literature_Reasoning': lit_reason,
            'Final_Prediction': final_prediction,
            'Final_Confidence': final_confidence,
            'Recommendation': recommendation
        }
    
    def add_validation_data(self, step11_df: pd.DataFrame, disease: str = "migraine") -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Add validation data to all candidates (not just high priority)"""
        
        gene_col = self._find_gene_column(step11_df)
        
        # Process ALL genes for metrics calculation
        logger.info(f"\nValidating all {len(step11_df)} candidates for performance metrics...")
        logger.info("="*80)
        
        all_validation_data = []
        all_predictions = []
        
        for idx, (_, row) in enumerate(step11_df.iterrows(), 1):
            gene = row[gene_col]
            
            # Calculate ground truth for all genes
            ground_truth = self.calculate_literature_druggability_ground_truth(row)
            
            # For high priority, do full validation
            if row.get('Step11_Priority') == 'High Priority':
                validation_data = self.validate_candidate(row, gene, disease)
            else:
                # For others, just add ground truth
                validation_data = {
                    'gene': gene,
                    'composite_score': row.get('Step11_Composite_Score', 0),
                    'priority': row.get('Step11_Priority', 'Unknown'),
                    'literature_publications': 0,
                    'has_literature_evidence': False,
                    'total_trials': 0,
                    'active_trials': 0,
                    'recruiting_trials': 0,
                    'has_clinical_trials': False,
                    'has_tool_compounds': False,
                    'tool_compound_count': 0,
                    'probe_quality': 'Not Assessed',
                    'has_combination_potential': False,
                    'combination_strategies': [],
                    'synthetic_lethality_candidates': [],
                    'pathway_combinations': [],
                    'validation_score': 0,
                    'validation_level': 'Not Assessed'
                }
                validation_data.update(ground_truth)
            
            # Step12 columns
            result = {
                gene_col: gene,
                'Step12_Literature_Publications': validation_data.get('literature_publications', 0),
                'Step12_Has_Literature_Evidence': validation_data.get('has_literature_evidence', False),
                'Step12_Total_Clinical_Trials': validation_data.get('total_trials', 0),
                'Step12_Active_Trials': validation_data.get('active_trials', 0),
                'Step12_Recruiting_Trials': validation_data.get('recruiting_trials', 0),
                'Step12_Has_Clinical_Trials': validation_data.get('has_clinical_trials', False),
                'Step12_Has_Tool_Compounds': validation_data.get('has_tool_compounds', False),
                'Step12_Tool_Compound_Count': validation_data.get('tool_compound_count', 0),
                'Step12_Probe_Quality': validation_data.get('probe_quality', 'Not Assessed'),
                'Step12_Has_Combination_Potential': validation_data.get('has_combination_potential', False),
                'Step12_Validation_Score': validation_data.get('validation_score', 0),
                'Step12_Validation_Level': validation_data.get('validation_level', 'Not Assessed'),
                'Step12_Combination_Strategies': '; '.join(validation_data.get('combination_strategies', [])[:3]),
                'Step12_Synthetic_Lethality': '; '.join(validation_data.get('synthetic_lethality_candidates', [])[:2]),
                'Step12_Pathway_Combinations': '; '.join(validation_data.get('pathway_combinations', [])[:2]),
                'SM_Literature_Score': ground_truth['SM_Literature_Score'],
                'Bio_Literature_Score': ground_truth['Bio_Literature_Score'],
                'SM_Is_Druggable': ground_truth['SM_Is_Druggable'],
                'Bio_Is_Druggable': ground_truth['Bio_Is_Druggable']
            }
            
            all_validation_data.append(result)
            
            # Generate final prediction
            prediction = self.generate_final_prediction(row, validation_data)
            all_predictions.append(prediction)
            
            if idx % 10 == 0:
                logger.info(f"  Processed {idx}/{len(step11_df)} genes...")
        
        validation_df = pd.DataFrame(all_validation_data)
        predictions_df = pd.DataFrame(all_predictions)
        
        combined_df = step11_df.merge(validation_df, on=gene_col, how='left')
        
        logger.info(f"\nAdded validation data for all {len(step11_df)} candidates")
        return combined_df, predictions_df
    
    def generate_verification_protocols(self, step11_df: pd.DataFrame, output_dir: str = 'verification_protocols'):
        """Generate individual verification protocol files for each high-priority target"""
        
        Path(output_dir).mkdir(exist_ok=True)
        
        gene_col = self._find_gene_column(step11_df)
        high_priority = step11_df[step11_df['Step11_Priority'] == 'High Priority']
        
        if len(high_priority) == 0:
            logger.warning("No high-priority candidates to generate protocols for")
            return
        
        logger.info(f"\nGenerating verification protocols for {len(high_priority)} targets...")
        
        for _, row in high_priority.iterrows():
            gene = row[gene_col]
            protocol = self.generate_literature_verification_protocol(row, gene)
            
            output_file = Path(output_dir) / f"{gene}_verification_protocol.txt"
            
            with open(output_file, 'w') as f:
                f.write("="*80 + "\n")
                f.write(f"LITERATURE VERIFICATION PROTOCOL: {gene}\n")
                f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("="*80 + "\n\n")
                
                f.write(f"PREDICTED MODALITY: {protocol['predicted_modality']}\n\n")
                
                f.write("LITERATURE SEARCH QUERIES:\n")
                for query in protocol['literature_search_queries']:
                    f.write(f"  {query}\n")
                f.write("\n")
                
                f.write("KEY EVIDENCE TO FIND:\n")
                for evidence in protocol['key_evidence_to_find']:
                    f.write(f"  {evidence}\n")
                f.write("\n")
                
                f.write("VALIDATION EXPERIMENTS:\n")
                for experiment in protocol['validation_experiments']:
                    f.write(f"  {experiment}\n")
                f.write("\n")
                
                f.write("RED FLAGS:\n")
                for flag in protocol['red_flags_to_check']:
                    f.write(f"   {flag}\n")
                f.write("\n")
                
                f.write("EXPECTED DRUGGABILITY MARKERS:\n")
                for marker in protocol['expected_druggability_markers']:
                    f.write(f"   {marker}\n")
                f.write("\n")
            
            logger.info(f"  Created protocol: {output_file}")
        
        logger.info(f"? Generated {len(high_priority)} verification protocols in {output_dir}/")
    
    def generate_validation_report(self, df: pd.DataFrame, output_file: str = 'Step12_validation_report.txt'):
        """Generate detailed validation report"""
        
        gene_col = self._find_gene_column(df)
        validated = df[df['Step12_Validation_Level'].notna()].sort_values(
            'Step12_Validation_Score', ascending=False
        )
        
        if len(validated) == 0:
            logger.warning("No validated candidates to report")
            return
        
        with open(output_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("HIGH-PRIORITY CANDIDATE VALIDATION REPORT\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Validated Candidates: {len(validated)}\n")
            f.write("="*80 + "\n\n")
            
            for _, row in validated.iterrows():
                gene = row[gene_col]
                
                f.write(f"\nTARGET: {gene}\n")
                f.write("-"*80 + "\n\n")
                
                f.write("SCORES:\n")
                f.write(f"  Composite: {row.get('Step11_Composite_Score', 0):.1f}/100\n")
                f.write(f"  Small Molecule: {row.get('Step11_SmallMolecule_Score', 0):.1f}/100 (Literature: {row.get('SM_Literature_Score', 0):.0f}/100)\n")
                f.write(f"  Biologic: {row.get('Step11_Biologic_Score', 0):.1f}/100 (Literature: {row.get('Bio_Literature_Score', 0):.0f}/100)\n")
                f.write(f"  Best Modality: {row.get('Step11_Best_Modality', 'Unknown')}\n\n")
                
                f.write("VALIDATION:\n")
                f.write(f"  Level: {row.get('Step12_Validation_Level', 'Unknown')}\n")
                f.write(f"  Score: {row.get('Step12_Validation_Score', 0)}/10\n")
                f.write(f"  Publications: {row.get('Step12_Literature_Publications', 0)}\n")
                f.write(f"  Clinical Trials: {row.get('Step12_Total_Clinical_Trials', 0)}\n\n")
        
        logger.info(f"Validation report saved: {output_file}")


def main():
    """Main execution function"""
    
    STEP11_INPUT = 'Step11.csv'
    STEP12_OUTPUT = 'Step12.csv'
    PREDICTIONS_OUTPUT = 'Step12_final_predictions.csv'
    DISEASE = "migraine"
    
    try:
        logger.info("="*80)
        logger.info("STEP 12: CANDIDATE VALIDATION & PERFORMANCE METRICS")
        logger.info("="*80)
        
        validator = CandidateValidator()
        step11_df = validator.load_step11_data(STEP11_INPUT)
        
        logger.info(f"\nDisease context: {DISEASE}")
        
        # Add validation data for all genes
        step12_df, predictions_df = validator.add_validation_data(step11_df, disease=DISEASE)
        
        # Calculate performance metrics
        metrics = validator.calculate_performance_metrics(step12_df)
        
        # Generate metrics report
        validator.generate_metrics_report(metrics)
        
        logger.info(f"\nSaving to {STEP12_OUTPUT}...")
        step12_df.to_csv(STEP12_OUTPUT, index=False)
        logger.info(f"Saved {len(step12_df)} rows")
        
        if len(predictions_df) > 0:
            logger.info(f"\nSaving final predictions to {PREDICTIONS_OUTPUT}...")
            predictions_df.to_csv(PREDICTIONS_OUTPUT, index=False)
            logger.info(f"Saved {len(predictions_df)} predictions")
        
        validator.generate_validation_report(step12_df)
        validator.generate_verification_protocols(step12_df)
        
        logger.info("\n" + "="*80)
        logger.info("STEP 12 COMPLETE!")
        logger.info("="*80)
        logger.info(f"\nGenerated files:")
        logger.info(f"   {STEP12_OUTPUT} - Complete validation data with ground truth")
        logger.info(f"   {PREDICTIONS_OUTPUT} - Final predictions")
        logger.info(f"   Step12_performance_metrics.txt - AUC, Accuracy, Correlation")
        logger.info(f"   Step12_validation_report.txt - Validation summary")
        logger.info(f"   verification_protocols/ - Individual verification protocols")
        
        # Print summary of metrics
        if metrics:
            logger.info("\n" + "="*80)
            logger.info("PERFORMANCE SUMMARY")
            logger.info("="*80)
            
            if 'Small_Molecule' in metrics and 'Error' not in metrics['Small_Molecule']:
                sm = metrics['Small_Molecule']
                logger.info(f"\nSmall Molecule:")
                logger.info(f"  AUC: {sm['AUC']:.3f}")
                logger.info(f"  Accuracy: {sm['Accuracy']:.3f} ({sm['Accuracy']*100:.1f}%)")
                logger.info(f"  Pearson r: {sm['Pearson_r']:.3f}")
            
            if 'Biologic' in metrics and 'Error' not in metrics['Biologic']:
                bio = metrics['Biologic']
                logger.info(f"\nBiologic:")
                logger.info(f"  AUC: {bio['AUC']:.3f}")
                logger.info(f"  Accuracy: {bio['Accuracy']:.3f} ({bio['Accuracy']*100:.1f}%)")
                logger.info(f"  Pearson r: {bio['Pearson_r']:.3f}")
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
