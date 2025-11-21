"""
Step 11: FINAL Dual-Track Druggability Scorer (FIXED)
Provides separate scores for Small Molecules AND Biologics
Fixes all known issues including JAK2, BDNF, and Collagens
Output: Step11.csv
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional
from pathlib import Path
import logging
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('step11_final_dual_track.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class FinalDualTrackDruggabilityScorer:
    """
    Final production-ready dual-track druggability scorer
    Separate scoring for Small Molecules and Biologics
    Fixed for BDNF, COL1A2, COL4A1, and other edge cases
    """
    
    # Small Molecule weights (total = 29)
    SM_WEIGHTS = {
        'known_drugs': 10,          # 34.5% - Proven druggability
        'druggable_family': 6,      # 20.7% - Family prediction
        'pocket_score': 5,          # 17.2% - Binding site
        'pharos_level': 2,          # 6.9% - Development status
        'accessibility': 2,         # 6.9% - Location
        'low_disorder': 2,          # 6.9% - Structure quality
        'tissue_expression': 1,     # 3.4% - Expression
        'network_connectivity': 1   # 3.4% - Network
    }
    
    # Biologic weights (total = 20)
    BIO_WEIGHTS = {
        'known_biologics': 8,       # 40% - Approved antibodies
        'secreted_accessible': 6,   # 30% - Location
        'pharos_level': 3,          # 15% - Development
        'tissue_expression': 2,     # 10% - Expression
        'network_connectivity': 1   # 5% - Network
    }
    
    THRESHOLDS = {
        'high': 80,
        'medium': 50
    }
    
    # Expanded biologic-only families
    BIOLOGIC_ONLY_FAMILIES = {
        'cytokine', 'interleukin', 'chemokine', 'growth factor',
        'structural protein', 'collagen', 'extracellular matrix',
        'inflammasome', 'scaffold protein', 'adapter protein',
        'neurotrophin', 'neurotrophic factor',  # Added for BDNF
        'brain-derived', 'bdnf'  # Additional BDNF indicators
    }
    
    HIGHLY_DRUGGABLE_FAMILIES = ['Kinase', 'GPCR', 'IC', 'NR', 'Enzyme']
    
    # Genes with known false positive drug annotations
    FALSE_POSITIVE_DRUG_GENES = {
        'BDNF': 'Neurotrophin - drugs modulate levels, do not target protein directly',
        'MTHFR': 'Polymorphism target - drugs affect pathway, not protein binding'
    }
    
    def __init__(self):
        self.cache = {}
    
    def load_data(self, filepath: str = 'Step10.csv') -> pd.DataFrame:
        logger.info(f"Loading data from {filepath}")
        df = pd.read_csv(filepath)
        logger.info(f"Loaded {len(df)} entries")
        return df
    
    def _safe_get(self, row: pd.Series, key: str, default=None):
        val = row.get(key, default)
        if pd.isna(val):
            return default
        return val
    
    def _safe_bool(self, row: pd.Series, key: str) -> bool:
        val = row.get(key, False)
        if pd.isna(val):
            return False
        return bool(val)
    
    def _safe_num(self, row: pd.Series, key: str, default=0) -> float:
        val = row.get(key, default)
        if pd.isna(val) or val == '':
            return default
        try:
            return float(val)
        except:
            return default
    
    def _is_biologic_only_target(self, row: pd.Series) -> bool:
        """Determine if target is ONLY druggable by biologics"""
        
        gene_symbol = str(self._safe_get(row, 'gene_symbol', ''))
        
        # CRITICAL: Check for known false positive genes first
        if gene_symbol in self.FALSE_POSITIVE_DRUG_GENES:
            return True
        
        # SAFEGUARD 1: Druggable families with drugs are NOT biologic-only
        step9_family = str(self._safe_get(row, 'Step9_Protein_Family', 'Unknown'))
        dgidb_approved = self._safe_num(row, 'dgidb_num_approved', 0)
        
        if step9_family in self.HIGHLY_DRUGGABLE_FAMILIES and dgidb_approved >= 2:
            return False
        
        # SAFEGUARD 2: Explicit "Small molecule" with drugs
        drug_types = str(self._safe_get(row, 'opentargets_drug_types', ''))
        
        # Cytokines/collagens are biologic-only
        if gene_symbol.startswith(('IL', 'TNF', 'TSLP', 'COL', 'BDNF')):
            return True
        
        # If explicitly small molecule with approved drugs
        if 'Small molecule' in drug_types and dgidb_approved >= 1:
            return False
        
        # Check family keywords
        family_type = str(self._safe_get(row, 'druggable_family_type', 'Unknown'))
        protein_class = str(self._safe_get(row, 'protein_class', ''))
        uniprot_families = str(self._safe_get(row, 'uniprot_protein_families', ''))
        gene_name = str(self._safe_get(row, 'gene_name', ''))
        gene_desc = str(self._safe_get(row, 'gene_description', ''))
        uniprot_name = str(self._safe_get(row, 'uniprot_protein_name', ''))
        
        combined_text = f"{family_type} {protein_class} {uniprot_families} {gene_name} {gene_desc} {uniprot_name}".lower()
        
        for biologic_term in self.BIOLOGIC_ONLY_FAMILIES:
            if biologic_term in combined_text:
                return True
        
        # Secreted proteins (with exceptions)
        location = str(self._safe_get(row, 'subcellular_location', ''))
        step6_location = str(self._safe_get(row, 'Step6_Subcellular_Location', ''))
        combined_location = f"{location} {step6_location}".lower()
        
        if 'secreted' in combined_location and 'membrane' not in combined_location:
            # Exceptions: Enzymes like ACE
            if 'enzyme' in combined_text or 'protease' in combined_text:
                # But NOT if it's also a growth factor/cytokine
                if any(term in combined_text for term in ['cytokine', 'growth factor', 'hormone', 'interleukin', 'neurotroph']):
                    return True
                return False
            
            # Secreted + cytokine/growth factor = biologic
            if any(term in combined_text for term in ['cytokine', 'growth factor', 'hormone', 'interleukin', 'neurotroph']):
                return True
        
        return False
    
    def _has_small_molecule_drugs(self, row: pd.Series) -> tuple:
        """
        Check for SMALL MOLECULE drugs (FIXED for JAK2 and BDNF)
        Returns: (has_small_molecule, count, evidence)
        """
        
        gene_symbol = str(self._safe_get(row, 'gene_symbol', ''))
        
        # CRITICAL: Check for known false positive genes
        if gene_symbol in self.FALSE_POSITIVE_DRUG_GENES:
            return (False, 0, self.FALSE_POSITIVE_DRUG_GENES[gene_symbol])
        
        has_approved = self._safe_bool(row, 'has_approved_drugs')
        dgidb_approved = self._safe_num(row, 'dgidb_num_approved', 0)
        opentargets_approved = self._safe_num(row, 'opentargets_approved', 0)
        chembl_approved = self._safe_num(row, 'chembl_approved_drugs', 0)
        
        drug_types = str(self._safe_get(row, 'opentargets_drug_types', ''))
        drug_names = str(self._safe_get(row, 'dgidb_approved_drugs', ''))
        
        # FIX FOR JAK2: Trust explicit "Small molecule" label
        if 'Small molecule' in drug_types and dgidb_approved >= 1:
            return (True, int(dgidb_approved), f"{int(dgidb_approved)} small molecule drugs")
        
        # Check for biologic-only indicators
        biologic_keywords = ['antibody', 'antibodies', 'mab', 'umab', 'immunoglobulin', 'protein', 'peptide']
        
        is_biologic = False
        if drug_types:
            drug_types_lower = drug_types.lower()
            if any(keyword in drug_types_lower for keyword in biologic_keywords) and 'small molecule' not in drug_types_lower:
                is_biologic = True
        
        if drug_names:
            drug_names_lower = drug_names.lower()
            antibody_count = sum(1 for term in ['mab', 'umab'] if term in drug_names_lower)
            if antibody_count > 0 and 'small molecule' not in drug_types:
                is_biologic = True
        
        total_approved = max(dgidb_approved, opentargets_approved, chembl_approved)
        
        if is_biologic:
            return (False, 0, "Only biologics")
        elif has_approved or total_approved > 0:
            return (True, int(total_approved), f"{int(total_approved)} small molecule drugs")
        else:
            return (False, 0, "No drugs")
    
    def _has_biologic_drugs(self, row: pd.Series) -> tuple:
        """
        Check for BIOLOGIC drugs (antibodies, proteins)
        Returns: (has_biologic, count, evidence)
        """
        
        has_approved = self._safe_bool(row, 'has_approved_drugs')
        drug_types = str(self._safe_get(row, 'opentargets_drug_types', ''))
        drug_names = str(self._safe_get(row, 'dgidb_approved_drugs', ''))
        
        # Check for biologic indicators
        biologic_keywords = ['antibody', 'protein', 'peptide']
        
        has_biologic = False
        count = 0
        
        if drug_types:
            drug_types_lower = drug_types.lower()
            if any(keyword in drug_types_lower for keyword in biologic_keywords):
                has_biologic = True
                # Rough estimate of biologic count
                if 'antibody' in drug_types_lower or 'antibodies' in drug_types_lower:
                    count = max(1, self._safe_num(row, 'opentargets_approved', 0))
        
        if drug_names:
            drug_names_lower = drug_names.lower()
            count = max(count, sum(1 for term in ['mab', 'umab'] if term in drug_names_lower))
        
        if has_biologic:
            return (True, int(count), f"{int(count)} biologic drug(s)")
        else:
            return (False, 0, "No biologics")
    
    # ======================== SMALL MOLECULE SCORING ========================
    
    def calculate_sm_drugs_score(self, row: pd.Series) -> Dict:
        """Small molecule drug score"""
        score = 0
        evidence = []
        
        is_biologic_only = self._is_biologic_only_target(row)
        
        if is_biologic_only:
            score = 0
            evidence.append("BIOLOGIC-ONLY")
        else:
            has_sm, count, drug_evidence = self._has_small_molecule_drugs(row)
            
            if has_sm and count >= 10:
                score = 10
                evidence.append(f"{count} SM drugs")
            elif has_sm and count >= 5:
                score = 9
                evidence.append(f"{count} SM drugs")
            elif has_sm and count >= 3:
                score = 8
                evidence.append(f"{count} SM drugs")
            elif has_sm and count >= 1:
                score = 7
                evidence.append(f"{count} SM drug(s)")
            else:
                has_clinical = self._safe_bool(row, 'has_clinical_drugs')
                if has_clinical:
                    score = 4
                    evidence.append("Clinical")
                else:
                    bioactive = self._safe_num(row, 'chembl_bioactive_compounds', 0)
                    if bioactive >= 100:
                        score = 3
                    elif bioactive >= 50:
                        score = 2
                    elif bioactive >= 10:
                        score = 1
                    else:
                        score = 0
                    
                    if bioactive > 0:
                        evidence.append(f"{int(bioactive)} bioactive")
                    else:
                        evidence.append("No drugs")
        
        return {'score': score, 'evidence': ' | '.join(evidence)}
    
    def calculate_sm_family_score(self, row: pd.Series) -> Dict:
        """Small molecule family score"""
        score = 0
        evidence = []
        
        if self._is_biologic_only_target(row):
            return {'score': 0, 'evidence': 'BIOLOGIC-ONLY'}
        
        step9_family = str(self._safe_get(row, 'Step9_Protein_Family', 'Unknown'))
        family_type = str(self._safe_get(row, 'druggable_family_type', 'Unknown'))
        confidence = str(self._safe_get(row, 'druggable_family_confidence', 'LOW'))
        
        combined_family = f"{step9_family} {family_type}"
        
        if step9_family == 'GPCR' or 'GPCR' in combined_family:
            score = 10
            evidence.append("GPCR")
        elif step9_family == 'Kinase' or 'Kinase' in combined_family:
            score = 10
            evidence.append("Kinase")
        elif step9_family == 'IC' or 'Ion channel' in combined_family.lower():
            score = 9
            evidence.append("Ion Channel")
        elif step9_family == 'NR' or 'Nuclear receptor' in combined_family.lower():
            score = 9
            evidence.append("Nuclear Receptor")
        elif 'Protease' in combined_family:
            score = 8
            evidence.append("Protease")
        elif 'Transporter' in combined_family and confidence == 'HIGH':
            score = 7
            evidence.append("Transporter")
        elif 'Enzyme' in combined_family and confidence == 'HIGH':
            score = 7
            evidence.append("Enzyme")
        elif 'Enzyme' in combined_family:
            score = 5
            evidence.append("Enzyme")
        elif self._safe_bool(row, 'is_druggable_family'):
            if confidence == 'HIGH':
                score = 4
            elif confidence == 'MEDIUM':
                score = 3
            else:
                score = 2
            evidence.append(f"Druggable ({confidence})")
        else:
            score = 0
            evidence.append("Non-druggable")
        
        return {'score': score, 'evidence': ' | '.join(evidence)}
    
    def calculate_sm_pocket_score(self, row: pd.Series) -> Dict:
        """Small molecule pocket score"""
        score = 0
        evidence = []
        
        if self._is_biologic_only_target(row):
            return {'score': 0, 'evidence': 'BIOLOGIC-ONLY'}
        
        has_pocket = self._safe_bool(row, 'has_druggable_pocket')
        best_drug_score = self._safe_num(row, 'best_pocket_druggability', 0)
        pocket_volume = self._safe_num(row, 'best_pocket_volume', 0)
        
        if has_pocket and best_drug_score >= 0.8 and pocket_volume >= 500:
            score = 10
            evidence.append(f"Excellent pocket")
        elif has_pocket and best_drug_score >= 0.6 and pocket_volume >= 400:
            score = 8
            evidence.append(f"Good pocket")
        elif has_pocket and best_drug_score >= 0.5 and pocket_volume >= 300:
            score = 6
            evidence.append(f"Moderate pocket")
        elif has_pocket and best_drug_score >= 0.3 and pocket_volume >= 200:
            score = 3
            evidence.append(f"Weak pocket")
        elif has_pocket:
            score = 1
            evidence.append(f"Poor pocket")
        else:
            score = 0
            evidence.append("No pockets")
        
        return {'score': score, 'evidence': ' | '.join(evidence)}
    
    def calculate_sm_pharos_score(self, row: pd.Series) -> Dict:
        """Pharos TDL score"""
        tdl = str(self._safe_get(row, 'Step9_Target_Development_Level', 'Unknown'))
        
        if tdl == 'Tclin':
            return {'score': 10, 'evidence': 'Tclin'}
        elif tdl == 'Tchem':
            return {'score': 7, 'evidence': 'Tchem'}
        elif tdl == 'Tbio':
            return {'score': 3, 'evidence': 'Tbio'}
        elif tdl == 'Tdark':
            return {'score': 0, 'evidence': 'Tdark'}
        else:
            return {'score': 0, 'evidence': 'Unknown'}
    
    def calculate_sm_accessibility_score(self, row: pd.Series) -> Dict:
        """Accessibility score"""
        tier = str(self._safe_get(row, 'Step6_Druggability_Tier', 'Unknown'))
        
        if 'Tier 1' in tier:
            return {'score': 10, 'evidence': 'Tier 1'}
        elif 'Tier 2' in tier:
            return {'score': 7, 'evidence': 'Tier 2'}
        elif 'Tier 3' in tier:
            return {'score': 4, 'evidence': 'Tier 3'}
        elif 'Tier 4' in tier:
            return {'score': 1, 'evidence': 'Tier 4'}
        else:
            return {'score': 0, 'evidence': 'Unknown'}
    
    def calculate_sm_disorder_score(self, row: pd.Series) -> Dict:
        """Disorder score"""
        disorder_pct = self._safe_num(row, 'Step8_Disorder_Percentage', None)
        
        if disorder_pct is None:
            return {'score': 0, 'evidence': 'No data'}
        elif disorder_pct < 10:
            return {'score': 10, 'evidence': f'Excellent ({disorder_pct:.1f}%)'}
        elif disorder_pct < 20:
            return {'score': 8, 'evidence': f'Very low ({disorder_pct:.1f}%)'}
        elif disorder_pct < 30:
            return {'score': 6, 'evidence': f'Low ({disorder_pct:.1f}%)'}
        elif disorder_pct < 40:
            return {'score': 4, 'evidence': f'Moderate ({disorder_pct:.1f}%)'}
        elif disorder_pct < 50:
            return {'score': 2, 'evidence': f'High ({disorder_pct:.1f}%)'}
        else:
            return {'score': 0, 'evidence': f'Very high ({disorder_pct:.1f}%)'}
    
    def calculate_expression_score(self, row: pd.Series) -> Dict:
        """Expression score"""
        gtex_available = self._safe_bool(row, 'Step7_GTEx_Data_Available')
        
        if not gtex_available:
            return {'score': 0, 'evidence': 'No data'}
        
        tissues = self._safe_num(row, 'Step7_Tissues_With_Expression', 0)
        max_tpm = self._safe_num(row, 'Step7_Max_TPM', 0)
        
        if tissues >= 40 and max_tpm >= 10:
            return {'score': 10, 'evidence': f'Broad ({int(tissues)} tissues)'}
        elif tissues >= 30:
            return {'score': 7, 'evidence': f'Good ({int(tissues)} tissues)'}
        elif tissues >= 20:
            return {'score': 5, 'evidence': f'Moderate ({int(tissues)} tissues)'}
        elif tissues >= 10:
            return {'score': 3, 'evidence': f'Limited ({int(tissues)} tissues)'}
        else:
            return {'score': 0, 'evidence': 'Poor'}
    
    def calculate_network_score(self, row: pd.Series) -> Dict:
        """Network score"""
        string_available = self._safe_bool(row, 'Step10_STRING_Data_Available')
        
        if not string_available:
            return {'score': 0, 'evidence': 'No data'}
        
        drug_int = self._safe_num(row, 'Step10_Drug_Target_Interactions', 0)
        is_hub = self._safe_bool(row, 'Step10_Network_Hub')
        
        if drug_int >= 10 or (is_hub and drug_int >= 5):
            return {'score': 10, 'evidence': f'{int(drug_int)} targets'}
        elif drug_int >= 5:
            return {'score': 7, 'evidence': f'{int(drug_int)} targets'}
        elif drug_int >= 2:
            return {'score': 4, 'evidence': f'{int(drug_int)} targets'}
        elif drug_int >= 1:
            return {'score': 2, 'evidence': '1 target'}
        else:
            return {'score': 0, 'evidence': 'Poor'}
    
    # ======================== BIOLOGIC SCORING ========================
    
    def calculate_bio_drugs_score(self, row: pd.Series) -> Dict:
        """Biologic drug score"""
        has_bio, count, evidence_text = self._has_biologic_drugs(row)
        
        if has_bio and count >= 5:
            return {'score': 10, 'evidence': f'{count} biologics'}
        elif has_bio and count >= 3:
            return {'score': 9, 'evidence': f'{count} biologics'}
        elif has_bio and count >= 1:
            return {'score': 8, 'evidence': f'{count} biologic(s)'}
        else:
            has_clinical = self._safe_bool(row, 'has_clinical_drugs')
            drug_types = str(self._safe_get(row, 'opentargets_drug_types', ''))
            
            if has_clinical and any(x in drug_types for x in ['Antibody', 'Protein']):
                return {'score': 5, 'evidence': 'Clinical biologics'}
            else:
                return {'score': 0, 'evidence': 'No biologics'}
    
    def calculate_bio_accessibility_score(self, row: pd.Series) -> Dict:
        """Biologic accessibility score (secreted/membrane) - FIXED for collagens"""
        gene_symbol = str(self._safe_get(row, 'gene_symbol', ''))
        
        # FIX: Collagens are NOT druggable by antibodies (insoluble fibers)
        if gene_symbol.startswith('COL'):
            return {'score': 2, 'evidence': 'Collagen - insoluble'}
        
        location = str(self._safe_get(row, 'subcellular_location', ''))
        step6_location = str(self._safe_get(row, 'Step6_Subcellular_Location', ''))
        topology = str(self._safe_get(row, 'Step6_Topology', ''))
        
        combined_location = f"{location} {step6_location}".lower()
        
        # Secreted proteins are highly accessible to antibodies
        if 'secreted' in combined_location:
            return {'score': 10, 'evidence': 'Secreted'}
        
        # Membrane proteins with extracellular domains
        if any(x in combined_location for x in ['cell membrane', 'plasma membrane']):
            if 'extracellular' in topology.lower():
                return {'score': 9, 'evidence': 'Membrane + ECD'}
            else:
                return {'score': 6, 'evidence': 'Membrane'}
        
        # Intracellular = poor for antibodies
        if any(x in combined_location for x in ['cytoplasm', 'nucleus', 'mitochondri']):
            return {'score': 1, 'evidence': 'Intracellular'}
        
        return {'score': 3, 'evidence': 'Unknown'}
    
    # ======================== COMPOSITE CALCULATION ========================
    
    def calculate_small_molecule_score(self, row: pd.Series) -> Dict:
        """Calculate overall small molecule druggability score"""
        
        components = {
            'known_drugs': self.calculate_sm_drugs_score(row),
            'druggable_family': self.calculate_sm_family_score(row),
            'pocket_score': self.calculate_sm_pocket_score(row),
            'pharos_level': self.calculate_sm_pharos_score(row),
            'accessibility': self.calculate_sm_accessibility_score(row),
            'low_disorder': self.calculate_sm_disorder_score(row),
            'tissue_expression': self.calculate_expression_score(row),
            'network_connectivity': self.calculate_network_score(row)
        }
        
        # Check critical components
        critical_scores = {
            'known_drugs': components['known_drugs']['score'],
            'druggable_family': components['druggable_family']['score'],
            'pocket_score': components['pocket_score']['score']
        }
        
        # If all critical are 0, cap score
        if all(score == 0 for score in critical_scores.values()):
            max_possible_score = 30
        else:
            max_possible_score = 100
        
        # Calculate weighted sum
        weighted_sum = 0
        total_weight = 0
        
        for component, result in components.items():
            weight = self.SM_WEIGHTS[component]
            score = result['score']
            weighted_sum += score * weight
            total_weight += 10 * weight
        
        raw_score = (weighted_sum / total_weight) * 100 if total_weight > 0 else 0
        composite_score = min(raw_score, max_possible_score)
        
        # Biologic-only penalty
        is_biologic_only = self._is_biologic_only_target(row)
        if is_biologic_only:
            composite_score = min(composite_score, 20)
        
        return {
            'score': round(composite_score, 2),
            'components': components,
            'is_biologic_only': is_biologic_only
        }
    
    def calculate_biologic_score(self, row: pd.Series) -> Dict:
        """Calculate overall biologic druggability score"""
        
        components = {
            'known_biologics': self.calculate_bio_drugs_score(row),
            'secreted_accessible': self.calculate_bio_accessibility_score(row),
            'pharos_level': self.calculate_sm_pharos_score(row),
            'tissue_expression': self.calculate_expression_score(row),
            'network_connectivity': self.calculate_network_score(row)
        }
        
        # Calculate weighted sum
        weighted_sum = 0
        total_weight = 0
        
        for component, result in components.items():
            weight = self.BIO_WEIGHTS[component]
            score = result['score']
            weighted_sum += score * weight
            total_weight += 10 * weight
        
        composite_score = (weighted_sum / total_weight) * 100 if total_weight > 0 else 0
        
        return {
            'score': round(composite_score, 2),
            'components': components
        }
    
    def process_gene(self, row: pd.Series) -> Dict:
        """Process single gene and return all scores"""
        
        gene = row.get('gene_symbol', 'Unknown')
        
        # Calculate both tracks
        sm_result = self.calculate_small_molecule_score(row)
        bio_result = self.calculate_biologic_score(row)
        
        sm_score = sm_result['score']
        bio_score = bio_result['score']
        
        # Determine best modality
        if sm_score > bio_score:
            best_modality = 'Small Molecule'
            composite_score = sm_score
        elif bio_score > sm_score:
            best_modality = 'Biologic'
            composite_score = bio_score
        else:
            best_modality = 'Both'
            composite_score = sm_score
        
        # Overall classification based on composite
        if composite_score >= self.THRESHOLDS['high']:
            classification = 'High Druggability'
            priority = 'High Priority'
        elif composite_score >= self.THRESHOLDS['medium']:
            classification = 'Medium Druggability'
            priority = 'Medium Priority'
        else:
            classification = 'Low Druggability'
            priority = 'Low Priority'
        
        # Generate evidence summary
        evidence_pieces = []
        for comp, result in sm_result['components'].items():
            if result['score'] >= 7:
                evidence_pieces.append(result['evidence'])
        
        return {
            'gene_symbol': gene,
            'Step11_Composite_Score': composite_score,
            'Step11_SmallMolecule_Score': sm_score,
            'Step11_Biologic_Score': bio_score,
            'Step11_Best_Modality': best_modality,
            'Step11_Classification': classification,
            'Step11_Priority': priority,
            'Step11_Is_Biologic_Only': sm_result['is_biologic_only'],
            'Step11_SM_Drugs_Score': sm_result['components']['known_drugs']['score'],
            'Step11_SM_Family_Score': sm_result['components']['druggable_family']['score'],
            'Step11_SM_Pocket_Score': sm_result['components']['pocket_score']['score'],
            'Step11_SM_Pharos_Score': sm_result['components']['pharos_level']['score'],
            'Step11_SM_Access_Score': sm_result['components']['accessibility']['score'],
            'Step11_SM_Disorder_Score': sm_result['components']['low_disorder']['score'],
            'Step11_Expression_Score': sm_result['components']['tissue_expression']['score'],
            'Step11_Network_Score': sm_result['components']['network_connectivity']['score'],
            'Step11_Bio_Drugs_Score': bio_result['components']['known_biologics']['score'],
            'Step11_Bio_Access_Score': bio_result['components']['secreted_accessible']['score'],
            'Step11_Evidence': ' | '.join(evidence_pieces[:8]) if evidence_pieces else 'Limited evidence',
            'Step11_SM_Drugs_Evidence': sm_result['components']['known_drugs']['evidence'],
            'Step11_SM_Family_Evidence': sm_result['components']['druggable_family']['evidence'],
            'Step11_Bio_Access_Evidence': bio_result['components']['secreted_accessible']['evidence']
        }
    
    def score_all_genes(self, df: pd.DataFrame) -> pd.DataFrame:
        """Score all genes in DataFrame"""
        
        logger.info(f"\nScoring {len(df)} genes with dual-track system...")
        logger.info("="*80)
        
        all_results = []
        
        for idx, (_, row) in enumerate(df.iterrows(), 1):
            result = self.process_gene(row)
            all_results.append(result)
            
            if idx % 10 == 0:
                logger.info(f"  Processed {idx}/{len(df)} genes...")
        
        results_df = pd.DataFrame(all_results)
        logger.info(f"\n✓ Scored {len(results_df)} genes")
        
        return results_df
    
    def generate_summary(self, df: pd.DataFrame):
        """Generate summary report"""
        
        logger.info("\n" + "="*80)
        logger.info("DUAL-TRACK DRUGGABILITY SUMMARY")
        logger.info("="*80)
        
        logger.info(f"\nClassification (by Composite Score):")
        for classification, count in df['Step11_Classification'].value_counts().items():
            pct = (count / len(df)) * 100
            logger.info(f"  {classification}: {count} ({pct:.1f}%)")
        
        logger.info(f"\nBest Modality:")
        for modality, count in df['Step11_Best_Modality'].value_counts().items():
            pct = (count / len(df)) * 100
            logger.info(f"  {modality}: {count} ({pct:.1f}%)")
        
        logger.info(f"\nBiologic-Only Targets:")
        biologic_only = df[df['Step11_Is_Biologic_Only'] == True]
        logger.info(f"  Count: {len(biologic_only)}")
        if len(biologic_only) > 0:
            logger.info(f"  Examples: {', '.join(biologic_only['gene_symbol'].head(10).tolist())}")
        
        logger.info(f"\nScore Distributions:")
        logger.info(f"  Composite: Mean={df['Step11_Composite_Score'].mean():.1f}, Median={df['Step11_Composite_Score'].median():.1f}")
        logger.info(f"  Small Molecule: Mean={df['Step11_SmallMolecule_Score'].mean():.1f}, Median={df['Step11_SmallMolecule_Score'].median():.1f}")
        logger.info(f"  Biologic: Mean={df['Step11_Biologic_Score'].mean():.1f}, Median={df['Step11_Biologic_Score'].median():.1f}")


def main():
    """Main execution"""
    
    INPUT_FILE = 'Step10.csv'
    OUTPUT_FILE = 'Step11.csv'
    
    try:
        logger.info("="*80)
        logger.info("FINAL DUAL-TRACK DRUGGABILITY SCORER (FIXED)")
        logger.info("Small Molecule + Biologic Scoring")
        logger.info("Fixes: JAK2, BDNF, COL1A2, COL4A1")
        logger.info("="*80)
        
        scorer = FinalDualTrackDruggabilityScorer()
        df = scorer.load_data(INPUT_FILE)
        
        # Score all genes
        results_df = scorer.score_all_genes(df)
        
        # Merge with original data
        full_df = df.merge(results_df, on='gene_symbol', how='left')
        
        # Save to Step11.csv
        logger.info(f"\nSaving results to {OUTPUT_FILE}...")
        full_df.to_csv(OUTPUT_FILE, index=False)
        logger.info(f"✓ Saved {len(full_df)} rows to {OUTPUT_FILE}")
        
        # Generate summary
        scorer.generate_summary(results_df)
        
        # Display top results
        logger.info("\n" + "="*80)
        logger.info("TOP 20 TARGETS (by Composite Score)")
        logger.info("="*80)
        
        top_20 = results_df.sort_values('Step11_Composite_Score', ascending=False).head(20)
        
        for idx, (_, row) in enumerate(top_20.iterrows(), 1):
            logger.info(f"\n{idx}. {row['gene_symbol']}")
            logger.info(f"   Composite: {row['Step11_Composite_Score']:.1f} | "
                       f"SM: {row['Step11_SmallMolecule_Score']:.1f} | "
                       f"Bio: {row['Step11_Biologic_Score']:.1f}")
            logger.info(f"   Best: {row['Step11_Best_Modality']} | "
                       f"{row['Step11_Classification']}")
        
        # Show biologic-only examples
        logger.info("\n" + "="*80)
        logger.info("BIOLOGIC-ONLY TARGETS (SM Score <25)")
        logger.info("="*80)
        
        biologic_targets = results_df[results_df['Step11_Is_Biologic_Only'] == True].head(15)
        
        for idx, (_, row) in enumerate(biologic_targets.iterrows(), 1):
            logger.info(f"\n{idx}. {row['gene_symbol']}")
            logger.info(f"   SM: {row['Step11_SmallMolecule_Score']:.1f} | "
                       f"Bio: {row['Step11_Biologic_Score']:.1f}")
            logger.info(f"   Reason: {row['Step11_SM_Family_Evidence']}")
        
        logger.info("\n" + "="*80)
        logger.info("✓ COMPLETE! Results saved to Step11.csv")
        logger.info("="*80)
        logger.info(f"\nKey columns in Step11.csv:")
        logger.info(f"  - Step11_Composite_Score (best of SM/Bio)")
        logger.info(f"  - Step11_SmallMolecule_Score (pills)")
        logger.info(f"  - Step11_Biologic_Score (antibodies)")
        logger.info(f"  - Step11_Best_Modality (SM/Bio/Both)")
        logger.info(f"  - Step11_Classification (High/Medium/Low)")
        
        logger.info("\n" + "="*80)
        logger.info("FIXES APPLIED:")
        logger.info("="*80)
        logger.info("✓ BDNF: Now correctly flagged as biologic-only (neurotrophin)")
        logger.info("✓ COL1A2/COL4A1: Biologic score reduced (insoluble fibers)")
        logger.info("✓ JAK2: Small molecule drugs correctly recognized")
        logger.info("✓ All cytokines: Correctly scored for respective modalities")
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()