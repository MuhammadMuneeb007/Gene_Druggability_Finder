"""
Druggability Assessment Pipeline - Step 3: Comprehensive Protein Family Classification
======================================================================================

This script includes ALL major druggable protein families:
- Classic druggable families (GPCRs, Kinases, Ion Channels, etc.)
- Emerging therapeutic targets (E3 ligases, Epigenetic regulators, etc.)
- Immunotherapy targets (Immune checkpoints)
- Structural proteins (Tubulin, Actin)

OUTPUT: Step3.csv with all Step2 fields + comprehensive protein family data
"""

import requests
import pandas as pd
import time
from typing import List, Dict, Optional, Set
from pathlib import Path
import json
import re


class ComprehensiveProteinFamilyClassifier:
    """Classify proteins into comprehensive druggable families including emerging targets."""
    
    def __init__(self):
        self.interpro_url = "https://www.ebi.ac.uk/interpro/api"
        self.uniprot_url = "https://rest.uniprot.org/uniprotkb"
        
        # Comprehensive druggable family definitions
        self.druggable_families = {
            # ===== CLASSIC DRUGGABLE FAMILIES =====
            'GPCR': {
                'keywords': ['g-protein coupled receptor', 'gpcr', '7tm', 'seven transmembrane',
                           'rhodopsin', 'secretin receptor', 'metabotropic glutamate',
                           'frizzled', 'smoothened', 'olfactory receptor'],
                'interpro': ['IPR000276', 'IPR017452', 'IPR000832', 'IPR001879'],
                'pfam': ['PF00001', 'PF00002', 'PF00003'],
                'category': 'GPCR'
            },
            'Kinase': {
                'keywords': ['kinase', 'phosphotransferase', 'protein kinase', 'tyrosine kinase',
                           'serine/threonine kinase', 'receptor tyrosine kinase', 'mapk', 'cdk'],
                'interpro': ['IPR000719', 'IPR011009', 'IPR017441', 'IPR008271', 'IPR020635'],
                'pfam': ['PF00069', 'PF07714'],
                'category': 'Kinase'
            },
            'Ion Channel': {
                'keywords': ['ion channel', 'voltage-gated', 'ligand-gated', 'potassium channel',
                           'sodium channel', 'calcium channel', 'chloride channel',
                           'transient receptor potential', 'trp channel', 'aquaporin',
                           'ionotropic receptor', 'nicotinic receptor'],
                'interpro': ['IPR005821', 'IPR006029', 'IPR006078', 'IPR001622', 'IPR002394'],
                'pfam': ['PF00520', 'PF00027', 'PF00060'],
                'category': 'Ion Channel'
            },
            'Nuclear Receptor': {
                'keywords': ['nuclear receptor', 'nuclear hormone receptor', 'steroid receptor',
                           'thyroid hormone receptor', 'retinoic acid receptor', 'peroxisome proliferator',
                           'estrogen receptor', 'androgen receptor', 'glucocorticoid receptor'],
                'interpro': ['IPR001628', 'IPR008946', 'IPR000536'],
                'pfam': ['PF00104', 'PF00105'],
                'category': 'Nuclear Receptor'
            },
            'Protease': {
                'keywords': ['protease', 'peptidase', 'proteinase', 'metalloprotease',
                           'serine protease', 'cysteine protease', 'aspartic protease',
                           'threonine protease', 'matrix metalloproteinase', 'mmp',
                           'cathepsin', 'caspase', 'thrombin', 'factor xa'],
                'interpro': ['IPR001254', 'IPR001969', 'IPR001563', 'IPR015500', 'IPR029058'],
                'pfam': ['PF00089', 'PF00413', 'PF00082', 'PF00112'],
                'category': 'Protease'
            },
            'Phosphatase': {
                'keywords': ['phosphatase', 'protein phosphatase', 'tyrosine phosphatase',
                           'serine/threonine phosphatase', 'dual specificity phosphatase',
                           'ptp', 'protein tyrosine phosphatase'],
                'interpro': ['IPR000340', 'IPR016130', 'IPR000387', 'IPR029021'],
                'pfam': ['PF00102', 'PF00481'],
                'category': 'Phosphatase'
            },
            'Transporter': {
                'keywords': ['transporter', 'solute carrier', 'abc transporter',
                           'sodium-dependent', 'glucose transporter', 'amino acid transporter',
                           'neurotransmitter transporter', 'drug transporter'],
                'interpro': ['IPR011527', 'IPR001757', 'IPR003439', 'IPR001066'],
                'pfam': ['PF00005', 'PF00664'],
                'category': 'Transporter'
            },
            
            # ===== EMERGING THERAPEUTIC TARGETS =====
            'E3 Ubiquitin Ligase': {
                'keywords': ['e3 ubiquitin ligase', 'ubiquitin ligase', 'ring finger',
                           'hect domain', 'cullin', 'f-box', 'scf complex',
                           'proteasome', 'ubiquitination', 'protein degradation'],
                'interpro': ['IPR001841', 'IPR013083', 'IPR003613', 'IPR001810', 'IPR019474'],
                'pfam': ['PF00097', 'PF12678', 'PF00646', 'PF13639'],
                'category': 'E3 Ubiquitin Ligase'
            },
            'Deubiquitinase': {
                'keywords': ['deubiquitinase', 'deubiquitinating enzyme', 'dub',
                           'ubiquitin-specific protease', 'usp', 'ubiquitin thioesterase',
                           'ubiquitin carboxyl-terminal hydrolase'],
                'interpro': ['IPR001394', 'IPR028889', 'IPR018200'],
                'pfam': ['PF00443', 'PF02338'],
                'category': 'Deubiquitinase'
            },
            'Histone Methyltransferase': {
                'keywords': ['histone methyltransferase', 'histone-lysine methyltransferase',
                           'set domain', 'hmt', 'lysine methyltransferase', 'dot1'],
                'interpro': ['IPR001214', 'IPR003616', 'IPR013128'],
                'pfam': ['PF00856'],
                'category': 'Histone Methyltransferase'
            },
            'Histone Demethylase': {
                'keywords': ['histone demethylase', 'lysine demethylase', 'kdm',
                           'jmjc domain', 'jumonji'],
                'interpro': ['IPR003347', 'IPR013129'],
                'pfam': ['PF02373'],
                'category': 'Histone Demethylase'
            },
            'Histone Acetyltransferase': {
                'keywords': ['histone acetyltransferase', 'hat', 'acetyltransferase',
                           'lysine acetyltransferase', 'kat', 'gcn5', 'pcaf'],
                'interpro': ['IPR000182', 'IPR016181'],
                'pfam': ['PF00583'],
                'category': 'Histone Acetyltransferase'
            },
            'Histone Deacetylase': {
                'keywords': ['histone deacetylase', 'hdac', 'deacetylase',
                           'sirtuin', 'class i hdac', 'class ii hdac'],
                'interpro': ['IPR000286', 'IPR003000'],
                'pfam': ['PF00850'],
                'category': 'Histone Deacetylase'
            },
            'Bromodomain': {
                'keywords': ['bromodomain', 'bet protein', 'brd2', 'brd3', 'brd4',
                           'acetyl-lysine binding', 'chromatin reader'],
                'interpro': ['IPR001487'],
                'pfam': ['PF00439'],
                'category': 'Bromodomain Protein'
            },
            'Chromodomain': {
                'keywords': ['chromodomain', 'heterochromatin protein', 'hp1',
                           'methyl-lysine binding', 'chromatin reader'],
                'interpro': ['IPR000953'],
                'pfam': ['PF00385'],
                'category': 'Chromodomain Protein'
            },
            'DNA Methyltransferase': {
                'keywords': ['dna methyltransferase', 'dnmt', 'cytosine methyltransferase',
                           'dna methylation'],
                'interpro': ['IPR001525', 'IPR018117'],
                'pfam': ['PF00145'],
                'category': 'DNA Methyltransferase'
            },
            
            # ===== STRUCTURAL/CYTOSKELETAL PROTEINS =====
            'Tubulin': {
                'keywords': ['tubulin', 'alpha-tubulin', 'beta-tubulin', 'microtubule',
                           'gamma-tubulin'],
                'interpro': ['IPR000217', 'IPR003008'],
                'pfam': ['PF00091', 'PF03953'],
                'category': 'Tubulin'
            },
            'Actin': {
                'keywords': ['actin', 'alpha-actin', 'beta-actin', 'gamma-actin',
                           'cytoskeletal protein', 'microfilament'],
                'interpro': ['IPR004000'],
                'pfam': ['PF00022'],
                'category': 'Actin'
            },
            'Motor Protein': {
                'keywords': ['kinesin', 'dynein', 'myosin', 'motor protein',
                           'microtubule motor', 'molecular motor'],
                'interpro': ['IPR001752', 'IPR027417', 'IPR001609'],
                'pfam': ['PF00225', 'PF00063'],
                'category': 'Motor Protein'
            },
            
            # ===== IMMUNE CHECKPOINT & IMMUNOTHERAPY TARGETS =====
            'Immune Checkpoint': {
                'keywords': ['pd-1', 'pd-l1', 'ctla-4', 'cd274', 'pdcd1', 'lag-3',
                           'tim-3', 'b7', 'immune checkpoint', 'programmed death',
                           'cytotoxic t-lymphocyte'],
                'interpro': ['IPR003599', 'IPR013783'],
                'pfam': ['PF00047'],
                'category': 'Immune Checkpoint'
            },
            'Cytokine Receptor': {
                'keywords': ['cytokine receptor', 'interleukin receptor', 'interferon receptor',
                           'tnf receptor', 'chemokine receptor', 'colony stimulating factor receptor'],
                'interpro': ['IPR003531', 'IPR001859', 'IPR003033'],
                'pfam': ['PF00019', 'PF01030'],
                'category': 'Cytokine Receptor'
            },
            'T-cell Receptor': {
                'keywords': ['t-cell receptor', 'tcr', 't cell antigen receptor',
                           'cd3', 'cd4', 'cd8'],
                'interpro': ['IPR003599', 'IPR007110'],
                'pfam': ['PF00047'],
                'category': 'T-cell Receptor'
            },
            
            # ===== OTHER IMPORTANT TARGETS =====
            'Chaperone': {
                'keywords': ['heat shock protein', 'hsp90', 'hsp70', 'hsp60',
                           'chaperone', 'chaperonin', 'protein folding'],
                'interpro': ['IPR001404', 'IPR029047', 'IPR013126'],
                'pfam': ['PF00012', 'PF00183'],
                'category': 'Chaperone'
            },
            'Cyclooxygenase': {
                'keywords': ['cyclooxygenase', 'cox-1', 'cox-2', 'prostaglandin',
                           'prostaglandin-endoperoxide synthase', 'ptgs'],
                'interpro': ['IPR000898'],
                'pfam': ['PF00008'],
                'category': 'Cyclooxygenase'
            },
            'Lipid Kinase': {
                'keywords': ['phosphatidylinositol 3-kinase', 'pi3k', 'pi3-kinase',
                           'lipid kinase', 'phosphoinositide kinase'],
                'interpro': ['IPR000403', 'IPR015433'],
                'pfam': ['PF00454'],
                'category': 'Lipid Kinase'
            },
            'GTPase': {
                'keywords': ['gtpase', 'ras protein', 'ras gtpase', 'rho gtpase',
                           'rab protein', 'ran protein', 'small gtpase'],
                'interpro': ['IPR001806', 'IPR005225', 'IPR027417'],
                'pfam': ['PF00071'],
                'category': 'GTPase'
            },
            'Integrin': {
                'keywords': ['integrin', 'cell adhesion', 'extracellular matrix receptor',
                           'integrin alpha', 'integrin beta'],
                'interpro': ['IPR013517', 'IPR002369'],
                'pfam': ['PF00357', 'PF08441'],
                'category': 'Integrin'
            },
            'Metalloenzyme': {
                'keywords': ['metalloenzyme', 'zinc finger', 'zinc-binding',
                           'carbonic anhydrase', 'alcohol dehydrogenase'],
                'interpro': ['IPR001353', 'IPR013149', 'IPR013154'],
                'pfam': ['PF00107', 'PF00096'],
                'category': 'Metalloenzyme'
            },
            'Transcription Factor': {
                'keywords': ['transcription factor', 'dna-binding', 'helix-turn-helix',
                           'zinc finger transcription', 'homeobox', 'leucine zipper'],
                'interpro': ['IPR001356', 'IPR009057', 'IPR017970'],
                'pfam': ['PF00010', 'PF00046'],
                'category': 'Transcription Factor'
            },
            
            # ===== GENERAL ENZYME CLASSES =====
            'Oxidoreductase': {
                'keywords': ['oxidoreductase', 'dehydrogenase', 'oxidase', 'reductase',
                           'peroxidase', 'monooxygenase', 'dioxygenase'],
                'interpro': ['IPR016040', 'IPR001327'],
                'pfam': [],
                'category': 'Oxidoreductase'
            },
            'Transferase': {
                'keywords': ['transferase', 'methyltransferase', 'glycosyltransferase',
                           'acyltransferase', 'aminotransferase'],
                'interpro': ['IPR016461'],
                'pfam': [],
                'category': 'Transferase'
            },
            'Hydrolase': {
                'keywords': ['hydrolase', 'esterase', 'lipase', 'amidase',
                           'glycosidase', 'nuclease'],
                'interpro': [],
                'pfam': [],
                'category': 'Hydrolase'
            },
            'Lyase': {
                'keywords': ['lyase', 'decarboxylase', 'aldolase', 'synthase'],
                'interpro': [],
                'pfam': [],
                'category': 'Lyase'
            },
            'Isomerase': {
                'keywords': ['isomerase', 'racemase', 'epimerase', 'mutase'],
                'interpro': [],
                'pfam': [],
                'category': 'Isomerase'
            },
            'Ligase': {
                'keywords': ['ligase', 'synthetase', 'carboxylase', 'dna ligase'],
                'interpro': [],
                'pfam': [],
                'category': 'Ligase'
            }
        }
    
    def query_interpro_protein(self, uniprot_id: str) -> Dict:
        """Query InterPro for protein domain and family information."""
        result = {
            'found': False,
            'interpro_entries': [],
            'entry_types': [],
            'domains': [],
            'families': [],
            'signatures': []
        }
        
        if not uniprot_id or pd.isna(uniprot_id):
            return result
        
        url = f"{self.interpro_url}/entry/interpro/protein/UniProt/{uniprot_id}"
        
        response = requests.get(url, timeout=15)
        
        if response.status_code == 200:
            data = response.json()
            
            if 'results' in data and len(data['results']) > 0:
                result['found'] = True
                
                for entry in data['results']:
                    metadata = entry.get('metadata', {})
                    entry_id = metadata.get('accession')
                    
                    # Handle name field which can be string or dict
                    entry_name = ''
                    name_field = metadata.get('name', '')
                    if isinstance(name_field, dict):
                        entry_name = name_field.get('name', '')
                    elif isinstance(name_field, str):
                        entry_name = name_field
                    
                    entry_type = metadata.get('type', '')
                    
                    if entry_id:
                        result['interpro_entries'].append(entry_id)
                    
                    if entry_type:
                        result['entry_types'].append(entry_type)
                        
                        if entry_type == 'Domain' and entry_name:
                            result['domains'].append(f"{entry_id}:{entry_name}")
                        elif entry_type == 'Family' and entry_name:
                            result['families'].append(f"{entry_id}:{entry_name}")
        
        return result
    
    def query_uniprot_features(self, uniprot_id: str) -> Dict:
        """Query UniProt for detailed protein features and annotations."""
        result = {
            'found': False,
            'protein_name': None,
            'gene_name': None,
            'protein_families': [],
            'protein_description': None,
            'keywords': [],
            'pfam_domains': [],
            'interpro_domains': [],
            'function': None
        }
        
        if not uniprot_id or pd.isna(uniprot_id):
            return result
        
        url = f"{self.uniprot_url}/{uniprot_id}"
        
        params = {
            'format': 'json',
            'fields': 'protein_name,gene_names,cc_function,keyword,xref_pfam,xref_interpro,ft_domain,protein_families'
        }
        
        response = requests.get(url, params=params, timeout=15)
        
        if response.status_code == 200:
            data = response.json()
            result['found'] = True
            
            # Protein name
            if 'proteinDescription' in data:
                prot_desc = data['proteinDescription']
                if 'recommendedName' in prot_desc:
                    rec_name = prot_desc['recommendedName'].get('fullName', {})
                    if isinstance(rec_name, dict):
                        result['protein_name'] = rec_name.get('value')
                    else:
                        result['protein_name'] = rec_name
            
            # Gene name
            if 'genes' in data and len(data['genes']) > 0:
                gene = data['genes'][0]
                if 'geneName' in gene:
                    result['gene_name'] = gene['geneName'].get('value')
            
            # Protein families from comments
            if 'comments' in data:
                for comment in data['comments']:
                    if comment.get('commentType') == 'SIMILARITY':
                        for text in comment.get('texts', []):
                            family_text = text.get('value', '')
                            if family_text:
                                result['protein_families'].append(family_text)
                    elif comment.get('commentType') == 'FUNCTION':
                        for text in comment.get('texts', []):
                            result['function'] = text.get('value', '')
            
            # Keywords
            if 'keywords' in data:
                for kw in data['keywords']:
                    result['keywords'].append(kw.get('name', ''))
            
            # Pfam domains from cross-references
            if 'uniProtKBCrossReferences' in data:
                for xref in data['uniProtKBCrossReferences']:
                    if xref.get('database') == 'Pfam':
                        pfam_id = xref.get('id')
                        pfam_name = xref.get('properties', [{}])[0].get('value', '') if xref.get('properties') else ''
                        if pfam_id:
                            result['pfam_domains'].append(f"{pfam_id}:{pfam_name}" if pfam_name else pfam_id)
                    elif xref.get('database') == 'InterPro':
                        interpro_id = xref.get('id')
                        interpro_name = xref.get('properties', [{}])[0].get('value', '') if xref.get('properties') else ''
                        if interpro_id:
                            result['interpro_domains'].append(f"{interpro_id}:{interpro_name}" if interpro_name else interpro_id)
        
        return result
    
    def classify_druggable_family(self, protein_data: Dict, interpro_data: Dict) -> Dict:
        """Classify protein into comprehensive druggable families."""
        result = {
            'is_druggable_family': False,
            'druggable_family_types': [],
            'confidence': 'NONE',
            'evidence': []
        }
        
        # Combine all text for searching
        search_text = ' '.join([
            protein_data.get('protein_name', '') or '',
            protein_data.get('function', '') or '',
            ' '.join(protein_data.get('protein_families', [])),
            ' '.join(protein_data.get('keywords', [])),
        ]).lower()
        
        # Get all domain/family IDs
        interpro_ids = set(interpro_data.get('interpro_entries', []))
        pfam_ids = set()
        for domain in protein_data.get('pfam_domains', []):
            pfam_id = domain.split(':')[0] if ':' in domain else domain
            pfam_ids.add(pfam_id)
        
        matched_families = []
        evidence_list = []
        
        # Check each druggable family
        for family_name, family_info in self.druggable_families.items():
            matched = False
            match_evidence = []
            
            # Check keywords
            for keyword in family_info['keywords']:
                if keyword.lower() in search_text:
                    matched = True
                    match_evidence.append(f"Keyword: {keyword}")
            
            # Check InterPro IDs
            for interpro_id in family_info.get('interpro', []):
                if interpro_id in interpro_ids:
                    matched = True
                    match_evidence.append(f"InterPro: {interpro_id}")
            
            # Check Pfam IDs
            for pfam_id in family_info.get('pfam', []):
                if pfam_id in pfam_ids:
                    matched = True
                    match_evidence.append(f"Pfam: {pfam_id}")
            
            if matched:
                matched_families.append(family_info['category'])
                evidence_list.extend([f"{family_info['category']}: {ev}" for ev in match_evidence])
        
        # Set results
        if matched_families:
            result['is_druggable_family'] = True
            result['druggable_family_types'] = list(set(matched_families))
            
            # Determine confidence based on number of matches
            if len(evidence_list) >= 3:
                result['confidence'] = 'HIGH'
            elif len(evidence_list) >= 2:
                result['confidence'] = 'MEDIUM'
            else:
                result['confidence'] = 'LOW'
            
            result['evidence'] = evidence_list
        
        return result
    
    def process_gene_row(self, row: pd.Series) -> Dict:
        """Process a single gene and classify its protein family."""
        gene_symbol = row['gene_symbol']
        uniprot_id = row.get('uniprot_swissprot')
        
        print(f"  {gene_symbol}...", end=' ')
        
        result = {
            'gene_symbol': gene_symbol,
            'uniprot_protein_name': None,
            'uniprot_function': None,
            'uniprot_keywords': None,
            'uniprot_protein_families': None,
            'pfam_domains': None,
            'pfam_domain_count': 0,
            'interpro_domains': None,
            'interpro_domain_count': 0,
            'interpro_families': None,
            'interpro_entry_types': None,
            'is_druggable_family': False,
            'druggable_family_type': None,
            'druggable_family_confidence': 'NONE',
            'druggable_family_evidence': None,
            'protein_class': None
        }
        
        if not uniprot_id or pd.isna(uniprot_id):
            print("✗ No UniProt ID")
            return result
        
        # Query UniProt
        protein_data = self.query_uniprot_features(uniprot_id)
        
        if protein_data['found']:
            result['uniprot_protein_name'] = protein_data.get('protein_name')
            result['uniprot_function'] = protein_data.get('function')
            
            keywords = protein_data.get('keywords', [])
            if keywords:
                result['uniprot_keywords'] = '|'.join(keywords[:20])
            
            families = protein_data.get('protein_families', [])
            if families:
                result['uniprot_protein_families'] = '|'.join(families)
            
            pfam_domains = protein_data.get('pfam_domains', [])
            if pfam_domains:
                result['pfam_domains'] = '|'.join(pfam_domains[:15])
                result['pfam_domain_count'] = len(pfam_domains)
            
            interpro_domains = protein_data.get('interpro_domains', [])
            if interpro_domains:
                result['interpro_domains'] = '|'.join(interpro_domains[:15])
                result['interpro_domain_count'] = len(interpro_domains)
            
            print(f"Pfam:{len(pfam_domains)} InterPro:{len(interpro_domains)}", end=' ')
        else:
            print("UniProt:✗", end=' ')
        
        time.sleep(0.3)
        
        # Query InterPro
        interpro_data = self.query_interpro_protein(uniprot_id)
        
        if interpro_data['found']:
            families = interpro_data.get('families', [])
            if families:
                result['interpro_families'] = '|'.join(families[:10])
            
            entry_types = interpro_data.get('entry_types', [])
            if entry_types:
                result['interpro_entry_types'] = '|'.join(set(entry_types))
        
        time.sleep(0.3)
        
        # Classify into druggable families
        classification = self.classify_druggable_family(protein_data, interpro_data)
        
        result['is_druggable_family'] = classification['is_druggable_family']
        result['druggable_family_confidence'] = classification['confidence']
        
        if classification['druggable_family_types']:
            result['druggable_family_type'] = '|'.join(classification['druggable_family_types'])
        
        if classification['evidence']:
            result['druggable_family_evidence'] = '|'.join(classification['evidence'][:10])
        
        # Determine protein class
        if classification['is_druggable_family']:
            result['protein_class'] = classification['druggable_family_types'][0]
        elif keywords:
            for kw in keywords[:5]:
                if any(term in kw.lower() for term in ['receptor', 'enzyme', 'channel', 'transporter']):
                    result['protein_class'] = kw
                    break
        
        # Print result
        if result['is_druggable_family']:
            print(f"✓ DRUGGABLE: {result['druggable_family_type']} ({result['druggable_family_confidence']})")
        else:
            print("✗ Not in druggable family")
        
        return result
    
    def process_step2_file(self, input_file: str = "Step2.csv") -> pd.DataFrame:
        """Process Step2.csv and add comprehensive protein family classification."""
        print("="*70)
        print("COMPREHENSIVE PROTEIN FAMILY CLASSIFICATION - STEP 3")
        print("="*70)
        print("\nIncluding:")
        print("  • Classic targets: GPCRs, Kinases, Ion Channels, Nuclear Receptors")
        print("  • Proteases & Phosphatases")
        print("  • Epigenetic regulators: HMTs, HDMs, HATs, HDACs, Bromodomains")
        print("  • E3 Ligases & Deubiquitinases (PROTAC targets)")
        print("  • Immune checkpoints & Cytokine receptors")
        print("  • Structural proteins: Tubulin, Actin, Motor proteins")
        print("  • And more...")
        
        # Read Step2.csv
        print(f"\nReading {input_file}...")
        step2_df = pd.read_csv(input_file)
        print(f"  ✓ Loaded {len(step2_df)} genes")
        
        # Process each gene
        print(f"\nClassifying protein families for {len(step2_df)} genes...\n")
        
        family_results = []
        
        for idx, row in step2_df.iterrows():
            print(f"[{idx+1}/{len(step2_df)}]", end=' ')
            result = self.process_gene_row(row)
            family_results.append(result)
        
        # Create family data DataFrame
        family_df = pd.DataFrame(family_results)
        
        # Merge with Step2 data
        family_df = family_df.drop(columns=['gene_symbol'])
        step3_df = pd.concat([step2_df, family_df], axis=1)
        
        # Print summary
        self.print_summary(step3_df)
        
        return step3_df
    
    def print_summary(self, df: pd.DataFrame):
        """Print comprehensive summary of protein family classification."""
        print("\n" + "="*70)
        print("COMPREHENSIVE PROTEIN FAMILY CLASSIFICATION SUMMARY")
        print("="*70)
        
        total_genes = len(df)
        druggable = (df['is_druggable_family'] == True).sum()
        not_druggable = total_genes - druggable
        
        print(f"\nTotal genes analyzed: {total_genes}")
        print(f"\nDruggable Family Classification:")
        print(f"  ✓ In druggable family:     {druggable} ({druggable/total_genes*100:.1f}%)")
        print(f"  ✗ Not in druggable family: {not_druggable} ({not_druggable/total_genes*100:.1f}%)")
        
        # Confidence levels
        if druggable > 0:
            high_conf = (df['druggable_family_confidence'] == 'HIGH').sum()
            medium_conf = (df['druggable_family_confidence'] == 'MEDIUM').sum()
            low_conf = (df['druggable_family_confidence'] == 'LOW').sum()
            
            print(f"\nConfidence Levels (for druggable genes):")
            print(f"  HIGH:   {high_conf} ({high_conf/druggable*100:.1f}%)")
            print(f"  MEDIUM: {medium_conf} ({medium_conf/druggable*100:.1f}%)")
            print(f"  LOW:    {low_conf} ({low_conf/druggable*100:.1f}%)")
        
        # Count by family type
        print(f"\nDruggable Family Types Found:")
        family_counts = {}
        
        for idx, row in df[df['is_druggable_family'] == True].iterrows():
            family_types = row.get('druggable_family_type', '')
            if family_types and not pd.isna(family_types):
                for family in str(family_types).split('|'):
                    family = family.strip()
                    family_counts[family] = family_counts.get(family, 0) + 1
        
        # Group by category
        classic = {}
        emerging = {}
        immune = {}
        structural = {}
        other = {}
        
        classic_families = ['GPCR', 'Kinase', 'Ion Channel', 'Nuclear Receptor', 'Protease', 'Phosphatase', 'Transporter']
        emerging_families = ['E3 Ubiquitin Ligase', 'Deubiquitinase', 'Histone Methyltransferase', 'Histone Demethylase', 
                           'Histone Acetyltransferase', 'Histone Deacetylase', 'Bromodomain Protein', 'Chromodomain Protein', 'DNA Methyltransferase']
        immune_families = ['Immune Checkpoint', 'Cytokine Receptor', 'T-cell Receptor']
        structural_families = ['Tubulin', 'Actin', 'Motor Protein']
        
        for family, count in family_counts.items():
            if family in classic_families:
                classic[family] = count
            elif family in emerging_families:
                emerging[family] = count
            elif family in immune_families:
                immune[family] = count
            elif family in structural_families:
                structural[family] = count
            else:
                other[family] = count
        
        if classic:
            print("\n  🎯 Classic Druggable Families:")
            for family, count in sorted(classic.items(), key=lambda x: x[1], reverse=True):
                print(f"     {family:30s}: {count:3d} genes")
        
        if emerging:
            print("\n  🧬 Emerging Therapeutic Targets:")
            for family, count in sorted(emerging.items(), key=lambda x: x[1], reverse=True):
                print(f"     {family:30s}: {count:3d} genes")
        
        if immune:
            print("\n  🛡️  Immunotherapy Targets:")
            for family, count in sorted(immune.items(), key=lambda x: x[1], reverse=True):
                print(f"     {family:30s}: {count:3d} genes")
        
        if structural:
            print("\n  🏗️  Structural/Cytoskeletal Proteins:")
            for family, count in sorted(structural.items(), key=lambda x: x[1], reverse=True):
                print(f"     {family:30s}: {count:3d} genes")
        
        if other:
            print("\n  ⚗️  Other Targets:")
            for family, count in sorted(other.items(), key=lambda x: x[1], reverse=True):
                print(f"     {family:30s}: {count:3d} genes")
        
        # Domain statistics
        print(f"\n📊 Domain Coverage:")
        print(f"  Genes with Pfam domains:     {(df['pfam_domain_count'] > 0).sum()} ({(df['pfam_domain_count'] > 0).sum()/total_genes*100:.1f}%)")
        print(f"  Genes with InterPro domains: {(df['interpro_domain_count'] > 0).sum()} ({(df['interpro_domain_count'] > 0).sum()/total_genes*100:.1f}%)")
        print(f"  Average Pfam domains/gene:   {df['pfam_domain_count'].mean():.1f}")
        print(f"  Average InterPro domains/gene: {df['interpro_domain_count'].mean():.1f}")
        
        # Top druggable genes
        print(f"\n🎯 Top 10 Druggable Genes:")
        druggable_genes = df[df['is_druggable_family'] == True].sort_values(
            'druggable_family_confidence', 
            ascending=False
        ).head(10)
        
        if len(druggable_genes) > 0:
            for idx, row in druggable_genes.iterrows():
                family_type = row.get('druggable_family_type', 'Unknown')
                confidence = row.get('druggable_family_confidence', 'NONE')
                print(f"  • {row['gene_symbol']:12s} - {family_type:35s} ({confidence})")
        else:
            print("  No genes in druggable families found.")
        
        print("="*70)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    # Check if Step2.csv exists
    if not Path("Step2.csv").exists():
        print("❌ Error: Step2.csv not found!")
        print("Please run Step 2 first to generate Step2.csv")
        exit(1)
    
    # Create comprehensive protein family classifier
    classifier = ComprehensiveProteinFamilyClassifier()
    
    # Process Step2.csv
    step3_df = classifier.process_step2_file("Step2.csv")
    
    # Save to Step3.csv
    output_path = Path("Step3.csv")
    step3_df.to_csv(output_path, index=False)
    
    print(f"\n✅ Results saved to: {output_path.absolute()}")
    print(f"\nTotal columns: {len(step3_df.columns)}")
    print(f"Total rows: {len(step3_df)}")
    
    # Display sample
    print("\n📋 Sample Results (first 5 genes):")
    sample_cols = ['gene_symbol', 'is_druggable_family', 'druggable_family_type',
                   'druggable_family_confidence']
    if all(col in step3_df.columns for col in sample_cols):
        print(step3_df[sample_cols].head(5).to_string(index=False))
    
    print("\n✅ Step 3 Complete!")
    print("\n📚 Comprehensive druggable family classification including:")
    print("   • All classic druggable families")
    print("   • Epigenetic regulators (HMTs, HDACs, Bromodomains, etc.)")
    print("   • E3 Ligases & DUBs (PROTAC targets)")
    print("   • Immune checkpoints & immunotherapy targets")
    print("   • Structural proteins (Tubulin, Actin)")
    print("   • And many more therapeutic target classes")
    print("\nNext: Proceed to Step 4 (Protein Structure Analysis)")