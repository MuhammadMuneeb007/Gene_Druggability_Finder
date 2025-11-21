"""
Druggability Assessment Pipeline - Step 2: Drug-Gene Interactions
==================================================================

This script takes Step1.csv and queries drug databases:
- DGIdb (Drug-Gene Interaction database) using GraphQL
- ChEMBL (bioactive compounds)
- Open Targets (known drugs)

OUTPUT: Step2.csv with all Step1 fields + drug interaction data
"""

import requests
import pandas as pd
import time
from typing import List, Dict, Optional
from pathlib import Path
import json


class DrugGeneInteractionFinder:
    """Query multiple drug databases for gene-drug interactions."""
    
    def __init__(self):
        self.dgidb_url = "https://dgidb.org/api/graphql"
        self.dgidb_headers = {"Content-Type": "application/json"}
        self.chembl_url = "https://www.ebi.ac.uk/chembl/api/data"
        self.opentargets_url = "https://api.platform.opentargets.org/api/v4/graphql"
        
    def query_dgidb(self, gene_symbol: str) -> Dict:
        """
        Query DGIdb for drug-gene interactions using GraphQL.
        
        Args:
            gene_symbol: Gene symbol
            
        Returns:
            Dictionary with interaction data
        """
        query = """
        query {
          geneMatches(searchTerms: ["%s"]) {
            directMatches {
              matches {
                name
                longName
                geneCategories {
                  name
                }
                interactions {
                  drug {
                    name
                    approved
                  }
                  interactionScore
                  interactionTypes {
                    type
                  }
                  publications {
                    pmid
                  }
                  sources {
                    fullName
                  }
                }
              }
            }
          }
        }
        """ % gene_symbol
        
        response = requests.post(
            self.dgidb_url, 
            json={"query": query}, 
            headers=self.dgidb_headers, 
            timeout=15
        )
        data = response.json()
        
        result = {
            'gene': gene_symbol,
            'found': False,
            'long_name': None,
            'categories': [],
            'num_interactions': 0,
            'approved_drugs': [],
            'all_drugs': [],
            'top_5_drugs': [],
            'interaction_types': [],
            'sources': []
        }
        
        if 'data' in data and data['data']['geneMatches']['directMatches']:
            matches = data['data']['geneMatches']['directMatches']
            if matches and matches[0]['matches']:
                gene_data = matches[0]['matches'][0]
                interactions = gene_data['interactions']
                
                result['found'] = True
                result['long_name'] = gene_data.get('longName')
                result['num_interactions'] = len(interactions)
                
                # Get categories
                categories = [cat['name'] for cat in gene_data.get('geneCategories', [])]
                result['categories'] = categories
                
                # Get approved drugs
                approved_drugs = [i['drug']['name'] for i in interactions if i['drug']['approved']]
                result['approved_drugs'] = approved_drugs
                
                # Get all drugs
                all_drugs = [i['drug']['name'] for i in interactions]
                result['all_drugs'] = all_drugs
                
                # Get top 5 by interaction score
                sorted_interactions = sorted(
                    interactions, 
                    key=lambda x: x['interactionScore'], 
                    reverse=True
                )
                top_drugs = [i['drug']['name'] for i in sorted_interactions[:5]]
                result['top_5_drugs'] = top_drugs
                
                # Get interaction types
                interaction_types = set()
                for i in interactions:
                    for int_type in i.get('interactionTypes', []):
                        if int_type.get('type'):
                            interaction_types.add(int_type['type'])
                result['interaction_types'] = list(interaction_types)
                
                # Get sources
                sources = set()
                for i in interactions:
                    for source in i.get('sources', []):
                        if source.get('fullName'):
                            sources.add(source['fullName'])
                result['sources'] = list(sources)
        
        return result
    
    def query_chembl_targets(self, gene_symbol: str) -> Dict:
        """
        Query ChEMBL for bioactive compounds targeting the gene.
        
        Args:
            gene_symbol: Gene symbol
            
        Returns:
            Dictionary with ChEMBL data
        """
        # First, search for target by gene name
        search_url = f"{self.chembl_url}/target/search.json"
        
        params = {
            'q': gene_symbol,
            'format': 'json'
        }
        
        response = requests.get(search_url, params=params, timeout=15)
        data = response.json()
        
        result = {
            'found': False,
            'target_chembl_id': None,
            'target_type': None,
            'num_approved_drugs': 0,
            'num_clinical_compounds': 0,
            'num_bioactive_compounds': 0,
            'total_unique_compounds': 0,
            'total_activities': 0
        }
        
        if 'targets' in data and len(data['targets']) > 0:
            # Get the first matching target
            target = data['targets'][0]
            target_chembl_id = target.get('target_chembl_id')
            
            if target_chembl_id:
                result['found'] = True
                result['target_chembl_id'] = target_chembl_id
                result['target_type'] = target.get('target_type')
                
                # Now get activities/compounds for this target
                activity_url = f"{self.chembl_url}/activity.json"
                
                activity_params = {
                    'target_chembl_id': target_chembl_id,
                    'limit': 1000,
                    'format': 'json'
                }
                
                activity_response = requests.get(activity_url, params=activity_params, timeout=20)
                activity_data = activity_response.json()
                
                activities = activity_data.get('activities', [])
                
                # Count compounds by phase
                approved_drugs = set()
                clinical_compounds = set()
                bioactive_compounds = set()
                unique_compounds = set()
                
                for activity in activities:
                    molecule_chembl_id = activity.get('molecule_chembl_id')
                    if molecule_chembl_id:
                        unique_compounds.add(molecule_chembl_id)
                        
                        # Check max phase (4 = approved, 3 = phase 3, etc.)
                        max_phase = activity.get('molecule_max_phase')
                        if max_phase == 4:
                            approved_drugs.add(molecule_chembl_id)
                        elif max_phase in [1, 2, 3]:
                            clinical_compounds.add(molecule_chembl_id)
                        else:
                            bioactive_compounds.add(molecule_chembl_id)
                
                result['num_approved_drugs'] = len(approved_drugs)
                result['num_clinical_compounds'] = len(clinical_compounds)
                result['num_bioactive_compounds'] = len(bioactive_compounds)
                result['total_unique_compounds'] = len(unique_compounds)
                result['total_activities'] = len(activities)
        
        return result
    
    def query_opentargets(self, ensembl_id: str) -> Dict:
        """
        Query Open Targets Platform for known drugs.
        
        Args:
            ensembl_id: Ensembl gene ID
            
        Returns:
            Dictionary with drug data
        """
        result = {
            'found': False,
            'total_known_drugs': 0,
            'unique_drugs': 0,
            'approved_drugs': 0,
            'phase3_drugs': 0,
            'phase2_drugs': 0,
            'phase1_drugs': 0,
            'preclinical': 0,
            'drug_names': [],
            'drug_types': [],
            'mechanisms': []
        }
        
        if not ensembl_id or pd.isna(ensembl_id):
            return result
        
        # GraphQL query for known drugs
        query = """
        query targetKnownDrugs($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            id
            approvedSymbol
            knownDrugs {
              uniqueDrugs
              uniqueTargets
              count
              rows {
                drug {
                  id
                  name
                  maximumClinicalTrialPhase
                  drugType
                }
                drugType
                mechanismOfAction
                phase
                status
              }
            }
          }
        }
        """
        
        variables = {"ensemblId": ensembl_id}
        
        response = requests.post(
            self.opentargets_url,
            json={'query': query, 'variables': variables},
            timeout=20
        )
        data = response.json()
        
        if 'data' in data and data['data'].get('target'):
            target_data = data['data']['target']
            known_drugs = target_data.get('knownDrugs')
            
            if known_drugs:
                result['found'] = True
                result['total_known_drugs'] = known_drugs.get('count', 0)
                result['unique_drugs'] = known_drugs.get('uniqueDrugs', 0)
                
                # Count drugs by phase
                approved = 0
                phase3 = 0
                phase2 = 0
                phase1 = 0
                preclinical = 0
                
                drug_names = []
                drug_types = set()
                mechanisms = set()
                
                for row in known_drugs.get('rows', []):
                    phase = row.get('phase')
                    if phase == 4:
                        approved += 1
                    elif phase == 3:
                        phase3 += 1
                    elif phase == 2:
                        phase2 += 1
                    elif phase == 1:
                        phase1 += 1
                    else:
                        preclinical += 1
                    
                    drug = row.get('drug', {})
                    drug_name = drug.get('name')
                    if drug_name and len(drug_names) < 10:
                        drug_names.append(drug_name)
                    
                    drug_type = row.get('drugType')
                    if drug_type:
                        drug_types.add(drug_type)
                    
                    moa = row.get('mechanismOfAction')
                    if moa:
                        mechanisms.add(moa)
                
                result['approved_drugs'] = approved
                result['phase3_drugs'] = phase3
                result['phase2_drugs'] = phase2
                result['phase1_drugs'] = phase1
                result['preclinical'] = preclinical
                result['drug_names'] = drug_names
                result['drug_types'] = list(drug_types)
                result['mechanisms'] = list(mechanisms)
        
        return result
    
    def process_gene_row(self, row: pd.Series) -> Dict:
        """
        Process a single gene and query all drug databases.
        
        Args:
            row: Row from Step1 DataFrame
            
        Returns:
            Dictionary with drug interaction data
        """
        gene_symbol = row['gene_symbol']
        ensembl_id = row.get('ensembl_gene_id')
        
        print(f"  {gene_symbol}...", end=' ')
        
        result = {
            'gene_symbol': gene_symbol,
            
            # DGIdb results
            'dgidb_found': False,
            'dgidb_long_name': None,
            'dgidb_num_interactions': 0,
            'dgidb_num_approved': 0,
            'dgidb_drugs': None,
            'dgidb_approved_drugs': None,
            'dgidb_top_5_drugs': None,
            'dgidb_categories': None,
            'dgidb_interaction_types': None,
            'dgidb_sources': None,
            
            # ChEMBL results
            'chembl_found': False,
            'chembl_target_id': None,
            'chembl_target_type': None,
            'chembl_approved_drugs': 0,
            'chembl_clinical_compounds': 0,
            'chembl_bioactive_compounds': 0,
            'chembl_total_compounds': 0,
            'chembl_total_activities': 0,
            
            # Open Targets results
            'opentargets_found': False,
            'opentargets_total_drugs': 0,
            'opentargets_approved': 0,
            'opentargets_phase3': 0,
            'opentargets_phase2': 0,
            'opentargets_phase1': 0,
            'opentargets_drug_names': None,
            'opentargets_drug_types': None,
            'opentargets_mechanisms': None,
            
            # Summary flags
            'has_approved_drugs': False,
            'has_clinical_drugs': False,
            'has_bioactive_compounds': False,
            'druggability_evidence': None
        }
        
        # Query DGIdb
        dgidb_data = self.query_dgidb(gene_symbol)
        if dgidb_data['found']:
            result['dgidb_found'] = True
            result['dgidb_long_name'] = dgidb_data.get('long_name')
            result['dgidb_num_interactions'] = dgidb_data.get('num_interactions', 0)
            result['dgidb_num_approved'] = len(dgidb_data.get('approved_drugs', []))
            
            # Store drug lists as pipe-separated strings
            all_drugs = dgidb_data.get('all_drugs', [])
            if all_drugs:
                result['dgidb_drugs'] = '|'.join(all_drugs[:20])  # Limit to 20
            
            approved_drugs = dgidb_data.get('approved_drugs', [])
            if approved_drugs:
                result['dgidb_approved_drugs'] = '|'.join(approved_drugs[:10])
            
            top_5 = dgidb_data.get('top_5_drugs', [])
            if top_5:
                result['dgidb_top_5_drugs'] = '|'.join(top_5)
            
            categories = dgidb_data.get('categories', [])
            if categories:
                result['dgidb_categories'] = '|'.join(categories)
            
            int_types = dgidb_data.get('interaction_types', [])
            if int_types:
                result['dgidb_interaction_types'] = '|'.join(int_types[:10])
            
            sources = dgidb_data.get('sources', [])
            if sources:
                result['dgidb_sources'] = '|'.join(sources[:10])
            
            print(f"DGI:{result['dgidb_num_interactions']}", end=' ')
        else:
            print("DGI:0", end=' ')
        
        time.sleep(0.5)  # Rate limiting
        
        # Query ChEMBL
        chembl_data = self.query_chembl_targets(gene_symbol)
        if chembl_data['found']:
            result['chembl_found'] = True
            result['chembl_target_id'] = chembl_data.get('target_chembl_id')
            result['chembl_target_type'] = chembl_data.get('target_type')
            result['chembl_approved_drugs'] = chembl_data.get('num_approved_drugs', 0)
            result['chembl_clinical_compounds'] = chembl_data.get('num_clinical_compounds', 0)
            result['chembl_bioactive_compounds'] = chembl_data.get('num_bioactive_compounds', 0)
            result['chembl_total_compounds'] = chembl_data.get('total_unique_compounds', 0)
            result['chembl_total_activities'] = chembl_data.get('total_activities', 0)
            
            print(f"ChEMBL:{result['chembl_approved_drugs']}approved", end=' ')
        else:
            print("ChEMBL:0", end=' ')
        
        time.sleep(0.5)  # Rate limiting
        
        # Query Open Targets
        if ensembl_id and not pd.isna(ensembl_id):
            ot_data = self.query_opentargets(ensembl_id)
            if ot_data['found']:
                result['opentargets_found'] = True
                result['opentargets_total_drugs'] = ot_data.get('total_known_drugs', 0)
                result['opentargets_approved'] = ot_data.get('approved_drugs', 0)
                result['opentargets_phase3'] = ot_data.get('phase3_drugs', 0)
                result['opentargets_phase2'] = ot_data.get('phase2_drugs', 0)
                result['opentargets_phase1'] = ot_data.get('phase1_drugs', 0)
                
                drug_names = ot_data.get('drug_names', [])
                if drug_names:
                    result['opentargets_drug_names'] = '|'.join(drug_names)
                
                drug_types = ot_data.get('drug_types', [])
                if drug_types:
                    result['opentargets_drug_types'] = '|'.join(drug_types)
                
                mechanisms = ot_data.get('mechanisms', [])
                if mechanisms:
                    result['opentargets_mechanisms'] = '|'.join(mechanisms[:10])
                
                print(f"OT:{result['opentargets_approved']}approved", end=' ')
            else:
                print("OT:0", end=' ')
            
            time.sleep(0.5)  # Rate limiting
        
        # Calculate summary flags
        total_approved = (
            result['dgidb_num_approved'] +
            result['chembl_approved_drugs'] + 
            result['opentargets_approved']
        )
        
        total_clinical = (
            result['chembl_clinical_compounds'] +
            result['opentargets_phase1'] +
            result['opentargets_phase2'] +
            result['opentargets_phase3']
        )
        
        result['has_approved_drugs'] = total_approved > 0
        result['has_clinical_drugs'] = total_clinical > 0
        result['has_bioactive_compounds'] = result['chembl_bioactive_compounds'] > 0
        
        # Determine druggability evidence level
        if result['has_approved_drugs']:
            result['druggability_evidence'] = 'HIGH - Approved drugs exist'
        elif result['has_clinical_drugs']:
            result['druggability_evidence'] = 'MEDIUM - Clinical stage compounds exist'
        elif result['has_bioactive_compounds'] or result['dgidb_num_interactions'] > 0:
            result['druggability_evidence'] = 'LOW - Bioactive compounds or interactions reported'
        else:
            result['druggability_evidence'] = 'NONE - No known drug interactions'
        
        print(f"✓ {result['druggability_evidence'].split('-')[0].strip()}")
        
        return result
    
    def process_step1_file(self, input_file: str = "Step1.csv") -> pd.DataFrame:
        """
        Process Step1.csv and add drug interaction data.
        
        Args:
            input_file: Path to Step1.csv
            
        Returns:
            DataFrame with all Step1 columns + drug data
        """
        print("="*70)
        print("DRUG-GENE INTERACTION ANALYSIS - STEP 2")
        print("="*70)
        
        # Read Step1.csv
        print(f"\nReading {input_file}...")
        step1_df = pd.read_csv(input_file)
        print(f"  ✓ Loaded {len(step1_df)} genes")
        
        # Process each gene
        print(f"\nQuerying drug databases for {len(step1_df)} genes...")
        print("Format: Gene... DGI:X ChEMBL:Y OT:Z ✓ Evidence\n")
        
        drug_results = []
        
        for idx, row in step1_df.iterrows():
            print(f"[{idx+1}/{len(step1_df)}]", end=' ')
            result = self.process_gene_row(row)
            drug_results.append(result)
        
        # Create drug data DataFrame
        drug_df = pd.DataFrame(drug_results)
        
        # Merge with Step1 data (keep gene_symbol from Step1)
        drug_df = drug_df.drop(columns=['gene_symbol'])
        
        # Merge on index
        step2_df = pd.concat([step1_df, drug_df], axis=1)
        
        # Print summary
        self.print_summary(step2_df)
        
        return step2_df
    
    def print_summary(self, df: pd.DataFrame):
        """Print summary of drug interaction findings."""
        print("\n" + "="*70)
        print("DRUG INTERACTION SUMMARY")
        print("="*70)
        
        total_genes = len(df)
        
        # Count genes with evidence
        high_evidence = (df['druggability_evidence'].str.contains('HIGH', na=False)).sum()
        medium_evidence = (df['druggability_evidence'].str.contains('MEDIUM', na=False)).sum()
        low_evidence = (df['druggability_evidence'].str.contains('LOW', na=False)).sum()
        no_evidence = (df['druggability_evidence'].str.contains('NONE', na=False)).sum()
        
        print(f"\nTotal genes analyzed: {total_genes}")
        print(f"\nDruggability Evidence Levels:")
        print(f"  🟢 HIGH (Approved drugs):        {high_evidence} ({high_evidence/total_genes*100:.1f}%)")
        print(f"  🟡 MEDIUM (Clinical compounds):  {medium_evidence} ({medium_evidence/total_genes*100:.1f}%)")
        print(f"  🟠 LOW (Bioactive compounds):    {low_evidence} ({low_evidence/total_genes*100:.1f}%)")
        print(f"  ⚪ NONE (No interactions):       {no_evidence} ({no_evidence/total_genes*100:.1f}%)")
        
        print(f"\nDatabase Coverage:")
        print(f"  DGIdb genes found:               {df['dgidb_found'].sum()} genes")
        print(f"  ChEMBL targets found:            {df['chembl_found'].sum()} genes")
        print(f"  Open Targets drugs found:        {df['opentargets_found'].sum()} genes")
        
        print(f"\nDrug Counts:")
        total_dgidb_approved = df['dgidb_num_approved'].sum()
        total_chembl_approved = df['chembl_approved_drugs'].sum()
        total_ot_approved = df['opentargets_approved'].sum()
        
        print(f"  Total approved drugs (DGIdb):    {total_dgidb_approved}")
        print(f"  Total approved drugs (ChEMBL):   {total_chembl_approved}")
        print(f"  Total approved drugs (OT):       {total_ot_approved}")
        print(f"  Total clinical compounds:        {df['chembl_clinical_compounds'].sum()}")
        print(f"  Total bioactive compounds:       {df['chembl_bioactive_compounds'].sum()}")
        print(f"  Total DGIdb interactions:        {df['dgidb_num_interactions'].sum()}")
        
        # Show top druggable genes
        print(f"\n🎯 Top 10 Druggable Genes (by approved drugs):")
        
        # Calculate total approved drugs per gene
        df_temp = df.copy()
        df_temp['total_approved'] = (
            df_temp['dgidb_num_approved'].fillna(0) +
            df_temp['chembl_approved_drugs'].fillna(0) +
            df_temp['opentargets_approved'].fillna(0)
        )
        
        top_druggable = df_temp[df_temp['total_approved'] > 0].sort_values(
            'total_approved', ascending=False
        ).head(10)
        
        if len(top_druggable) > 0:
            for idx, row in top_druggable.iterrows():
                print(f"  • {row['gene_symbol']:12s} "
                      f"Total: {int(row['total_approved']):3d} "
                      f"(DGI:{int(row['dgidb_num_approved']):2d}, "
                      f"ChEMBL:{int(row['chembl_approved_drugs']):2d}, "
                      f"OT:{int(row['opentargets_approved']):2d})")
        else:
            print("  No genes with approved drugs found.")
        
        print("="*70)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    # Check if Step1.csv exists
    if not Path("Step1.csv").exists():
        print("❌ Error: Step1.csv not found!")
        print("Please run Step 1 first to generate Step1.csv")
        exit(1)
    
    # Create drug interaction finder
    finder = DrugGeneInteractionFinder()
    
    # Process Step1.csv and add drug data
    step2_df = finder.process_step1_file("Step1.csv")
    
    # Save to Step2.csv
    output_path = Path("Step2.csv")
    step2_df.to_csv(output_path, index=False)
    
    print(f"\n✅ Results saved to: {output_path.absolute()}")
    print(f"\nTotal columns: {len(step2_df.columns)}")
    print(f"Total rows: {len(step2_df)}")
    
    # Display sample
    print("\n📋 Sample Drug Data (first 5 genes):")
    drug_columns = ['gene_symbol', 'dgidb_num_interactions', 'dgidb_num_approved',
                   'chembl_approved_drugs', 'opentargets_approved', 'druggability_evidence']
    if all(col in step2_df.columns for col in drug_columns):
        print(step2_df[drug_columns].head(5).to_string(index=False))
    
    print("\n✅ Step 2 Complete!")
    print("\nNext: Proceed to Step 3 (Protein Family Classification)")