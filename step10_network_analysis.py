"""
Step 10: Network Analysis (STRING Protein-Protein Interactions)
Reads Step9.csv, queries STRING for protein interactions,
identifies pathway context and drug target connectivity, saves to Step10.csv
"""

import pandas as pd
import requests
import time
from typing import Dict, List, Optional, Set
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('step10_network_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class NetworkAnalyzer:
    """Analyze protein-protein interaction networks using STRING database"""
    
    # Known drug target families for interaction analysis
    KNOWN_DRUG_TARGETS = {
        'GPCRs': ['HTR1A', 'HTR1B', 'HTR1D', 'HTR1F', 'HTR2A', 'HTR2C', 'ADRA1A', 'ADRA2A', 'DRD2', 'DRD3'],
        'Ion Channels': ['CACNA1A', 'CACNA1C', 'SCN1A', 'SCN9A', 'KCNMA1', 'KCNK18', 'TRPV1', 'TRPM8'],
        'Kinases': ['MAPK1', 'MAPK3', 'AKT1', 'PRKCA', 'PRKCE', 'CDK1', 'CDK2'],
        'Enzymes': ['NOS1', 'NOS2', 'NOS3', 'ACE', 'ACE2', 'MTHFR', 'COMT'],
        'Cytokines': ['TNF', 'IL1B', 'IL6', 'IL10', 'IFNG', 'TGFB1'],
        'Transporters': ['SLC6A4', 'SLC6A3', 'ATP1A2', 'ATP1A3']
    }
    
    def __init__(self):
        self.string_base_url = "https://string-db.org/api"
        self.species = 9606  # Human
        self.cache = {}
        
    def load_step9_data(self, filepath: str = 'Step9.csv') -> pd.DataFrame:
        """Load Step 9 data"""
        logger.info(f"Loading Step 9 data from {filepath}")
        df = pd.read_csv(filepath)
        logger.info(f"Successfully loaded {len(df)} entries")
        return df
    
    def _find_gene_column(self, df: pd.DataFrame) -> str:
        """Find gene identifier column"""
        for col in ['Gene', 'gene', 'Gene_Symbol', 'gene_symbol', 'Gene_Name', 'gene_name']:
            if col in df.columns:
                return col
        raise ValueError(f"Could not find gene column. Available: {df.columns.tolist()}")
    
    def get_string_interactions(self, gene_symbol: str, limit: int = 10) -> Optional[Dict]:
        """
        Get protein-protein interactions from STRING database
        
        Args:
            gene_symbol: Gene symbol
            limit: Maximum number of interactions to retrieve
        
        Returns:
            Dictionary with interaction data
        """
        
        logger.info(f"  Querying STRING for {gene_symbol}...")
        
        try:
            # STRING API endpoint for network interactions
            url = f"{self.string_base_url}/json/interaction_partners"
            
            params = {
                'identifiers': gene_symbol,
                'species': self.species,
                'limit': limit,
                'required_score': 400  # Medium confidence (0-1000)
            }
            
            response = requests.get(url, params=params, timeout=15)
            
            if response.status_code == 200:
                data = response.json()
                
                if data and len(data) > 0:
                    logger.info(f"    ? Found {len(data)} interactions")
                    return data
                else:
                    logger.info(f"    ? No interactions found")
                    return None
            else:
                logger.warning(f"    ? STRING API returned {response.status_code}")
                return None
                
        except requests.Timeout:
            logger.warning(f"    ? Timeout")
            return None
        except Exception as e:
            logger.warning(f"    ? Error: {e}")
            return None
    
    def get_functional_enrichment(self, gene_symbol: str) -> Optional[Dict]:
        """
        Get functional enrichment/pathway information from STRING
        
        Args:
            gene_symbol: Gene symbol
        
        Returns:
            Dictionary with enrichment data
        """
        
        try:
            url = f"{self.string_base_url}/json/enrichment"
            
            params = {
                'identifiers': gene_symbol,
                'species': self.species
            }
            
            response = requests.get(url, params=params, timeout=15)
            
            if response.status_code == 200:
                data = response.json()
                
                if data and len(data) > 0:
                    # Filter for relevant categories
                    pathways = [item for item in data if item.get('category') in ['KEGG', 'Reactome', 'WikiPathways']]
                    processes = [item for item in data if item.get('category') == 'Process']
                    
                    return {
                        'pathways': pathways[:5],  # Top 5 pathways
                        'processes': processes[:5]  # Top 5 processes
                    }
                
                return None
            else:
                return None
                
        except Exception as e:
            logger.debug(f"    Enrichment query failed: {e}")
            return None
    
    def analyze_network(self, gene_symbol: str) -> Dict:
        """Complete network analysis for a protein"""
        
        # Check cache
        cache_key = gene_symbol.upper()
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        logger.info(f"[Network Analysis] {gene_symbol}")
        
        network_data = {
            'Gene': gene_symbol,
            'STRING_Data_Available': False,
            'Total_Interactions': 0,
            'High_Confidence_Interactions': 0,
            'Interacting_Partners': [],
            'Drug_Target_Interactions': 0,
            'Drug_Target_Partners': [],
            'Known_Drug_Target_Families': [],
            'Network_Connectivity': 'Unknown',
            'Pathway_Context': [],
            'Biological_Processes': [],
            'Network_Hub': False
        }
        
        # Get interactions
        interactions = self.get_string_interactions(gene_symbol, limit=20)
        
        if interactions:
            network_data['STRING_Data_Available'] = True
            network_data['Total_Interactions'] = len(interactions)
            
            # Process interactions
            partners = []
            high_conf_count = 0
            drug_target_partners = []
            drug_target_families = set()
            
            for interaction in interactions:
                partner_name = interaction.get('preferredName_B', interaction.get('stringId_B', ''))
                score = interaction.get('score', 0)
                
                partners.append(partner_name)
                
                # High confidence interactions (score > 0.7)
                if score > 700:
                    high_conf_count += 1
                
                # Check if partner is a known drug target
                for family, targets in self.KNOWN_DRUG_TARGETS.items():
                    if partner_name in targets:
                        drug_target_partners.append(partner_name)
                        drug_target_families.add(family)
            
            network_data['Interacting_Partners'] = partners[:10]  # Top 10
            network_data['High_Confidence_Interactions'] = high_conf_count
            network_data['Drug_Target_Interactions'] = len(drug_target_partners)
            network_data['Drug_Target_Partners'] = drug_target_partners[:5]  # Top 5
            network_data['Known_Drug_Target_Families'] = list(drug_target_families)
            
            # Determine network connectivity
            if network_data['Total_Interactions'] >= 15:
                network_data['Network_Connectivity'] = 'Highly Connected (Hub)'
                network_data['Network_Hub'] = True
            elif network_data['Total_Interactions'] >= 10:
                network_data['Network_Connectivity'] = 'Well Connected'
            elif network_data['Total_Interactions'] >= 5:
                network_data['Network_Connectivity'] = 'Moderately Connected'
            else:
                network_data['Network_Connectivity'] = 'Poorly Connected'
            
            logger.info(f"  Interactions: {network_data['Total_Interactions']}, " +
                       f"Drug targets: {network_data['Drug_Target_Interactions']}")
        
        # Get pathway enrichment
        enrichment = self.get_functional_enrichment(gene_symbol)
        
        if enrichment:
            if enrichment.get('pathways'):
                pathways = [p.get('description', p.get('term', '')) for p in enrichment['pathways']]
                network_data['Pathway_Context'] = pathways[:3]  # Top 3
            
            if enrichment.get('processes'):
                processes = [p.get('description', p.get('term', '')) for p in enrichment['processes']]
                network_data['Biological_Processes'] = processes[:3]  # Top 3
        
        self.cache[cache_key] = network_data
        time.sleep(0.5)  # Rate limiting
        
        return network_data
    
    def calculate_network_scores(self, network_data: Dict) -> Dict:
        """Calculate network-based druggability scores"""
        
        scores = {
            'Step10_Connectivity_Score': 0,
            'Step10_Drug_Target_Proximity_Score': 0,
            'Step10_Pathway_Relevance_Score': 0,
            'Step10_Overall_Network_Score': 0,
            'Step10_Network_Confidence': 'Unknown'
        }
        
        if not network_data['STRING_Data_Available']:
            scores['Step10_Network_Confidence'] = 'Low - No Data'
            return scores
        
        # Connectivity score (0-5)
        total_interactions = network_data['Total_Interactions']
        if total_interactions >= 15:
            scores['Step10_Connectivity_Score'] = 5
        elif total_interactions >= 10:
            scores['Step10_Connectivity_Score'] = 4
        elif total_interactions >= 5:
            scores['Step10_Connectivity_Score'] = 3
        elif total_interactions >= 2:
            scores['Step10_Connectivity_Score'] = 2
        else:
            scores['Step10_Connectivity_Score'] = 1
        
        # Drug target proximity score (0-5)
        drug_target_count = network_data['Drug_Target_Interactions']
        if drug_target_count >= 5:
            scores['Step10_Drug_Target_Proximity_Score'] = 5
        elif drug_target_count >= 3:
            scores['Step10_Drug_Target_Proximity_Score'] = 4
        elif drug_target_count >= 2:
            scores['Step10_Drug_Target_Proximity_Score'] = 3
        elif drug_target_count >= 1:
            scores['Step10_Drug_Target_Proximity_Score'] = 2
        else:
            scores['Step10_Drug_Target_Proximity_Score'] = 1
        
        # Pathway relevance score (0-5)
        pathway_count = len(network_data['Pathway_Context'])
        process_count = len(network_data['Biological_Processes'])
        
        if pathway_count >= 3 or process_count >= 3:
            scores['Step10_Pathway_Relevance_Score'] = 5
        elif pathway_count >= 2 or process_count >= 2:
            scores['Step10_Pathway_Relevance_Score'] = 4
        elif pathway_count >= 1 or process_count >= 1:
            scores['Step10_Pathway_Relevance_Score'] = 3
        else:
            scores['Step10_Pathway_Relevance_Score'] = 2
        
        # Overall network score (0-10)
        scores['Step10_Overall_Network_Score'] = min(
            scores['Step10_Connectivity_Score'] + scores['Step10_Drug_Target_Proximity_Score'],
            10
        )
        
        # Confidence
        high_conf = network_data['High_Confidence_Interactions']
        if high_conf >= 5:
            scores['Step10_Network_Confidence'] = 'High - Strong Evidence'
        elif high_conf >= 2:
            scores['Step10_Network_Confidence'] = 'Medium - Moderate Evidence'
        else:
            scores['Step10_Network_Confidence'] = 'Low - Limited Evidence'
        
        return scores
    
    def assess_network_potential(self, network_data: Dict, scores: Dict) -> Dict:
        """Assess druggability based on network properties"""
        
        assessment = {
            'Step10_Network_Quality': 'Unknown',
            'Step10_Druggability_Via_Network': 'Unknown',
            'Step10_Key_Interactions': '',
            'Step10_Pathway_Summary': '',
            'Step10_Network_Strategy': ''
        }
        
        if not network_data['STRING_Data_Available']:
            return assessment
        
        overall_score = scores['Step10_Overall_Network_Score']
        
        # Network quality
        if overall_score >= 8:
            assessment['Step10_Network_Quality'] = 'Excellent'
        elif overall_score >= 6:
            assessment['Step10_Network_Quality'] = 'Very Good'
        elif overall_score >= 4:
            assessment['Step10_Network_Quality'] = 'Good'
        else:
            assessment['Step10_Network_Quality'] = 'Moderate'
        
        # Druggability via network
        drug_target_count = network_data['Drug_Target_Interactions']
        if drug_target_count >= 3:
            assessment['Step10_Druggability_Via_Network'] = 'High - Connected to drug targets'
        elif drug_target_count >= 1:
            assessment['Step10_Druggability_Via_Network'] = 'Moderate - Some drug target links'
        else:
            assessment['Step10_Druggability_Via_Network'] = 'Low - No direct drug target links'
        
        # Key interactions
        if network_data['Drug_Target_Partners']:
            partners_str = ', '.join(network_data['Drug_Target_Partners'][:3])
            families_str = ', '.join(network_data['Known_Drug_Target_Families'][:3])
            assessment['Step10_Key_Interactions'] = f"Drug targets: {partners_str} ({families_str})"
        else:
            top_partners = network_data['Interacting_Partners'][:3]
            assessment['Step10_Key_Interactions'] = f"Top partners: {', '.join(top_partners)}"
        
        # Pathway summary
        if network_data['Pathway_Context']:
            assessment['Step10_Pathway_Summary'] = '; '.join(network_data['Pathway_Context'][:2])
        elif network_data['Biological_Processes']:
            assessment['Step10_Pathway_Summary'] = '; '.join(network_data['Biological_Processes'][:2])
        else:
            assessment['Step10_Pathway_Summary'] = 'No pathway data'
        
        # Network strategy
        if network_data['Network_Hub']:
            assessment['Step10_Network_Strategy'] = 'Target hub protein - broad pathway impact'
        elif drug_target_count >= 2:
            assessment['Step10_Network_Strategy'] = 'Target alongside known drugs - combination therapy'
        elif network_data['Total_Interactions'] >= 5:
            assessment['Step10_Network_Strategy'] = 'Target within pathway - modulate network'
        else:
            assessment['Step10_Network_Strategy'] = 'Direct targeting - limited network effects'
        
        return assessment
    
    def add_network_data(self, step9_df: pd.DataFrame) -> pd.DataFrame:
        """Add network analysis to Step 9 DataFrame"""
        
        gene_col = self._find_gene_column(step9_df)
        unique_genes = step9_df[gene_col].unique().tolist()
        
        logger.info(f"\nAnalyzing networks for {len(unique_genes)} genes...")
        logger.info("="*80)
        
        all_network_data = []
        
        for idx, gene in enumerate(unique_genes, 1):
            network_data = self.analyze_network(gene)
            scores = self.calculate_network_scores(network_data)
            assessment = self.assess_network_potential(network_data, scores)
            
            combined = {
                gene_col: gene,
                'Step10_STRING_Data_Available': network_data['STRING_Data_Available'],
                'Step10_Total_Interactions': network_data['Total_Interactions'],
                'Step10_High_Confidence_Interactions': network_data['High_Confidence_Interactions'],
                'Step10_Drug_Target_Interactions': network_data['Drug_Target_Interactions'],
                'Step10_Drug_Target_Partners': '; '.join(network_data['Drug_Target_Partners']),
                'Step10_Drug_Target_Families': '; '.join(network_data['Known_Drug_Target_Families']),
                'Step10_Network_Connectivity': network_data['Network_Connectivity'],
                'Step10_Network_Hub': network_data['Network_Hub'],
                'Step10_Top_Partners': '; '.join(network_data['Interacting_Partners'][:5]),
                'Step10_Pathways': '; '.join(network_data['Pathway_Context']),
                'Step10_Biological_Processes': '; '.join(network_data['Biological_Processes']),
                **scores,
                **assessment
            }
            
            all_network_data.append(combined)
            
            logger.info(f"  [{idx}/{len(unique_genes)}] {gene}: " +
                       f"Interactions={network_data['Total_Interactions']}, " +
                       f"Network Score={scores['Step10_Overall_Network_Score']}/10")
        
        network_df = pd.DataFrame(all_network_data)
        combined_df = step9_df.merge(network_df, on=gene_col, how='left')
        
        logger.info(f"\n? Added network data to {len(combined_df)} rows")
        return combined_df
    
    def generate_summary_report(self, df: pd.DataFrame):
        """Generate summary statistics"""
        
        gene_col = self._find_gene_column(df)
        protein_level = df.drop_duplicates(subset=[gene_col])
        
        logger.info("\n" + "="*80)
        logger.info("STEP 10 NETWORK ANALYSIS SUMMARY")
        logger.info("="*80)
        
        logger.info(f"\nData Availability:")
        available = protein_level['Step10_STRING_Data_Available'].sum()
        logger.info(f"  STRING data: {available}/{len(protein_level)}")
        
        logger.info(f"\nNetwork Connectivity:")
        for conn, count in protein_level['Step10_Network_Connectivity'].value_counts().items():
            logger.info(f"  {conn}: {count}")
        
        logger.info(f"\nNetwork Hubs:")
        hubs = protein_level['Step10_Network_Hub'].sum()
        logger.info(f"  Hub proteins: {hubs}")
        
        logger.info(f"\nDrug Target Connectivity:")
        has_drug_targets = (protein_level['Step10_Drug_Target_Interactions'] > 0).sum()
        logger.info(f"  Connected to drug targets: {has_drug_targets}/{available}")
        
        logger.info(f"\nMost Connected Drug Target Families:")
        all_families = []
        for families in protein_level['Step10_Drug_Target_Families'].dropna():
            if families and families != '':
                all_families.extend(families.split('; '))
        
        if all_families:
            from collections import Counter
            family_counts = Counter(all_families)
            for family, count in family_counts.most_common(5):
                logger.info(f"  {family}: {count}")
        
        return protein_level


def main():
    """Main execution function"""
    
    STEP9_INPUT = 'Step9.csv'
    STEP10_OUTPUT = 'Step10.csv'
    
    try:
        logger.info("="*80)
        logger.info("STEP 10: NETWORK ANALYSIS (STRING)")
        logger.info("="*80)
        
        analyzer = NetworkAnalyzer()
        step9_df = analyzer.load_step9_data(STEP9_INPUT)
        step10_df = analyzer.add_network_data(step9_df)
        
        logger.info(f"\nSaving to {STEP10_OUTPUT}...")
        step10_df.to_csv(STEP10_OUTPUT, index=False)
        logger.info(f"? Saved {len(step10_df)} rows")
        
        protein_level = analyzer.generate_summary_report(step10_df)
        
        gene_col = analyzer._find_gene_column(step10_df)
        summary_cols = [
            gene_col, 'Step10_Total_Interactions', 'Step10_Drug_Target_Interactions',
            'Step10_Network_Connectivity', 'Step10_Overall_Network_Score'
        ]
        
        available_cols = [col for col in summary_cols if col in protein_level.columns]
        summary = protein_level[available_cols].sort_values(
            'Step10_Overall_Network_Score', ascending=False
        )
        
        summary.to_csv('Step10_protein_summary.csv', index=False)
        logger.info(f"\n? Summary saved")
        
        logger.info("\n" + "="*80)
        logger.info("STEP 10 COMPLETE!")
        logger.info("="*80)
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()