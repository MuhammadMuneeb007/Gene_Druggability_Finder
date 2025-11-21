"""
Step 9: Pharos Analysis using curl (subprocess)
Works around requests library being blocked by using curl directly
"""

import pandas as pd
import subprocess
import json
import time
from typing import Dict, Optional
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('step9_pharos_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class PharosAnalyzer:
    """Analyze target development level using Pharos GraphQL API via curl"""
    
    TDL_DESCRIPTIONS = {
        'Tclin': 'Clinical - Targets with approved drugs or clinical candidates',
        'Tchem': 'Chemical - Targets with bioactive compounds (IC50/Ki <30nM)',
        'Tbio': 'Biological - Targets with disease association but no compounds',
        'Tdark': 'Dark - Minimal knowledge, no compounds, no disease association'
    }
    
    TDL_DRUGGABILITY = {
        'Tclin': {'score': 5, 'quality': 'Excellent'},
        'Tchem': {'score': 4, 'quality': 'Very Good'},
        'Tbio': {'score': 3, 'quality': 'Moderate'},
        'Tdark': {'score': 2, 'quality': 'Poor'}
    }
    
    def __init__(self):
        self.pharos_url = "https://pharos-api.ncats.io/graphql"
        self.cache = {}
        
    def load_step8_data(self, filepath: str = 'Step8.csv') -> pd.DataFrame:
        """Load Step 8 data"""
        logger.info(f"Loading Step 8 data from {filepath}")
        df = pd.read_csv(filepath)
        logger.info(f"Successfully loaded {len(df)} entries")
        return df
    
    def _find_gene_column(self, df: pd.DataFrame) -> str:
        """Find gene identifier column"""
        for col in ['Gene', 'gene', 'Gene_Symbol', 'gene_symbol', 'Gene_Name', 'gene_name']:
            if col in df.columns:
                return col
        raise ValueError(f"Could not find gene column. Available: {df.columns.tolist()}")
    
    def query_pharos_with_curl(self, gene_symbol: str) -> Optional[Dict]:
        """
        Query Pharos GraphQL API using curl subprocess
        """
        
        logger.info(f"  Querying Pharos for {gene_symbol}...")
        
        # GraphQL query
        query = {
            "query": f'query {{ target(q: {{sym: "{gene_symbol.upper()}"}}) {{ name sym tdl fam description novelty }} }}'
        }
        
        query_json = json.dumps(query)
        
        # Build curl command
        curl_command = [
            'curl',
            '--request', 'POST',
            '--header', 'content-type: application/json',
            '--url', self.pharos_url,
            '--data', query_json,
            '--silent',  # Silent mode
            '--max-time', '15'  # Timeout after 15 seconds
        ]
        
        try:
            # Execute curl
            result = subprocess.run(
                curl_command,
                capture_output=True,
                text=True,
                timeout=20
            )
            
            if result.returncode == 0:
                # Parse JSON response
                response_data = json.loads(result.stdout)
                
                if 'data' in response_data and 'target' in response_data['data']:
                    target = response_data['data']['target']
                    
                    if target:
                        logger.info(f"    ? Found: {target.get('name', 'N/A')} " +
                                  f"(TDL: {target.get('tdl', 'Unknown')})")
                        return target
                    else:
                        logger.info(f"    ? No target found")
                        return None
                else:
                    logger.warning(f"    ? Unexpected response format")
                    return None
            else:
                logger.warning(f"    ? curl failed with code {result.returncode}")
                logger.warning(f"    Error: {result.stderr[:200]}")
                return None
                
        except subprocess.TimeoutExpired:
            logger.warning(f"    ? Timeout")
            return None
        except json.JSONDecodeError as e:
            logger.warning(f"    ? JSON parse error: {e}")
            return None
        except Exception as e:
            logger.warning(f"    ? Error: {e}")
            return None
    
    def analyze_pharos_target(self, gene_symbol: str) -> Dict:
        """Analyze target using Pharos data"""
        
        cache_key = gene_symbol.upper()
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        logger.info(f"[Pharos Analysis] {gene_symbol}")
        
        pharos_data = {
            'Gene': gene_symbol,
            'Pharos_Data_Available': False,
            'Target_Name': 'Unknown',
            'Target_Development_Level': 'Unknown',
            'TDL_Description': 'Unknown',
            'Protein_Family': 'Unknown',
            'Target_Description': 'N/A',
            'Novelty_Score': 0.5,
            'Development_Status': 'Unknown'
        }
        
        # Query Pharos via curl
        target_data = self.query_pharos_with_curl(gene_symbol)
        
        if target_data:
            pharos_data['Pharos_Data_Available'] = True
            pharos_data['Target_Name'] = target_data.get('name', 'Unknown')
            
            tdl = target_data.get('tdl', 'Unknown')
            pharos_data['Target_Development_Level'] = tdl
            pharos_data['TDL_Description'] = self.TDL_DESCRIPTIONS.get(tdl, 'Unknown')
            
            pharos_data['Protein_Family'] = target_data.get('fam') or 'Unknown'
            
            description = target_data.get('description', '')
            pharos_data['Target_Description'] = description[:200] if description else 'N/A'
            
            pharos_data['Novelty_Score'] = target_data.get('novelty', 0.5)
            
            if tdl == 'Tclin':
                pharos_data['Development_Status'] = 'Clinical stage or approved'
            elif tdl == 'Tchem':
                pharos_data['Development_Status'] = 'Chemical tool stage'
            elif tdl == 'Tbio':
                pharos_data['Development_Status'] = 'Biological validation stage'
            elif tdl == 'Tdark':
                pharos_data['Development_Status'] = 'Early/undeveloped'
            
            logger.info(f"  Summary: TDL={tdl}, Family={pharos_data['Protein_Family']}, " +
                       f"Novelty={pharos_data['Novelty_Score']:.5f}")
        else:
            logger.warning(f"  ! No Pharos data for {gene_symbol}")
        
        self.cache[cache_key] = pharos_data
        time.sleep(0.5)
        
        return pharos_data
    
    def calculate_pharos_scores(self, pharos_data: Dict) -> Dict:
        """Calculate druggability scores"""
        
        scores = {
            'Step9_TDL_Score': 0,
            'Step9_Novelty_Score': 0,
            'Step9_Overall_Pharos_Score': 0,
            'Step9_Druggability_Confidence': 'Unknown'
        }
        
        if not pharos_data['Pharos_Data_Available']:
            scores['Step9_Druggability_Confidence'] = 'Low - No Data'
            return scores
        
        tdl = pharos_data['Target_Development_Level']
        tdl_info = self.TDL_DRUGGABILITY.get(tdl, {'score': 0})
        scores['Step9_TDL_Score'] = tdl_info['score']
        
        # Novelty score (0-5, lower novelty = more studied = better)
        novelty = pharos_data['Novelty_Score']
        if novelty < 0.1:
            scores['Step9_Novelty_Score'] = 5
        elif novelty < 0.3:
            scores['Step9_Novelty_Score'] = 4
        elif novelty < 0.5:
            scores['Step9_Novelty_Score'] = 3
        elif novelty < 0.7:
            scores['Step9_Novelty_Score'] = 2
        else:
            scores['Step9_Novelty_Score'] = 1
        
        # Overall score (0-10)
        scores['Step9_Overall_Pharos_Score'] = min(
            scores['Step9_TDL_Score'] + scores['Step9_Novelty_Score'],
            10
        )
        
        if tdl in ['Tclin', 'Tchem']:
            scores['Step9_Druggability_Confidence'] = 'High - Validated Target'
        elif tdl == 'Tbio':
            scores['Step9_Druggability_Confidence'] = 'Medium - Biological Evidence'
        elif tdl == 'Tdark':
            scores['Step9_Druggability_Confidence'] = 'Low - Understudied'
        
        return scores
    
    def assess_target_potential(self, pharos_data: Dict, scores: Dict) -> Dict:
        """Assess target potential"""
        
        assessment = {
            'Step9_Target_Quality': 'Unknown',
            'Step9_Drug_Development_Potential': 'Unknown',
            'Step9_Recommended_Strategy': ''
        }
        
        if not pharos_data['Pharos_Data_Available']:
            return assessment
        
        tdl = pharos_data['Target_Development_Level']
        overall_score = scores['Step9_Overall_Pharos_Score']
        
        if overall_score >= 8:
            assessment['Step9_Target_Quality'] = 'Excellent'
        elif overall_score >= 6:
            assessment['Step9_Target_Quality'] = 'Very Good'
        elif overall_score >= 4:
            assessment['Step9_Target_Quality'] = 'Good'
        else:
            assessment['Step9_Target_Quality'] = 'Moderate'
        
        tdl_info = self.TDL_DRUGGABILITY.get(tdl, {'quality': 'Unknown'})
        assessment['Step9_Drug_Development_Potential'] = tdl_info['quality']
        
        if tdl == 'Tclin':
            assessment['Step9_Recommended_Strategy'] = 'Repurposing or improved drug design'
        elif tdl == 'Tchem':
            assessment['Step9_Recommended_Strategy'] = 'Lead optimization from known compounds'
        elif tdl == 'Tbio':
            assessment['Step9_Recommended_Strategy'] = 'HTS screening, probe development'
        elif tdl == 'Tdark':
            assessment['Step9_Recommended_Strategy'] = 'Target validation required'
        else:
            assessment['Step9_Recommended_Strategy'] = 'Further assessment needed'
        
        return assessment
    
    def add_pharos_data(self, step8_df: pd.DataFrame) -> pd.DataFrame:
        """Add Pharos data to Step 8"""
        
        gene_col = self._find_gene_column(step8_df)
        unique_genes = step8_df[gene_col].unique().tolist()
        
        logger.info(f"\nQuerying Pharos for {len(unique_genes)} genes...")
        logger.info("="*80)
        
        all_data = []
        
        for idx, gene in enumerate(unique_genes, 1):
            pharos_data = self.analyze_pharos_target(gene)
            scores = self.calculate_pharos_scores(pharos_data)
            assessment = self.assess_target_potential(pharos_data, scores)
            
            combined = {
                gene_col: gene,
                'Step9_Pharos_Data_Available': pharos_data['Pharos_Data_Available'],
                'Step9_Target_Name': pharos_data['Target_Name'],
                'Step9_Target_Development_Level': pharos_data['Target_Development_Level'],
                'Step9_TDL_Description': pharos_data['TDL_Description'],
                'Step9_Protein_Family': pharos_data['Protein_Family'],
                'Step9_Novelty_Score': round(pharos_data['Novelty_Score'], 5),
                'Step9_Development_Status': pharos_data['Development_Status'],
                'Step9_Target_Description': pharos_data['Target_Description'],
                **scores,
                **assessment
            }
            
            all_data.append(combined)
            
            logger.info(f"  [{idx}/{len(unique_genes)}] {gene}: TDL={pharos_data['Target_Development_Level']}, " +
                       f"Score={scores['Step9_Overall_Pharos_Score']}/10")
        
        result_df = pd.DataFrame(all_data)
        combined_df = step8_df.merge(result_df, on=gene_col, how='left')
        
        logger.info(f"\n? Added Pharos data to {len(combined_df)} rows")
        return combined_df
    
    def generate_summary_report(self, df: pd.DataFrame):
        """Generate summary"""
        
        gene_col = self._find_gene_column(df)
        protein_level = df.drop_duplicates(subset=[gene_col])
        
        logger.info("\n" + "="*80)
        logger.info("STEP 9 PHAROS ANALYSIS SUMMARY")
        logger.info("="*80)
        
        logger.info(f"\nData Availability:")
        available = protein_level['Step9_Pharos_Data_Available'].sum()
        logger.info(f"  Pharos data: {available}/{len(protein_level)}")
        
        logger.info(f"\nTarget Development Level (TDL):")
        for tdl, count in protein_level['Step9_Target_Development_Level'].value_counts().items():
            logger.info(f"  {tdl}: {count}")
            if tdl in self.TDL_DESCRIPTIONS:
                logger.info(f"    ? {self.TDL_DESCRIPTIONS[tdl]}")
        
        logger.info(f"\nProtein Families:")
        for fam, count in protein_level['Step9_Protein_Family'].value_counts().head(5).items():
            logger.info(f"  {fam}: {count}")
        
        return protein_level


def main():
    STEP8_INPUT = 'Step8.csv'
    STEP9_OUTPUT = 'Step9.csv'
    
    try:
        logger.info("="*80)
        logger.info("STEP 9: PHAROS/TCRD ANALYSIS (via curl)")
        logger.info("="*80)
        
        analyzer = PharosAnalyzer()
        step8_df = analyzer.load_step8_data(STEP8_INPUT)
        step9_df = analyzer.add_pharos_data(step8_df)
        
        logger.info(f"\nSaving to {STEP9_OUTPUT}...")
        step9_df.to_csv(STEP9_OUTPUT, index=False)
        logger.info(f"? Saved {len(step9_df)} rows")
        
        protein_level = analyzer.generate_summary_report(step9_df)
        
        gene_col = analyzer._find_gene_column(step9_df)
        summary = protein_level[[
            gene_col, 'Step9_Target_Development_Level', 
            'Step9_Overall_Pharos_Score', 'Step9_Protein_Family'
        ]].sort_values('Step9_Overall_Pharos_Score', ascending=False)
        
        summary.to_csv('Step9_protein_summary.csv', index=False)
        logger.info(f"\n? Summary saved")
        
        logger.info("\n" + "="*80)
        logger.info("STEP 9 COMPLETE!")
        logger.info("="*80)
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()