"""
Step 7: Comprehensive Tissue Expression Analysis (Disease-Agnostic)
Reads Step6.csv, queries GTEx for expression across ALL tissues,
provides complete tissue expression profile for druggability assessment
"""

import pandas as pd
import requests
import time
from typing import Dict, List, Optional
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('step7_expression_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class ComprehensiveTissueExpressionAnalyzer:
    """Analyze tissue expression patterns across all tissues for druggability"""
    
    # Tissue categories for organization
    TISSUE_CATEGORIES = {
        'Brain': ['Brain_'],
        'Blood': ['Whole_Blood'],
        'Vascular': ['Artery_', 'Heart_'],
        'Immune': ['Spleen', 'Lymph'],
        'Gastrointestinal': ['Stomach', 'Colon', 'Esophagus', 'Small_Intestine', 'Liver', 'Pancreas'],
        'Reproductive': ['Testis', 'Ovary', 'Uterus', 'Prostate', 'Vagina', 'Breast'],
        'Kidney': ['Kidney'],
        'Lung': ['Lung'],
        'Skin': ['Skin'],
        'Muscle': ['Muscle'],
        'Adipose': ['Adipose'],
        'Nerve': ['Nerve'],
        'Other': []
    }
    
    def __init__(self):
        self.gtex_base_url = "https://gtexportal.org/api/v2"
        self.cache = {}
        
    def load_step6_data(self, filepath: str = 'Step6.csv') -> pd.DataFrame:
        """Load Step 6 data"""
        logger.info(f"Loading Step 6 data from {filepath}")
        df = pd.read_csv(filepath)
        logger.info(f"Successfully loaded {len(df)} entries")
        return df
    
    def _find_gene_column(self, df: pd.DataFrame) -> str:
        """Find gene identifier column"""
        for col in ['Gene', 'gene', 'Gene_Symbol', 'gene_symbol', 'Gene_Name', 'gene_name']:
            if col in df.columns:
                return col
        raise ValueError(f"Could not find gene column. Available: {df.columns.tolist()}")
    
    def get_gencode_id(self, gene_symbol: str) -> Optional[str]:
        """Get versioned Gencode ID from gene symbol"""
        
        url = f"{self.gtex_base_url}/reference/gene"
        params = {
            'geneId': gene_symbol,
            'format': 'json'
        }
        
        try:
            response = requests.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                
                if 'data' in data and len(data['data']) > 0:
                    gene_info = data['data'][0]
                    gencode_id = gene_info.get('gencodeId')
                    
                    if gencode_id:
                        logger.info(f"    Found Gencode ID: {gencode_id}")
                        return gencode_id
            
            logger.warning(f"    No Gencode ID found for {gene_symbol}")
            return None
            
        except Exception as e:
            logger.warning(f"    Error getting Gencode ID: {e}")
            return None
    
    def categorize_tissue(self, tissue_id: str) -> str:
        """Categorize tissue by system"""
        
        for category, patterns in self.TISSUE_CATEGORIES.items():
            for pattern in patterns:
                if pattern in tissue_id:
                    return category
        return 'Other'
    
    def get_comprehensive_expression(self, gene_symbol: str) -> Dict:
        """Get comprehensive tissue expression profile"""
        
        # Check cache
        if gene_symbol in self.cache:
            return self.cache[gene_symbol]
        
        logger.info(f"  Querying GTEx for {gene_symbol}...")
        
        expression_data = {
            'Gene': gene_symbol,
            'GTEx_Data_Available': False,
            'Total_Tissues': 0,
            'Tissues_With_Expression': 0,
            'Tissues_High_Expression': 0,
            'Tissues_Medium_Expression': 0,
            'Tissues_Low_Expression': 0,
            'Max_TPM': 0.0,
            'Mean_TPM': 0.0,
            'Median_TPM': 0.0,
            'Tissue_Specificity': 'Unknown',
            'Top_Expressed_Tissues': [],
            'Expression_By_Category': {},
            'All_Tissue_Expression': {},
            'Brain_Expression': 'No',
            'Brain_Max_TPM': 0.0,
            'Blood_Expression': 'No',
            'Blood_TPM': 0.0,
            'Accessible_Expression': 'Unknown',
        }
        
        # Get Gencode ID
        gencode_id = self.get_gencode_id(gene_symbol)
        
        if not gencode_id:
            self.cache[gene_symbol] = expression_data
            return expression_data
        
        # Get expression data
        url = f"{self.gtex_base_url}/expression/medianGeneExpression"
        params = {
            'gencodeId': gencode_id,
            'datasetId': 'gtex_v8',
            'format': 'json'
        }
        
        try:
            response = requests.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                
                if 'data' in data and len(data['data']) > 0:
                    expression_data['GTEx_Data_Available'] = True
                    tissue_data = data['data']
                    
                    all_tpms = []
                    high_expression = []
                    medium_expression = []
                    low_expression = []
                    category_max_tpm = {}
                    tissue_expression = {}
                    
                    brain_tpms = []
                    blood_tpm = 0.0
                    
                    # Process all tissues
                    for tissue in tissue_data:
                        tissue_id = tissue.get('tissueSiteDetailId', '')
                        median_tpm = tissue.get('median', 0.0)
                        
                        tissue_expression[tissue_id] = median_tpm
                        
                        if median_tpm > 0:
                            all_tpms.append(median_tpm)
                        
                        # Categorize by expression level
                        if median_tpm >= 10:
                            high_expression.append(tissue_id)
                        elif median_tpm >= 1:
                            medium_expression.append(tissue_id)
                        elif median_tpm >= 0.1:
                            low_expression.append(tissue_id)
                        
                        # Track by category
                        category = self.categorize_tissue(tissue_id)
                        if category not in category_max_tpm:
                            category_max_tpm[category] = 0.0
                        category_max_tpm[category] = max(category_max_tpm[category], median_tpm)
                        
                        # Track specific systems
                        if tissue_id.startswith('Brain_'):
                            brain_tpms.append(median_tpm)
                        
                        if tissue_id == 'Whole_Blood':
                            blood_tpm = median_tpm
                    
                    # Calculate overall statistics
                    expression_data['Total_Tissues'] = len(tissue_data)
                    expression_data['Tissues_With_Expression'] = len(all_tpms)
                    expression_data['Tissues_High_Expression'] = len(high_expression)
                    expression_data['Tissues_Medium_Expression'] = len(medium_expression)
                    expression_data['Tissues_Low_Expression'] = len(low_expression)
                    
                    if all_tpms:
                        expression_data['Max_TPM'] = max(all_tpms)
                        expression_data['Mean_TPM'] = sum(all_tpms) / len(all_tpms)
                        sorted_tpms = sorted(all_tpms)
                        mid = len(sorted_tpms) // 2
                        expression_data['Median_TPM'] = sorted_tpms[mid]
                    
                    # Tissue specificity
                    tissues_expressed = len([tpm for tpm in all_tpms if tpm > 0.1])
                    if tissues_expressed <= 5:
                        expression_data['Tissue_Specificity'] = 'Highly Specific'
                    elif tissues_expressed <= 15:
                        expression_data['Tissue_Specificity'] = 'Moderately Specific'
                    elif tissues_expressed <= 30:
                        expression_data['Tissue_Specificity'] = 'Broadly Expressed'
                    else:
                        expression_data['Tissue_Specificity'] = 'Ubiquitously Expressed'
                    
                    # Top expressed tissues (top 5)
                    sorted_tissues = sorted(tissue_expression.items(), key=lambda x: x[1], reverse=True)
                    top_5 = [(tissue, tpm) for tissue, tpm in sorted_tissues[:5] if tpm > 0]
                    expression_data['Top_Expressed_Tissues'] = [f"{t}:{tpm:.1f}" for t, tpm in top_5]
                    
                    # Expression by category
                    expression_data['Expression_By_Category'] = category_max_tpm
                    expression_data['All_Tissue_Expression'] = tissue_expression
                    
                    # Brain expression
                    if brain_tpms and max(brain_tpms) > 0.1:
                        expression_data['Brain_Expression'] = 'Yes'
                        expression_data['Brain_Max_TPM'] = max(brain_tpms)
                    
                    # Blood expression
                    if blood_tpm > 0.1:
                        expression_data['Blood_Expression'] = 'Yes'
                        expression_data['Blood_TPM'] = blood_tpm
                    
                    # Accessible expression (blood/secreted/membrane)
                    if blood_tpm > 1 or expression_data['Brain_Expression'] == 'Yes':
                        expression_data['Accessible_Expression'] = 'Yes - Systemic/CNS'
                    elif blood_tpm > 0.1:
                        expression_data['Accessible_Expression'] = 'Low - Limited Systemic'
                    else:
                        expression_data['Accessible_Expression'] = 'No - Intracellular Only'
                    
                    logger.info(f"    ? Expressed in {tissues_expressed}/{len(tissue_data)} tissues")
                    logger.info(f"    ? Max TPM: {expression_data['Max_TPM']:.2f}")
                    logger.info(f"    ? Specificity: {expression_data['Tissue_Specificity']}")
                    logger.info(f"    ? Top tissue: {top_5[0][0] if top_5 else 'None'}")
                else:
                    logger.info(f"    ? No expression data returned")
            else:
                logger.warning(f"    ? GTEx API returned {response.status_code}")
        
        except requests.Timeout:
            logger.warning(f"    ? Request timed out")
        except Exception as e:
            logger.warning(f"    ? Error: {e}")
        
        self.cache[gene_symbol] = expression_data
        time.sleep(0.5)
        
        return expression_data
    
    def calculate_druggability_scores(self, expression_data: Dict) -> Dict:
        """Calculate druggability scores based on expression profile"""
        
        scores = {
            'Step7_Expression_Level_Score': 0,
            'Step7_Tissue_Accessibility_Score': 0,
            'Step7_Expression_Breadth_Score': 0,
            'Step7_Overall_Expression_Score': 0,
            'Step7_Expression_Confidence': 'Unknown'
        }
        
        # Expression level score (0-5) - higher is generally better for validation
        max_tpm = expression_data['Max_TPM']
        if max_tpm >= 50:
            scores['Step7_Expression_Level_Score'] = 5
        elif max_tpm >= 10:
            scores['Step7_Expression_Level_Score'] = 4
        elif max_tpm >= 1:
            scores['Step7_Expression_Level_Score'] = 3
        elif max_tpm >= 0.1:
            scores['Step7_Expression_Level_Score'] = 2
        elif max_tpm > 0:
            scores['Step7_Expression_Level_Score'] = 1
        
        # Tissue accessibility score (0-5) - easier to access = higher score
        if expression_data['Blood_Expression'] == 'Yes':
            if expression_data['Blood_TPM'] >= 10:
                scores['Step7_Tissue_Accessibility_Score'] = 5
            elif expression_data['Blood_TPM'] >= 1:
                scores['Step7_Tissue_Accessibility_Score'] = 4
            else:
                scores['Step7_Tissue_Accessibility_Score'] = 3
        elif expression_data['Brain_Expression'] == 'Yes':
            scores['Step7_Tissue_Accessibility_Score'] = 3  # CNS accessible
        else:
            # Check for other accessible tissues
            accessible_categories = ['Vascular', 'Lung', 'Skin']
            category_expr = expression_data.get('Expression_By_Category', {})
            if any(category_expr.get(cat, 0) > 1 for cat in accessible_categories):
                scores['Step7_Tissue_Accessibility_Score'] = 2
            else:
                scores['Step7_Tissue_Accessibility_Score'] = 1
        
        # Expression breadth score (0-5) - specificity can be good (less off-target)
        specificity = expression_data['Tissue_Specificity']
        if specificity == 'Highly Specific':
            scores['Step7_Expression_Breadth_Score'] = 5  # Specific is good
        elif specificity == 'Moderately Specific':
            scores['Step7_Expression_Breadth_Score'] = 4
        elif specificity == 'Broadly Expressed':
            scores['Step7_Expression_Breadth_Score'] = 3
        else:  # Ubiquitous
            scores['Step7_Expression_Breadth_Score'] = 2  # More off-target risk
        
        # Overall expression score (0-10)
        scores['Step7_Overall_Expression_Score'] = min(
            scores['Step7_Expression_Level_Score'] + 
            scores['Step7_Tissue_Accessibility_Score'],
            10
        )
        
        # Confidence
        if expression_data['GTEx_Data_Available']:
            scores['Step7_Expression_Confidence'] = 'High - GTEx V8'
        else:
            scores['Step7_Expression_Confidence'] = 'Low - No Data'
        
        return scores
    
    def assess_druggability(self, expression_data: Dict, scores: Dict) -> Dict:
        """Assess overall druggability based on expression"""
        
        assessment = {
            'Step7_Expression_Quality': 'Unknown',
            'Step7_Expression_Summary': '',
            'Step7_Target_Validation_Feasibility': '',
            'Step7_Drug_Delivery_Considerations': '',
            'Step7_Off_Target_Risk': ''
        }
        
        overall_score = scores['Step7_Overall_Expression_Score']
        max_tpm = expression_data['Max_TPM']
        specificity = expression_data['Tissue_Specificity']
        
        # Overall quality
        if overall_score >= 8 and expression_data['GTEx_Data_Available']:
            assessment['Step7_Expression_Quality'] = 'Excellent'
        elif overall_score >= 6:
            assessment['Step7_Expression_Quality'] = 'Good'
        elif overall_score >= 4:
            assessment['Step7_Expression_Quality'] = 'Moderate'
        elif overall_score >= 2:
            assessment['Step7_Expression_Quality'] = 'Poor'
        else:
            assessment['Step7_Expression_Quality'] = 'Very Poor'
        
        # Summary
        tissues_expressed = expression_data['Tissues_With_Expression']
        total_tissues = expression_data['Total_Tissues']
        top_tissues = expression_data['Top_Expressed_Tissues'][:3]
        
        assessment['Step7_Expression_Summary'] = (
            f"Expressed in {tissues_expressed}/{total_tissues} tissues; "
            f"Max TPM: {max_tpm:.1f}; {specificity}; "
            f"Top: {', '.join([t.split(':')[0] for t in top_tissues]) if top_tissues else 'None'}"
        )
        
        # Target validation feasibility
        if max_tpm >= 10:
            assessment['Step7_Target_Validation_Feasibility'] = 'Excellent - High expression enables robust validation'
        elif max_tpm >= 1:
            assessment['Step7_Target_Validation_Feasibility'] = 'Good - Sufficient for most validation methods'
        elif max_tpm >= 0.1:
            assessment['Step7_Target_Validation_Feasibility'] = 'Moderate - May require sensitive methods'
        else:
            assessment['Step7_Target_Validation_Feasibility'] = 'Poor - Low expression may limit validation'
        
        # Drug delivery considerations
        if expression_data['Blood_Expression'] == 'Yes':
            assessment['Step7_Drug_Delivery_Considerations'] = 'Systemic delivery feasible'
        elif expression_data['Brain_Expression'] == 'Yes':
            assessment['Step7_Drug_Delivery_Considerations'] = 'CNS delivery required (BBB penetration)'
        else:
            top_category = max(
                expression_data['Expression_By_Category'].items(),
                key=lambda x: x[1],
                default=('Unknown', 0)
            )[0]
            assessment['Step7_Drug_Delivery_Considerations'] = f'Targeted delivery to {top_category} may be needed'
        
        # Off-target risk
        if specificity == 'Highly Specific':
            assessment['Step7_Off_Target_Risk'] = 'Low - Tissue-specific expression'
        elif specificity == 'Moderately Specific':
            assessment['Step7_Off_Target_Risk'] = 'Moderate - Some off-target expression expected'
        else:
            assessment['Step7_Off_Target_Risk'] = 'High - Broad expression may cause off-target effects'
        
        return assessment
    
    def add_expression_data(self, step6_df: pd.DataFrame) -> pd.DataFrame:
        """Add comprehensive expression data to Step 6"""
        
        gene_col = self._find_gene_column(step6_df)
        unique_genes = step6_df[gene_col].unique().tolist()
        
        logger.info(f"\nQuerying GTEx for {len(unique_genes)} genes...")
        logger.info("="*80)
        
        all_data = []
        
        for idx, gene in enumerate(unique_genes, 1):
            logger.info(f"[{idx}/{len(unique_genes)}] {gene}")
            
            expression_data = self.get_comprehensive_expression(gene)
            scores = self.calculate_druggability_scores(expression_data)
            assessment = self.assess_druggability(expression_data, scores)
            
            # Prepare category expression as string
            category_expr_str = '; '.join([
                f"{cat}:{tpm:.1f}" 
                for cat, tpm in sorted(
                    expression_data['Expression_By_Category'].items(),
                    key=lambda x: x[1],
                    reverse=True
                )[:5]  # Top 5 categories
            ])
            
            combined = {
                gene_col: gene,
                'Step7_GTEx_Data_Available': expression_data['GTEx_Data_Available'],
                'Step7_Total_Tissues': expression_data['Total_Tissues'],
                'Step7_Tissues_With_Expression': expression_data['Tissues_With_Expression'],
                'Step7_Max_TPM': round(expression_data['Max_TPM'], 3),
                'Step7_Mean_TPM': round(expression_data['Mean_TPM'], 3),
                'Step7_Median_TPM': round(expression_data['Median_TPM'], 3),
                'Step7_Tissue_Specificity': expression_data['Tissue_Specificity'],
                'Step7_High_Expression_Tissues': expression_data['Tissues_High_Expression'],
                'Step7_Medium_Expression_Tissues': expression_data['Tissues_Medium_Expression'],
                'Step7_Low_Expression_Tissues': expression_data['Tissues_Low_Expression'],
                'Step7_Top_3_Tissues': '; '.join(expression_data['Top_Expressed_Tissues'][:3]),
                'Step7_Category_Expression': category_expr_str,
                'Step7_Brain_Expression': expression_data['Brain_Expression'],
                'Step7_Brain_Max_TPM': round(expression_data['Brain_Max_TPM'], 3),
                'Step7_Blood_Expression': expression_data['Blood_Expression'],
                'Step7_Blood_TPM': round(expression_data['Blood_TPM'], 3),
                'Step7_Accessible_Expression': expression_data['Accessible_Expression'],
                **scores,
                **assessment
            }
            
            all_data.append(combined)
            
            logger.info(f"    Overall Score: {scores['Step7_Overall_Expression_Score']}/10, " +
                       f"Quality: {assessment['Step7_Expression_Quality']}")
        
        expression_df = pd.DataFrame(all_data)
        combined_df = step6_df.merge(expression_df, on=gene_col, how='left')
        
        logger.info(f"\n? Added expression data to {len(combined_df)} rows")
        return combined_df


def main():
    STEP6_INPUT = 'Step6.csv'
    STEP7_OUTPUT = 'Step7.csv'
    
    try:
        logger.info("="*80)
        logger.info("STEP 7: COMPREHENSIVE TISSUE EXPRESSION ANALYSIS")
        logger.info("="*80)
        
        analyzer = ComprehensiveTissueExpressionAnalyzer()
        step6_df = analyzer.load_step6_data(STEP6_INPUT)
        step7_df = analyzer.add_expression_data(step6_df)
        
        logger.info(f"\nSaving to {STEP7_OUTPUT}...")
        step7_df.to_csv(STEP7_OUTPUT, index=False)
        logger.info(f"? Saved {len(step7_df)} rows")
        
        # Summary
        gene_col = analyzer._find_gene_column(step7_df)
        protein_level = step7_df.drop_duplicates(subset=[gene_col])
        
        logger.info("\n" + "="*80)
        logger.info("EXPRESSION ANALYSIS SUMMARY")
        logger.info("="*80)
        
        logger.info(f"\nData availability: {protein_level['Step7_GTEx_Data_Available'].sum()}/{len(protein_level)}")
        
        logger.info(f"\nTissue specificity distribution:")
        for spec, count in protein_level['Step7_Tissue_Specificity'].value_counts().items():
            logger.info(f"  {spec}: {count}")
        
        logger.info(f"\nExpression quality distribution:")
        for qual, count in protein_level['Step7_Expression_Quality'].value_counts().items():
            logger.info(f"  {qual}: {count}")
        
        # Top 10 by expression level
        logger.info(f"\nTop 10 genes by max expression level:")
        top_expr = protein_level.nlargest(10, 'Step7_Max_TPM')
        for idx, row in top_expr.iterrows():
            logger.info(f"  {row[gene_col]}: {row['Step7_Max_TPM']:.1f} TPM " +
                       f"({row['Step7_Tissue_Specificity']})")
        
        # Save comprehensive summary
        summary_cols = [
            gene_col, 'Step7_Max_TPM', 'Step7_Tissue_Specificity',
            'Step7_Tissues_With_Expression', 'Step7_Top_3_Tissues',
            'Step7_Overall_Expression_Score', 'Step7_Expression_Quality',
            'Step7_Expression_Summary', 'Step7_Target_Validation_Feasibility'
        ]
        
        available_cols = [col for col in summary_cols if col in protein_level.columns]
        summary = protein_level[available_cols].sort_values(
            'Step7_Overall_Expression_Score', ascending=False
        )
        
        summary.to_csv('Step7_protein_summary.csv', index=False)
        logger.info(f"\n? Summary saved to Step7_protein_summary.csv")
        
        logger.info("\n" + "="*80)
        logger.info("COMPLETE!")
        logger.info("="*80)
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()