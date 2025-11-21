"""
Druggability Assessment Pipeline - Step 5: Binding Pocket Prediction with Fpocket
==================================================================================

This script takes Step4.csv and predicts binding pockets using Fpocket:
- Runs fpocket on all protein structures
- Parses pocket information from {PDB_ID}_info.txt files
- Calculates pocket properties (volume, hydrophobicity, druggability)
- Saves comprehensive results

OUTPUT: Step5.csv with all Step4 fields + binding pocket data
"""

import pandas as pd
import subprocess
import os
import re
from typing import List, Dict, Optional
from pathlib import Path
import shutil


class FpocketBindingPocketPredictor:
    """Predict and analyze binding pockets using Fpocket."""
    
    def __init__(self, pockets_dir: str = "binding_pockets", structures_base: str = "protein_structures"):
        self.structures_base = Path(structures_base)
        self.pdb_dir = self.structures_base / "pdb_experimental"
        self.alphafold_dir = self.structures_base / "alphafold_predicted"
        self.pockets_dir = Path(pockets_dir)
        self.pockets_dir.mkdir(parents=True, exist_ok=True)
        
        self.fpocket_available = self.check_fpocket_installation()
        
        self.thresholds = {
            'volume_min': 200,
            'volume_optimal': 500,
            'druggability_min': 0.5,
            'druggability_good': 0.7
        }
    
    def check_fpocket_installation(self) -> bool:
        """Check if fpocket is installed."""
        try:
            result = subprocess.run(['fpocket', '-h'], capture_output=True, text=True, timeout=5)
            return result.returncode in [0, 1]
        except:
            return False
    
    def list_available_structures(self) -> Dict[str, Path]:
        """List all available structure files."""
        structure_files = {}
        
        if self.pdb_dir.exists():
            for pdb_file in self.pdb_dir.glob("*.pdb"):
                pdb_id = pdb_file.stem
                structure_files[pdb_id] = pdb_file
        
        if self.alphafold_dir.exists():
            for af_file in self.alphafold_dir.glob("*.pdb"):
                match = re.search(r'AF-([A-Z0-9]+)-', af_file.name)
                if match:
                    structure_files[match.group(1)] = af_file
        
        return structure_files
    
    def find_structure_for_gene(self, row: pd.Series, available_structures: Dict[str, Path]) -> Optional[Path]:
        """Find structure file for a gene."""
        # Try PDB IDs
        pdb_ids = row.get('pdb_ids')
        if pdb_ids and not pd.isna(pdb_ids):
            pdb_list = str(pdb_ids).split('|')
            for pdb_id in pdb_list:
                pdb_id = pdb_id.strip()
                if pdb_id in available_structures:
                    return available_structures[pdb_id]
        
        # Try best_pdb_id
        best_pdb = row.get('best_pdb_id')
        if best_pdb and not pd.isna(best_pdb):
            if best_pdb in available_structures:
                return available_structures[best_pdb]
        
        # Try UniProt for AlphaFold
        uniprot_id = row.get('uniprot_swissprot')
        if uniprot_id and not pd.isna(uniprot_id):
            if uniprot_id in available_structures:
                return available_structures[uniprot_id]
        
        return None
    
    def run_fpocket(self, pdb_file: Path, gene_symbol: str) -> Optional[Path]:
        """
        Run fpocket on a PDB file.
        Fpocket creates output directory next to the input file.
        """
        if not pdb_file.exists():
            return None
        
        # Run fpocket from the PDB file's directory
        cmd = ['fpocket', '-f', pdb_file.name]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(pdb_file.parent),
                timeout=120
            )
            
            # Fpocket creates {filename}_out directory
            output_dir = pdb_file.parent / f"{pdb_file.stem}_out"
            
            if output_dir.exists():
                # Move to organized location
                target_dir = self.pockets_dir / gene_symbol
                target_dir.mkdir(parents=True, exist_ok=True)
                
                final_output = target_dir / f"{pdb_file.stem}_out"
                
                # Remove old if exists
                if final_output.exists():
                    shutil.rmtree(final_output)
                
                # Move output directory
                shutil.move(str(output_dir), str(final_output))
                
                return final_output
            
        except Exception as e:
            print(f"Error: {str(e)[:50]}")
        
        return None
    
    def parse_fpocket_info_file(self, fpocket_dir: Path) -> List[Dict]:
        """
        Parse fpocket info file.
        File location: {PDB_ID}_out/{PDB_ID}_info.txt
        """
        pockets = []
        
        # Find info file - format: {PDB_ID}_info.txt
        info_files = list(fpocket_dir.glob("*_info.txt"))
        
        if not info_files:
            return pockets
        
        info_file = info_files[0]
        
        with open(info_file, 'r') as f:
            content = f.read()
        
        # Parse each pocket section
        sections = re.split(r'Pocket\s+(\d+)\s*:', content)
        
        for i in range(1, len(sections), 2):
            pocket_num = int(sections[i])
            pocket_text = sections[i + 1]
            
            pocket = {
                'pocket_id': pocket_num,
                'score': self._extract_value(pocket_text, r'Score\s*:\s*([\d.]+)'),
                'druggability_score': self._extract_value(pocket_text, r'Druggability Score\s*:\s*([\d.]+)'),
                'volume': self._extract_value(pocket_text, r'Volume\s*:\s*([\d.]+)'),
                'hydrophobicity': self._extract_value(pocket_text, r'Hydrophobicity score\s*:\s*([\d.]+)'),
                'polarity': self._extract_value(pocket_text, r'Polarity score\s*:\s*([\d.]+)'),
                'alpha_spheres': self._extract_value(pocket_text, r'Number of Alpha Spheres\s*:\s*(\d+)', is_int=True),
            }
            
            pockets.append(pocket)
        
        return pockets
    
    def _extract_value(self, text: str, pattern: str, is_int: bool = False) -> Optional[float]:
        """Extract numeric value from text using regex."""
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            try:
                value = match.group(1)
                return int(value) if is_int else float(value)
            except:
                pass
        return None
    
    def assess_pocket_druggability(self, pocket: Dict) -> Dict:
        """Assess pocket druggability based on multiple criteria."""
        assessment = {
            'is_druggable': False,
            'category': 'Poor',
            'confidence': 'Low',
            'favorable_features': [],
            'concerns': []
        }
        
        score_components = []
        
        # Volume assessment
        volume = pocket.get('volume')
        if volume:
            if volume >= self.thresholds['volume_optimal']:
                score_components.append(1.0)
                assessment['favorable_features'].append(f"Optimal volume ({volume:.0f}A3)")
            elif volume >= self.thresholds['volume_min']:
                score_components.append(0.7)
                assessment['favorable_features'].append(f"Adequate volume ({volume:.0f}A3)")
            else:
                score_components.append(0.3)
                assessment['concerns'].append(f"Small volume ({volume:.0f}A3)")
        
        # Druggability score
        drug_score = pocket.get('druggability_score')
        if drug_score:
            if drug_score >= self.thresholds['druggability_good']:
                score_components.append(1.0)
                assessment['favorable_features'].append(f"High druggability ({drug_score:.2f})")
            elif drug_score >= self.thresholds['druggability_min']:
                score_components.append(0.6)
            else:
                score_components.append(0.3)
                assessment['concerns'].append(f"Low druggability ({drug_score:.2f})")
        
        # Hydrophobicity
        hydro = pocket.get('hydrophobicity')
        if hydro:
            if hydro > 0.5:
                score_components.append(0.8)
            elif hydro > 0.3:
                score_components.append(0.5)
        
        # Overall score
        overall = pocket.get('score')
        if overall:
            if overall > 0.5:
                score_components.append(0.7)
            else:
                score_components.append(0.4)
        
        # Calculate final assessment
        if score_components:
            final_score = sum(score_components) / len(score_components)
            
            if final_score >= 0.8:
                assessment['is_druggable'] = True
                assessment['category'] = 'Excellent'
                assessment['confidence'] = 'High'
            elif final_score >= 0.65:
                assessment['is_druggable'] = True
                assessment['category'] = 'Good'
                assessment['confidence'] = 'Medium'
            elif final_score >= 0.5:
                assessment['is_druggable'] = True
                assessment['category'] = 'Moderate'
                assessment['confidence'] = 'Low'
        
        return assessment
    
    def process_gene_row(self, row: pd.Series, available_structures: Dict[str, Path]) -> Dict:
        """Process a single gene and predict binding pockets."""
        gene_symbol = row['gene_symbol']
        
        print(f"  {gene_symbol}...", end=' ')
        
        result = {
            'gene_symbol': gene_symbol,
            'fpocket_run': False,
            'pockets_found': 0,
            'best_pocket_id': None,
            'best_pocket_volume': None,
            'best_pocket_druggability': None,
            'best_pocket_hydrophobicity': None,
            'best_pocket_score': None,
            'has_druggable_pocket': False,
            'num_druggable_pockets': 0,
            'best_pocket_category': None,
            'best_pocket_confidence': None,
            'favorable_features': None,
            'concerns': None,
            'all_pockets_volumes': None,
            'all_pockets_druggability': None,
            'pocket_details': None
        }
        
        # Find structure file
        structure_file = self.find_structure_for_gene(row, available_structures)
        
        if not structure_file:
            print("X No structure")
            return result
        
        print(f"{structure_file.name}...", end=' ')
        
        # Run fpocket
        fpocket_output = self.run_fpocket(structure_file, gene_symbol)
        
        if not fpocket_output:
            print("X Failed")
            return result
        
        result['fpocket_run'] = True
        
        # Parse results
        pockets = self.parse_fpocket_info_file(fpocket_output)
        
        if not pockets:
            print("X No pockets")
            return result
        
        result['pockets_found'] = len(pockets)
        
        # Analyze pockets
        druggable_pockets = []
        volumes = []
        drug_scores = []
        details = []
        
        for pocket in pockets:
            assessment = self.assess_pocket_druggability(pocket)
            
            if assessment['is_druggable']:
                druggable_pockets.append(pocket)
            
            vol = pocket.get('volume')
            if vol:
                volumes.append(f"{vol:.0f}")
            
            drug = pocket.get('druggability_score')
            if drug:
                drug_scores.append(f"{drug:.2f}")
            
            details.append(
                f"P{pocket['pocket_id']}:Vol={vol if vol else 0:.0f}:Drug={drug if drug else 0:.2f}"
            )
        
        result['num_druggable_pockets'] = len(druggable_pockets)
        result['has_druggable_pocket'] = len(druggable_pockets) > 0
        
        if volumes:
            result['all_pockets_volumes'] = '|'.join(volumes)
        if drug_scores:
            result['all_pockets_druggability'] = '|'.join(drug_scores)
        if details:
            result['pocket_details'] = '|'.join(details[:10])
        
        # Get best pocket
        if pockets:
            best = max(pockets, key=lambda p: p.get('druggability_score') or 0)
            best_assessment = self.assess_pocket_druggability(best)
            
            result['best_pocket_id'] = best['pocket_id']
            result['best_pocket_volume'] = best.get('volume')
            result['best_pocket_druggability'] = best.get('druggability_score')
            result['best_pocket_hydrophobicity'] = best.get('hydrophobicity')
            result['best_pocket_score'] = best.get('score')
            result['best_pocket_category'] = best_assessment['category']
            result['best_pocket_confidence'] = best_assessment['confidence']
            
            if best_assessment['favorable_features']:
                result['favorable_features'] = '|'.join(best_assessment['favorable_features'])
            if best_assessment['concerns']:
                result['concerns'] = '|'.join(best_assessment['concerns'])
            
            print(f"OK {len(pockets)}p {len(druggable_pockets)}d "
                  f"Vol:{best.get('volume',0):.0f} Drug:{best.get('druggability_score',0):.2f} "
                  f"{best_assessment['category']}")
        
        return result
    
    def process_step4_file(self, input_file: str = "Step4.csv") -> pd.DataFrame:
        """Process Step4.csv and add binding pocket predictions."""
        print("="*70)
        print("BINDING POCKET PREDICTION - STEP 5")
        print("="*70)
        
        if not self.fpocket_available:
            print("\nERROR: Fpocket not installed")
            print("Install: sudo apt-get install fpocket")
            return None
        
        print("\nOK Fpocket available")
        
        # List available structures
        print("\nScanning for structures...")
        available_structures = self.list_available_structures()
        
        if self.pdb_dir.exists():
            pdb_count = len(list(self.pdb_dir.glob("*.pdb")))
            print(f"  PDB: {pdb_count} files")
        if self.alphafold_dir.exists():
            af_count = len(list(self.alphafold_dir.glob("*.pdb")))
            print(f"  AlphaFold: {af_count} files")
        
        # Read Step4.csv
        print(f"\nReading {input_file}...")
        step4_df = pd.read_csv(input_file)
        print(f"  OK Loaded {len(step4_df)} genes")
        
        # Process each gene
        print(f"\nRunning fpocket on all genes...\n")
        
        pocket_results = []
        
        for idx, row in step4_df.iterrows():
            print(f"[{idx+1}/{len(step4_df)}]", end=' ')
            result = self.process_gene_row(row, available_structures)
            pocket_results.append(result)
        
        # Create pocket DataFrame
        pocket_df = pd.DataFrame(pocket_results)
        
        # Merge with Step4 data
        pocket_df = pocket_df.drop(columns=['gene_symbol'])
        step5_df = pd.concat([step4_df, pocket_df], axis=1)
        
        # Print summary
        self.print_summary(step5_df)
        
        return step5_df
    
    def print_summary(self, df: pd.DataFrame):
        """Print comprehensive summary of pocket predictions."""
        print("\n" + "="*70)
        print("BINDING POCKET PREDICTION SUMMARY")
        print("="*70)
        
        total = len(df)
        fpocket_success = (df['fpocket_run'] == True).sum()
        pockets_found = (df['pockets_found'] > 0).sum()
        has_druggable = (df['has_druggable_pocket'] == True).sum()
        
        print(f"\nTotal genes: {total}")
        print(f"Fpocket success: {fpocket_success} ({fpocket_success/total*100:.1f}%)")
        print(f"Pockets found: {pockets_found} ({pockets_found/total*100:.1f}%)")
        print(f"Druggable pockets: {has_druggable} ({has_druggable/total*100:.1f}%)")
        
        if pockets_found > 0:
            total_pockets = df['pockets_found'].sum()
            avg_pockets = df[df['pockets_found'] > 0]['pockets_found'].mean()
            total_druggable = df['num_druggable_pockets'].sum()
            
            print(f"\nPocket Statistics:")
            print(f"  Total pockets: {int(total_pockets)}")
            print(f"  Avg pockets/protein: {avg_pockets:.1f}")
            print(f"  Total druggable: {int(total_druggable)}")
            
            if total_pockets > 0:
                print(f"  Druggability rate: {total_druggable/total_pockets*100:.1f}%")
        
        if has_druggable > 0:
            # Categories
            excellent = (df['best_pocket_category'] == 'Excellent').sum()
            good = (df['best_pocket_category'] == 'Good').sum()
            moderate = (df['best_pocket_category'] == 'Moderate').sum()
            
            print(f"\nDruggability Categories:")
            print(f"  Excellent: {excellent}")
            print(f"  Good: {good}")
            print(f"  Moderate: {moderate}")
            
            # Volume distribution
            volumes = df[df['best_pocket_volume'].notna()]['best_pocket_volume']
            if len(volumes) > 0:
                print(f"\nPocket Volumes:")
                print(f"  Average: {volumes.mean():.0f} A3")
                print(f"  Range: {volumes.min():.0f}-{volumes.max():.0f} A3")
            
            # Top 10
            print(f"\nTop 10 Druggable Genes:")
            top = df[df['has_druggable_pocket'] == True].sort_values(
                'best_pocket_druggability', ascending=False
            ).head(10)
            
            for _, row in top.iterrows():
                print(f"  {row['gene_symbol']:10s} P{row['best_pocket_id']}: "
                      f"Vol={row['best_pocket_volume']:.0f}A3 "
                      f"Drug={row['best_pocket_druggability']:.2f} "
                      f"{row['best_pocket_category']}")
        
        print("="*70)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    if not Path("Step4.csv").exists():
        print("ERROR: Step4.csv not found!")
        exit(1)
    
    predictor = FpocketBindingPocketPredictor(
        pockets_dir="binding_pockets",
        structures_base="protein_structures"
    )
    
    step5_df = predictor.process_step4_file("Step4.csv")
    
    if step5_df is not None:
        output_path = Path("Step5.csv")
        step5_df.to_csv(output_path, index=False)
        
        print(f"\nOK Results saved to: {output_path.absolute()}")
        print(f"Total columns: {len(step5_df.columns)}")
        print(f"Total rows: {len(step5_df)}")
        
        print("\nStep 5 Complete!")
        print(f"Pocket analysis files: {predictor.pockets_dir.absolute()}/")
    else:
        print("\nERROR: Step 5 Failed")