"""
Druggability Assessment Pipeline - Step 1: Comprehensive Gene ID Conversion
============================================================================
OUTPUT: Step1.csv with all gene identifiers
"""

import requests
import pandas as pd
import time
from typing import List, Dict, Optional
import json
from pathlib import Path


class ComprehensiveGeneIDConverter:
    """Get ALL available gene and protein identifiers."""
    
    def __init__(self):
        self.mygene_url = "https://mygene.info/v3"
        self.uniprot_url = "https://rest.uniprot.org/uniprotkb"
        self.ensembl_url = "https://rest.ensembl.org"
        
    def query_mygene_comprehensive(self, gene_symbols: List[str], species: str = "human") -> Dict:
        """Query MyGene.info with ALL available fields."""
        print(f"\nQuerying MyGene.info for {len(gene_symbols)} genes (comprehensive mode)...")
        
        url = f"{self.mygene_url}/query"
        
        # Request ALL available fields
        fields = [
            'ensembl.gene',
            'ensembl.transcript',
            'ensembl.protein',
            'ensembl.type_of_gene',
            'uniprot',
            'entrezgene',
            'refseq',
            'accession',
            'symbol',
            'name',
            'summary',
            'alias',
            'other_names',
            'type_of_gene',
            'genomic_pos',
            'HGNC',
            'MIM',
            'pharmgkb',
            'reagent',
            'interpro',
            'pfam',
            'prosite',
            'pdb',
            'ec',
            'pathway',
            'generif',
            'homologene',
            'taxid',
            'wikipedia'
        ]
        
        query_string = ','.join(gene_symbols)
        
        params = {
            'q': query_string,
            'scopes': 'symbol,alias',
            'fields': ','.join(fields),
            'species': species,
            'size': 1
        }
        
        try:
            response = requests.post(url, data=params, timeout=60)
            response.raise_for_status()
            data = response.json()
            
            results = {}
            for item in data:
                query = item.get('query', '')
                if query:
                    results[query] = item
            
            print(f"  ✓ Retrieved comprehensive data for {len(results)} genes")
            return results
            
        except Exception as e:
            print(f"  ✗ Error querying MyGene: {e}")
            return {}
    
    def query_uniprot_comprehensive(self, gene_symbol: str, organism: str = "9606") -> Optional[Dict]:
        """Query UniProt API for comprehensive protein information."""
        url = f"{self.uniprot_url}/search"
        
        params = {
            'query': f'gene:{gene_symbol} AND organism_id:{organism} AND reviewed:true',
            'format': 'json',
            'size': 1,
            'fields': 'accession,id,gene_names,protein_name,organism_name,length,xref_pdb,xref_refseq,xref_ensembl,xref_string,xref_geneid,cc_subcellular_location,ft_domain'
        }
        
        try:
            response = requests.get(url, params=params, timeout=15)
            response.raise_for_status()
            data = response.json()
            
            if data.get('results'):
                result = data['results'][0]
                return result
            
            # If no reviewed entry, try any entry
            params['query'] = f'gene:{gene_symbol} AND organism_id:{organism}'
            response = requests.get(url, params=params, timeout=15)
            data = response.json()
            
            if data.get('results'):
                return data['results'][0]
                
        except Exception as e:
            pass
        
        return None
    
    def extract_comprehensive_info(self, gene_symbol: str, mygene_data: Dict, 
                                   uniprot_data: Optional[Dict] = None) -> Dict:
        """Extract ALL available identifiers and information."""
        
        info = {
            # Basic identifiers
            'gene_symbol': gene_symbol,
            'gene_symbol_official': None,
            'gene_name': None,
            'gene_description': None,
            'gene_type': None,
            'aliases': None,
            'other_names': None,
            
            # Ensembl IDs
            'ensembl_gene_id': None,
            'ensembl_transcript_ids': None,
            'ensembl_protein_ids': None,
            
            # UniProt IDs
            'uniprot_swissprot': None,
            'uniprot_trembl': None,
            'uniprot_entry_name': None,
            
            # NCBI/Entrez
            'entrez_gene_id': None,
            
            # RefSeq IDs
            'refseq_mrna': None,
            'refseq_protein': None,
            'refseq_genomic': None,
            
            # Other database IDs
            'hgnc_id': None,
            'omim_id': None,
            'pharmgkb_id': None,
            'string_id': None,
            
            # Structure databases
            'pdb_ids': None,
            
            # Functional annotations
            'ec_number': None,
            'interpro_domains': None,
            'pfam_domains': None,
            'prosite_domains': None,
            
            # Genomic location
            'chromosome': None,
            'start_position': None,
            'end_position': None,
            'strand': None,
            
            # Protein information
            'protein_length': None,
            'subcellular_location': None,
            
            # Additional info
            'taxonomy_id': None,
            'homologene_id': None,
            'wikipedia_url': None,
            'pathways': None
        }
        
        if not mygene_data or 'notfound' in mygene_data:
            return info
        
        # ===== EXTRACT FROM MYGENE =====
        
        # Official symbol and names
        info['gene_symbol_official'] = mygene_data.get('symbol', gene_symbol)
        info['gene_name'] = mygene_data.get('name')
        info['gene_description'] = mygene_data.get('summary')
        info['gene_type'] = mygene_data.get('type_of_gene')
        
        # Aliases
        if 'alias' in mygene_data:
            aliases = mygene_data['alias']
            if isinstance(aliases, list):
                info['aliases'] = '|'.join(str(a) for a in aliases)
            else:
                info['aliases'] = str(aliases)
        
        if 'other_names' in mygene_data:
            other = mygene_data['other_names']
            if isinstance(other, list):
                info['other_names'] = '|'.join(str(o) for o in other[:5])
            else:
                info['other_names'] = str(other)
        
        # Ensembl IDs
        if 'ensembl' in mygene_data:
            ensembl = mygene_data['ensembl']
            if isinstance(ensembl, list):
                ensembl = ensembl[0] if ensembl else {}
            
            if isinstance(ensembl, dict):
                info['ensembl_gene_id'] = ensembl.get('gene')
                
                # Transcripts
                transcripts = ensembl.get('transcript')
                if transcripts:
                    if isinstance(transcripts, list):
                        info['ensembl_transcript_ids'] = '|'.join(transcripts[:10])
                    else:
                        info['ensembl_transcript_ids'] = str(transcripts)
                
                # Proteins
                proteins = ensembl.get('protein')
                if proteins:
                    if isinstance(proteins, list):
                        info['ensembl_protein_ids'] = '|'.join(proteins[:10])
                    else:
                        info['ensembl_protein_ids'] = str(proteins)
        
        # UniProt IDs
        if 'uniprot' in mygene_data:
            uniprot = mygene_data['uniprot']
            if isinstance(uniprot, dict):
                # Swiss-Prot (reviewed)
                swiss = uniprot.get('Swiss-Prot')
                if swiss:
                    info['uniprot_swissprot'] = swiss if isinstance(swiss, str) else (swiss[0] if swiss else None)
                
                # TrEMBL (unreviewed)
                trembl = uniprot.get('TrEMBL')
                if trembl:
                    if isinstance(trembl, list):
                        info['uniprot_trembl'] = '|'.join(trembl[:5])
                    else:
                        info['uniprot_trembl'] = str(trembl)
            elif isinstance(uniprot, str):
                info['uniprot_swissprot'] = uniprot
        
        # Entrez Gene ID
        info['entrez_gene_id'] = mygene_data.get('entrezgene')
        
        # RefSeq IDs
        if 'refseq' in mygene_data:
            refseq = mygene_data['refseq']
            if isinstance(refseq, dict):
                # mRNA
                rna = refseq.get('rna')
                if rna:
                    if isinstance(rna, list):
                        info['refseq_mrna'] = '|'.join(rna[:10])
                    else:
                        info['refseq_mrna'] = str(rna)
                
                # Protein
                protein = refseq.get('protein')
                if protein:
                    if isinstance(protein, list):
                        info['refseq_protein'] = '|'.join(protein[:10])
                    else:
                        info['refseq_protein'] = str(protein)
                
                # Genomic
                genomic = refseq.get('genomic')
                if genomic:
                    if isinstance(genomic, list):
                        info['refseq_genomic'] = '|'.join(genomic[:5])
                    else:
                        info['refseq_genomic'] = str(genomic)
        
        # HGNC ID
        if 'HGNC' in mygene_data:
            hgnc = mygene_data['HGNC']
            if isinstance(hgnc, list):
                info['hgnc_id'] = hgnc[0] if hgnc else None
            else:
                info['hgnc_id'] = str(hgnc)
        
        # OMIM ID
        if 'MIM' in mygene_data:
            mim = mygene_data['MIM']
            if isinstance(mim, list):
                info['omim_id'] = '|'.join(str(m) for m in mim)
            else:
                info['omim_id'] = str(mim)
        
        # PharmGKB ID
        if 'pharmgkb' in mygene_data:
            info['pharmgkb_id'] = str(mygene_data['pharmgkb'])
        
        # PDB IDs
        if 'pdb' in mygene_data:
            pdb = mygene_data['pdb']
            if isinstance(pdb, list):
                info['pdb_ids'] = '|'.join(pdb[:20])
            else:
                info['pdb_ids'] = str(pdb)
        
        # EC number
        if 'ec' in mygene_data:
            ec = mygene_data['ec']
            if isinstance(ec, list):
                info['ec_number'] = '|'.join(ec)
            else:
                info['ec_number'] = str(ec)
        
        # InterPro domains
        if 'interpro' in mygene_data:
            interpro = mygene_data['interpro']
            if isinstance(interpro, list):
                domains = [d.get('id', '') for d in interpro if isinstance(d, dict)]
                info['interpro_domains'] = '|'.join(domains[:10])
            elif isinstance(interpro, dict):
                info['interpro_domains'] = interpro.get('id')
        
        # Pfam domains
        if 'pfam' in mygene_data:
            pfam = mygene_data['pfam']
            if isinstance(pfam, list):
                domains = [p.get('id', '') if isinstance(p, dict) else str(p) for p in pfam]
                info['pfam_domains'] = '|'.join(domains[:10])
            else:
                info['pfam_domains'] = str(pfam)
        
        # PROSITE domains
        if 'prosite' in mygene_data:
            prosite = mygene_data['prosite']
            if isinstance(prosite, list):
                info['prosite_domains'] = '|'.join(prosite[:10])
            else:
                info['prosite_domains'] = str(prosite)
        
        # Genomic position
        if 'genomic_pos' in mygene_data:
            gpos = mygene_data['genomic_pos']
            if isinstance(gpos, dict):
                info['chromosome'] = gpos.get('chr')
                info['start_position'] = gpos.get('start')
                info['end_position'] = gpos.get('end')
                info['strand'] = gpos.get('strand')
            elif isinstance(gpos, list) and gpos:
                gpos = gpos[0]
                info['chromosome'] = gpos.get('chr')
                info['start_position'] = gpos.get('start')
                info['end_position'] = gpos.get('end')
                info['strand'] = gpos.get('strand')
        
        # Taxonomy
        info['taxonomy_id'] = mygene_data.get('taxid')
        
        # Homologene
        if 'homologene' in mygene_data:
            homolog = mygene_data['homologene']
            if isinstance(homolog, dict):
                info['homologene_id'] = homolog.get('id')
            else:
                info['homologene_id'] = str(homolog)
        
        # Wikipedia
        if 'wikipedia' in mygene_data:
            wiki = mygene_data['wikipedia']
            if isinstance(wiki, dict):
                info['wikipedia_url'] = wiki.get('url_stub')
        
        # Pathways
        if 'pathway' in mygene_data:
            pathway = mygene_data['pathway']
            if isinstance(pathway, dict):
                pathways = []
                for db, data in pathway.items():
                    if isinstance(data, list):
                        for item in data[:5]:
                            if isinstance(item, dict):
                                pathways.append(f"{db}:{item.get('id', '')}")
                if pathways:
                    info['pathways'] = '|'.join(pathways)
        
        # ===== EXTRACT FROM UNIPROT =====
        if uniprot_data:
            # Primary accession
            if not info['uniprot_swissprot']:
                info['uniprot_swissprot'] = uniprot_data.get('primaryAccession')
            
            # Entry name
            info['uniprot_entry_name'] = uniprot_data.get('uniProtkbId')
            
            # Protein length
            if 'sequence' in uniprot_data:
                info['protein_length'] = uniprot_data['sequence'].get('length')
            
            # Subcellular location
            if 'comments' in uniprot_data:
                for comment in uniprot_data['comments']:
                    if comment.get('commentType') == 'SUBCELLULAR LOCATION':
                        locations = comment.get('subcellularLocations', [])
                        if locations:
                            locs = []
                            for loc in locations[:3]:
                                loc_val = loc.get('location', {}).get('value')
                                if loc_val:
                                    locs.append(loc_val)
                            info['subcellular_location'] = '|'.join(locs)
                        break
            
            # PDB IDs from UniProt
            if 'uniProtKBCrossReferences' in uniprot_data and not info['pdb_ids']:
                pdb_refs = []
                for xref in uniprot_data['uniProtKBCrossReferences']:
                    if xref.get('database') == 'PDB':
                        pdb_refs.append(xref.get('id'))
                if pdb_refs:
                    info['pdb_ids'] = '|'.join(pdb_refs[:20])
            
            # STRING ID
            if 'uniProtKBCrossReferences' in uniprot_data:
                for xref in uniprot_data['uniProtKBCrossReferences']:
                    if xref.get('database') == 'STRING':
                        info['string_id'] = xref.get('id')
                        break
        
        return info
    
    def convert_genes_comprehensive(self, 
                                   gene_symbols: List[str],
                                   enrich_uniprot: bool = True,
                                   species: str = "human") -> pd.DataFrame:
        """
        Convert genes and get ALL available identifiers.
        
        Args:
            gene_symbols: List of gene symbols
            enrich_uniprot: Query UniProt for additional info (slower)
            species: Species name
            
        Returns:
            DataFrame with comprehensive gene information
        """
        print("="*70)
        print("COMPREHENSIVE GENE ID CONVERSION - STEP 1")
        print("="*70)
        print(f"\nInput: {len(gene_symbols)} gene symbols")
        print(f"Species: {species}")
        print(f"Enrich with UniProt: {enrich_uniprot}")
        
        # Query MyGene.info
        mygene_results = self.query_mygene_comprehensive(gene_symbols, species=species)
        
        # Process each gene
        all_results = []
        
        print(f"\nProcessing genes:")
        for idx, gene_symbol in enumerate(gene_symbols, 1):
            print(f"  [{idx}/{len(gene_symbols)}] {gene_symbol}...", end='')
            
            # Get MyGene data
            mygene_data = mygene_results.get(gene_symbol, {})
            
            # Enrich with UniProt if requested
            uniprot_data = None
            if enrich_uniprot:
                uniprot_data = self.query_uniprot_comprehensive(gene_symbol)
                time.sleep(0.3)
            
            # Extract all information
            info = self.extract_comprehensive_info(gene_symbol, mygene_data, uniprot_data)
            all_results.append(info)
            
            print(" ✓")
        
        # Create DataFrame
        df = pd.DataFrame(all_results)
        
        # Print summary
        self.print_summary(df)
        
        return df
    
    def print_summary(self, df: pd.DataFrame):
        """Print comprehensive summary statistics."""
        print("\n" + "="*70)
        print("COMPREHENSIVE CONVERSION SUMMARY")
        print("="*70)
        print(f"Total genes: {len(df)}")
        
        print("\n📊 Primary Identifiers:")
        print(f"  ✓ Ensembl Gene ID: {df['ensembl_gene_id'].notna().sum()} ({df['ensembl_gene_id'].notna().sum()/len(df)*100:.1f}%)")
        print(f"  ✓ UniProt (Swiss-Prot): {df['uniprot_swissprot'].notna().sum()} ({df['uniprot_swissprot'].notna().sum()/len(df)*100:.1f}%)")
        print(f"  ✓ Entrez Gene ID: {df['entrez_gene_id'].notna().sum()} ({df['entrez_gene_id'].notna().sum()/len(df)*100:.1f}%)")
        
        print("\n📊 Transcript & Protein IDs:")
        print(f"  ✓ Ensembl Transcripts: {df['ensembl_transcript_ids'].notna().sum()} ({df['ensembl_transcript_ids'].notna().sum()/len(df)*100:.1f}%)")
        print(f"  ✓ Ensembl Proteins: {df['ensembl_protein_ids'].notna().sum()} ({df['ensembl_protein_ids'].notna().sum()/len(df)*100:.1f}%)")
        print(f"  ✓ RefSeq mRNA: {df['refseq_mrna'].notna().sum()} ({df['refseq_mrna'].notna().sum()/len(df)*100:.1f}%)")
        print(f"  ✓ RefSeq Protein: {df['refseq_protein'].notna().sum()} ({df['refseq_protein'].notna().sum()/len(df)*100:.1f}%)")
        
        print("\n📊 Database Cross-references:")
        print(f"  ✓ HGNC ID: {df['hgnc_id'].notna().sum()} ({df['hgnc_id'].notna().sum()/len(df)*100:.1f}%)")
        print(f"  ✓ OMIM ID: {df['omim_id'].notna().sum()} ({df['omim_id'].notna().sum()/len(df)*100:.1f}%)")
        print(f"  ✓ PDB Structures: {df['pdb_ids'].notna().sum()} ({df['pdb_ids'].notna().sum()/len(df)*100:.1f}%)")
        
        print("\n📊 Functional Annotations:")
        print(f"  ✓ InterPro Domains: {df['interpro_domains'].notna().sum()} ({df['interpro_domains'].notna().sum()/len(df)*100:.1f}%)")
        print(f"  ✓ Pfam Domains: {df['pfam_domains'].notna().sum()} ({df['pfam_domains'].notna().sum()/len(df)*100:.1f}%)")
        print(f"  ✓ EC Number: {df['ec_number'].notna().sum()} ({df['ec_number'].notna().sum()/len(df)*100:.1f}%)")
        print(f"  ✓ Gene Description: {df['gene_description'].notna().sum()} ({df['gene_description'].notna().sum()/len(df)*100:.1f}%)")
        
        # Genes with no IDs
        missing_all = df[(df['ensembl_gene_id'].isna()) & 
                        (df['uniprot_swissprot'].isna()) & 
                        (df['entrez_gene_id'].isna())]['gene_symbol'].tolist()
        
        if missing_all:
            print(f"\n⚠️  Genes with NO IDs found: {', '.join(missing_all)}")
        
        print("="*70)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    # Load gene list from JSON file (created by main.py)
    gene_list_file = Path("gene_list.json")
    
    if not gene_list_file.exists():
        print("❌ Error: gene_list.json not found!")
        print("This file should be created by main.py")
        print("Please run main.py instead of running step01 directly")
        exit(1)
    
    try:
        with open(gene_list_file, 'r') as f:
            data = json.load(f)
            gene_list = data.get('genes', [])
        
        if not gene_list:
            print("❌ Error: No genes found in gene_list.json")
            exit(1)
        
        print(f"✓ Loaded {len(gene_list)} genes from gene_list.json")
    
    except Exception as e:
        print(f"❌ Error reading gene_list.json: {e}")
        exit(1)

    # Create converter
    converter = ComprehensiveGeneIDConverter()
    
    # Convert with full enrichment
    print("\n🔄 Running comprehensive conversion with UniProt enrichment...")
    results_df = converter.convert_genes_comprehensive(
        gene_symbols=gene_list,
        enrich_uniprot=True,
        species="human"
    )
    
    # Save as Step1.csv
    output_path = Path("Step1.csv")
    results_df.to_csv(output_path, index=False)
    print(f"\n✅ Results saved to: {output_path.absolute()}")
    
    # Display sample
    print("\n📋 Sample Results (first 5 genes, key columns):")
    key_columns = ['gene_symbol', 'ensembl_gene_id', 'uniprot_swissprot', 
                   'entrez_gene_id', 'hgnc_id', 'pdb_ids']
    print(results_df[key_columns].head().to_string())
    
    print(f"\n✅ Step 1 Complete! All gene IDs saved to Step1.csv")
    print(f"\nTotal columns: {len(results_df.columns)}")
    print(f"Total rows: {len(results_df)}")