import os
import pandas as pd
import numpy as np
from Bio import PDB
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from tqdm import tqdm
import warnings

# --- CONFIGURATION ---
INPUT_DIR = "data/structures_pdb"
OUTPUT_FILE = "results/structural_metrics.csv"

# Suppress annoying PDB warnings (common in BioPython)
warnings.filterwarnings('ignore')

# --- HELPER: EXTRACT pLDDT SCORE ---
def get_avg_plddt(pdb_path):
    """
    AlphaFold hides its confidence score (0-100) in the 'B-factor' column of the PDB file.
    High (>90) = Very confident.
    Low (<50) = Disordered/Unstructured.
    """
    parser = PDB.PDBParser()
    structure = parser.get_structure("temp", pdb_path)
    
    plddt_scores = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # The B-factor column holds the pLDDT score
                    plddt_scores.append(atom.get_bfactor())
    
    # Return the average score for the whole protein
    return np.mean(plddt_scores) if plddt_scores else 0

# --- HELPER: CALCULATE BIOCHEMICAL PROPERTIES ---
def get_biochem_props(pdb_path):
    """
    Extracts the amino acid sequence and calculates Stability & Hydrophobicity.
    """
    parser = PDB.PDBParser()
    structure = parser.get_structure("temp", pdb_path)
    
    # 1. Extract Sequence from 3D structure
    ppb = PDB.PPBuilder()
    peptides = ppb.build_peptides(structure)
    sequence = "".join([str(pp.get_sequence()) for pp in peptides])
    
    if not sequence:
        return None, None, 0 # Empty structure
        
    # 2. Analyze Sequence
    analyzed_seq = ProteinAnalysis(sequence)
    
    # Instability Index (>40 means unstable)
    try:
        instability = analyzed_seq.instability_index()
    except:
        instability = 0 # Short sequences fail this test
        
    # Aromaticity (How many ring structures? often relates to druggability)
    aromaticity = analyzed_seq.aromaticity()
    
    # Molecular Weight
    weight = analyzed_seq.molecular_weight()
    
    return instability, aromaticity, weight

# --- MAIN WORKFLOW ---
def main():
    print(f"📂 Analyzing structures in: {INPUT_DIR}")
    
    files = [f for f in os.listdir(INPUT_DIR) if f.endswith(".pdb")]
    results = []

    print(f"🚀 Processing {len(files)} files...")
    
    for filename in tqdm(files, unit="protein"):
        file_path = os.path.join(INPUT_DIR, filename)
        
        # Filename format is usually "SYMBOL_ID.pdb"
        # We split it to get the Gene Symbol back
        try:
            gene_symbol = filename.split('_')[0]
            uniprot_id = filename.split('_')[1].replace('.pdb', '')
        except:
            gene_symbol = filename
            uniprot_id = "Unknown"
            
        # 1. Get Structural Confidence (pLDDT)
        plddt = get_avg_plddt(file_path)
        
        # 2. Get Biochemical Properties
        instability, aromaticity, weight = get_biochem_props(file_path)
        
        results.append({
            "Gene_Symbol": gene_symbol,
            "UniProt_ID": uniprot_id,
            "pLDDT_Score": round(plddt, 2),
            "Instability_Index": round(instability, 2),
            "Aromaticity": round(aromaticity, 4),
            "Molecular_Weight": round(weight, 2)
        })

    # --- EXPORT ---
    df_results = pd.DataFrame(results)
    
    # Sort by Instability (Most Stable first)
    df_results = df_results.sort_values(by="Instability_Index")
    
    # Save
    df_results.to_csv(OUTPUT_FILE, index=False)
    
    print("-" * 60)
    print("✅ Analysis Complete!")
    print(f"📄 Results saved to: {OUTPUT_FILE}")
    print("\n🔍 TOP 5 MOST STABLE PROTEINS FOUND:")
    print(df_results[['Gene_Symbol', 'pLDDT_Score', 'Instability_Index']].head(5))

if __name__ == "__main__":
    main()