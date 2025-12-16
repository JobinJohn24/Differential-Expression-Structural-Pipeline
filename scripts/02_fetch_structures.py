import pandas as pd
import requests
import os
from tqdm import tqdm


INPUT_FILE = "results/phase1_targets.csv"
OUTPUT_DIR = "data/structures_pdb"
TARGET_COUNT = 100  

os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- HELPER 1: Get UniProt ID from Ensembl ---
def get_uniprot_id(ensembl_id):
    clean_id = str(ensembl_id).split('.')[0]
    url = f"https://mygene.info/v3/gene/{clean_id}?fields=symbol,uniprot"
    try:
        response = requests.get(url, timeout=5)
        data = response.json()
        symbol = data.get('symbol', clean_id)
        
        uniprot_id = None
        if 'uniprot' in data and 'Swiss-Prot' in data['uniprot']:
            result = data['uniprot']['Swiss-Prot']
            uniprot_id = result[0] if isinstance(result, list) else result
        
        return symbol, uniprot_id
    except:
        return clean_id, None


def get_alphafold_url(uniprot_id):
    """
    Asks the AlphaFold API for the correct download link.
    """
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(api_url, timeout=5)
        data = response.json()
        if isinstance(data, list) and len(data) > 0:
            return data[0]['pdbUrl']
    except:
        return None
    return None


def main():
    print(f"📂 Reading targets from: {INPUT_FILE}")
    try:
        df = pd.read_csv(INPUT_FILE)
    except FileNotFoundError:
        print("❌ ERROR: CSV not found.")
        return

    print(f"🚀 Scanning for {TARGET_COUNT} structures using Smart API...")
    
    downloaded = 0
    pbar = tqdm(total=TARGET_COUNT, unit="pdb")

    for index, row in df.iterrows():
        if downloaded >= TARGET_COUNT:
            break
            
        gene_ens_id = row['Ensembl_ID']
        symbol, uniprot_id = get_uniprot_id(gene_ens_id)
        
        if uniprot_id:
            # 1. Ask API for the link
            pdb_url = get_alphafold_url(uniprot_id)
            
            if pdb_url:
                # 2. Download the specific link provided by the API
                try:
                    r = requests.get(pdb_url, timeout=10)
                    if r.status_code == 200:
                        save_path = os.path.join(OUTPUT_DIR, f"{symbol}_{uniprot_id}.pdb")
                        with open(save_path, "w") as f:
                            f.write(r.text)
                        
                        downloaded += 1
                        pbar.update(1)
                        pbar.set_description(f"✅ Found: {symbol}")
                except:
                    pass
            else:
                pbar.set_description(f"Scanning... (No AF Entry)")
        else:
            pbar.set_description(f"Scanning...")

    pbar.close()
    print("\n" + "="*40)
    print(f"🎉 Mission Complete! Downloaded: {downloaded} structures.")
    print(f"📂 Saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
