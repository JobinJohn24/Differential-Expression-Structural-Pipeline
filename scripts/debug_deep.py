import pandas as pd
import requests
import json

# --- CONFIG ---
INPUT_FILE = "results/phase1_targets.csv"

# --- DEBUG SCRIPT ---
print(f"🔍 DIAGNOSTIC MODE: Checking {INPUT_FILE}...\n")

# 1. Load Data
try:
    df = pd.read_csv(INPUT_FILE)
    print(f"✅ CSV Loaded. Columns found: {list(df.columns)}")
    print(f"✅ Total Rows: {len(df)}")
except Exception as e:
    print(f"❌ FATAL: Could not read CSV. {e}")
    exit()

# 2. Check the First 3 Rows ONLY
print("\n🔍 detailed Inspection of Top 3 Genes:")
print("="*60)

for i in range(3):
    row = df.iloc[i]
    raw_id = row['Ensembl_ID']
    clean_id = str(raw_id).split('.')[0]
    
    print(f"GENE #{i+1}")
    print(f"  - Raw ID in CSV:  '{raw_id}'")
    print(f"  - Cleaned ID:     '{clean_id}'")
    
    # 3. Test MyGene.info API manually
    url = f"https://mygene.info/v3/gene/{clean_id}?fields=symbol,uniprot"
    print(f"  - Requesting:     {url}")
    
    try:
        response = requests.get(url, timeout=10)
        print(f"  - API Status:     {response.status_code}")
        
        data = response.json()
        # DUMP THE RAW JSON so we can see what's wrong
        print(f"  - RAW JSON:       {json.dumps(data, indent=2)}")
        
        # 4. Test the Extraction Logic
        uniprot_id = None
        if 'uniprot' in data and 'Swiss-Prot' in data['uniprot']:
            result = data['uniprot']['Swiss-Prot']
            uniprot_id = result[0] if isinstance(result, list) else result
            print(f"  - ✅ EXTRACTED:    {uniprot_id}")
        else:
            print(f"  - ❌ FAILED extraction. 'Swiss-Prot' key missing?")
            
    except Exception as e:
        print(f"  - ❌ NETWORK ERROR: {e}")

    print("-" * 60)