import pandas as pd
import requests

# 1. Load the CSV
try:
    df = pd.read_csv("results/phase1_targets.csv")
    print("✅ CSV Loaded Successfully.")
except Exception as e:
    print(f"❌ Error loading CSV: {e}")
    exit()

# 2. Grab the first ID
try:
    first_gene_raw = df.iloc[0]['Ensembl_ID']
    print(f"\n🔍 DEBUGGING GENE #1")
    print(f"   Raw ID in CSV:   '{first_gene_raw}'") # Quotes show us if there's hidden whitespace
except KeyError:
    print("❌ Error: Column 'Ensembl_ID' not found in CSV.")
    print(f"   Columns found: {df.columns.tolist()}")
    exit()

# 3. Test the Cleaning Logic
clean_id = str(first_gene_raw).split('.')[0]
print(f"   Cleaned ID:      '{clean_id}'")

# 4. Test the API Connection (MyGene.info)
url = f"https://mygene.info/v3/gene/{clean_id}?fields=uniprot"
print(f"   Querying URL:    {url}")

try:
    response = requests.get(url)
    print(f"   API Status:      {response.status_code}")
    print(f"   API Response:    {response.text}") # <--- This is the key clue
except Exception as e:
    print(f"❌ API Connection Failed: {e}")

# 5. Check AlphaFold URL (Hypothetical)
print(f"\n------------------------------------------------")
print("If the API response above is empty or missing 'Swiss-Prot', that is why it failed.")