import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import requests
from tqdm import tqdm


PHASE1_FILE = "results/phase1_targets.csv"       # Expression Data
PHASE3_FILE = "results/structural_metrics.csv"   # Structural Data
OUTPUT_PLOT = "plots/structural_volcano.png"


def add_symbols_to_expression_data(df):
    """
    We need to match the Ensembl IDs (Phase 1) to the Symbols (Phase 3).
    We ask MyGene.info to translate the Ensembl IDs for us.
    """
    print("🔄 Mapping Ensembl IDs to Gene Symbols (this links the two datasets)...")
    
    # Get list of IDs
    ids = df['Ensembl_ID'].astype(str).tolist()
    clean_ids = [i.split('.')[0] for i in ids]
    
    # Batch query MyGene
    url = "https://mygene.info/v3/query"
    headers = {'content-type': 'application/x-www-form-urlencoded'}
    params = {
        'q': ",".join(clean_ids[:1000]), 
        'scopes': 'ensembl.gene',
        'fields': 'symbol',
        'species': 'human'
    }
    
    response = requests.post(url, data=params, headers=headers)
    data = response.json()
    
    
    id_map = {}
    for item in data:
        if 'symbol' in item:
            id_map[item['query']] = item['symbol']
            
    # Apply map
    df['clean_id'] = df['Ensembl_ID'].apply(lambda x: str(x).split('.')[0])
    df['Gene_Symbol'] = df['clean_id'].map(id_map)
    
    return df


def main():
    print("📊 Loading Datasets...")
    
    
    exp_df = pd.read_csv(PHASE1_FILE)
    struct_df = pd.read_csv(PHASE3_FILE)
    
    print(f"   - Expression Data: {len(exp_df)} rows")
    print(f"   - Structural Data: {len(struct_df)} rows")

    
    exp_df = add_symbols_to_expression_data(exp_df)
    
    merged_df = pd.merge(struct_df, exp_df, on="Gene_Symbol", how="inner")
    
    print(f"✅ Successfully merged {len(merged_df)} targets for plotting.")
    
    
    sns.set_style("whitegrid")
    
    # Scatter Plot
    # X = Log2FoldChange (How much cancer loves it)
    # Y = Instability Index (Lower is Better/More Stable)
    # Color = pLDDT (Confidence)
    
    plot = sns.scatterplot(
        data=merged_df,
        x="log2FoldChange",
        y="Instability_Index",
        hue="pLDDT_Score",
        palette="viridis",
        s=100,
        edgecolor="black",
        alpha=0.8
    )
    
    
    plt.axhline(y=40, color='red', linestyle='--', label="Unstable (>40)")
    plt.axvline(x=2.0, color='blue', linestyle='--', label="High Expression (>2.0)")
    
    top_genes = merged_df.sort_values("Instability_Index").head(5)
    for i, row in top_genes.iterrows():
        plt.text(
            row['log2FoldChange'] + 0.1, 
            row['Instability_Index'], 
            row['Gene_Symbol'], 
            fontsize=10, 
            fontweight='bold'
        )

    plt.title("The Structural Consequence Map: LUAD Targets", fontsize=15)
    plt.xlabel("Gene Expression (Log2 Fold Change)", fontsize=12)
    plt.ylabel("Protein Instability Index (Lower is More Stable)", fontsize=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Save
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"🎉 Plot saved to: {OUTPUT_PLOT}")
    
    # Show the table of winners
    print("\n🏆 THE FINAL CANDIDATES (Stable & Highly Expressed):")
    print(merged_df[['Gene_Symbol', 'log2FoldChange', 'Instability_Index', 'pLDDT_Score']]
          .sort_values("Instability_Index")
          .head(10))

if __name__ == "__main__":
    main()
