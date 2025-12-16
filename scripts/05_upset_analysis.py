import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_contents
import requests


PHASE1_FILE = "results/phase1_targets.csv"
PHASE3_FILE = "results/structural_metrics.csv"
OUTPUT_PLOT = "plots/upset_consequence.png"


def add_symbols_to_expression_data(df):
    ids = df['Ensembl_ID'].astype(str).tolist()
    clean_ids = [i.split('.')[0] for i in ids]
    
    url = "https://mygene.info/v3/query"
    params = {
        'q': ",".join(clean_ids[:1000]),
        'scopes': 'ensembl.gene',
        'fields': 'symbol',
        'species': 'human'
    }
    try:
        response = requests.post(url, data=params, headers={'content-type': 'application/x-www-form-urlencoded'})
        data = response.json()
        id_map = {}
        for item in data:
            if 'symbol' in item:
                id_map[item['query']] = item['symbol']
        
        df['clean_id'] = df['Ensembl_ID'].apply(lambda x: str(x).split('.')[0])
        df['Gene_Symbol'] = df['clean_id'].map(id_map)
    except:
        pass
    return df


def main():
    print("📊 Loading Data for UpSet Analysis...")
    exp_df = pd.read_csv(PHASE1_FILE)
    struct_df = pd.read_csv(PHASE3_FILE)

    exp_df = add_symbols_to_expression_data(exp_df)
    merged_df = pd.merge(struct_df, exp_df, on="Gene_Symbol", how="inner")
    
    # --- CATEGORICAL SETS (HIGH EXPRESSION, STABLE AND HIGH CONFIDENCE)
    high_exp = set(merged_df[merged_df['log2FoldChange'] > 2.0]['Gene_Symbol'])
    stable = set(merged_df[merged_df['Instability_Index'] < 40]['Gene_Symbol'])
    confident = set(merged_df[merged_df['pLDDT_Score'] > 70]['Gene_Symbol'])

    content = {
        'High Expression (>2.0 FC)': high_exp,
        'Stable Structure (<40 Index)': stable,
        'High Confidence (>70 pLDDT)': confident
    }
    
    upset_data = from_contents(content)

    print("🎨 Generating Publication-Ready UpSet Plot...")
    
    fig = plt.figure(figsize=(12, 6))
    
    upset = UpSet(upset_data, subset_size='count', show_counts=True, sort_by='cardinality', element_size=40)
    
    # Subsets identified with necessary annotations
    upset.style_subsets(present=["High Expression (>2.0 FC)", "Stable Structure (<40 Index)", "High Confidence (>70 pLDDT)"],
                        facecolor="green", label="Goldilocks Targets")
    
    plot_dict = upset.plot(fig=fig)
 
    plot_dict['intersections'].legend(
        loc='upper left', 
        bbox_to_anchor=(1.05, 1.0), 
        title="Target Classification",
        frameon=False 
    )
    
    plt.suptitle("Intersection of Druggability Criteria: Finding the 'Perfect' Targets", fontsize=16, y=0.98)
    
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
    print(f"🎉 Polished plot saved to: {OUTPUT_PLOT}")

if __name__ == "__main__":
    main()
