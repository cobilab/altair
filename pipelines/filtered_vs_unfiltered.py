import os
from collections import Counter
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import pearsonr, spearmanr

# Add this at the beginning of your script
genetic_code = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

# Define key and secondary mutations for major variants
variant_mutations = {
    'Ancestral': {'key': [], 'secondary': []},
    'D614G': {'key': ['D614G'], 'secondary': []},
    'Alpha': {'key': ['N501Y', 'A570D', 'D614G', 'P681H'], 'secondary': ['T716I', 'S982A', 'D1118H', '69-70del', 'Y144del']},
    'Beta': {'key': ['K417N', 'E484K', 'N501Y', 'D614G'], 'secondary': ['L18F', 'D80A', 'D215G', 'A701V', '242-244del']},
    'Gamma': {'key': ['K417T', 'E484K', 'N501Y', 'D614G'], 'secondary': ['L18F', 'T20N', 'P26S', 'D138Y', 'R190S', 'H655Y', 'T1027I', 'V1176F']},
    'Delta': {'key': ['L452R', 'T478K', 'D614G', 'P681R'], 'secondary': ['T19R', 'G142D', 'R158G', '156-157del', 'F157del', 'R203M', 'D950N']},
    'Omicron BA.1': {'key': ['N501Y', 'T478K', 'D614G', 'P681H', 'K417N'], 'secondary': ['G339D', 'S371L', 'S373P', 'S375F', 'Q493R', 'G496S', 'Q498R', 'Y505H', 'N679K', 'N764K', '69-70del', '143-145del']},
    'Omicron BA.2': {'key': ['N501Y', 'T478K', 'D614G', 'P681H', 'G446S'], 'secondary': ['S371F', 'S373P', 'S375F', 'T376A', 'D405N', 'R408S', 'K417N', 'N440K', 'Q493R', 'N764K', 'Q954H', 'N969K']},
    'Omicron BA.4/5': {'key': ['N501Y', 'T478K', 'D614G', 'P681H', 'L452R', 'F486V'], 'secondary': ['S371F', 'S373P', 'S375F', 'Q498R', 'N440K', 'K417N', 'N764K', 'Q954H', 'N969K']},
    'Lambda': {'key': ['L452Q', 'F490S', 'D614G'], 'secondary': ['G75V', 'T76I', 'D253N', 'L452Q', 'F490S', 'T859N', '246-252del']},
    'Mu': {'key': ['E484K', 'N501Y', 'D614G', 'P681H'], 'secondary': ['T95I', 'Y144S', 'Y145N', 'R346K', 'N439K', 'K458R', 'T478K', '144-145del']},
    'Epsilon': {'key': ['L452R', 'D614G'], 'secondary': ['S13I', 'W152C', 'N501Y']},
    'Zeta': {'key': ['E484K', 'D614G'], 'secondary': ['L18F', 'T20N', 'P26S', 'D138Y', 'R190S', 'H655Y', 'V1176F']},
    'Eta': {'key': ['E484K', 'D614G'], 'secondary': ['Q52R', 'A67V', '69-70del', '144del', 'N501Y', 'D614G', 'Q677H']},
    'Theta': {'key': ['E484K', 'N501Y', 'D614G', 'P681H'], 'secondary': ['H66D', '141-143del', 'Y144del', 'K417N']},
    'Iota': {'key': ['L452R', 'D614G'], 'secondary': ['T95I', 'D253G', 'E484K', 'A701V']},
    'Kappa': {'key': ['L452R', 'E484Q', 'D614G', 'P681R'], 'secondary': ['T95I', 'G142D', 'E154K', 'Q1071H']},
}

# Define flexible mutation positions
def get_mutation_range(position, range_size=10):
    return range(position - range_size, position + range_size + 1)

mutation_positions = {
    'D614G': get_mutation_range(23403),
    'N501Y': get_mutation_range(23063),
    'A570D': get_mutation_range(23271),
    'P681H': get_mutation_range(23604),
    'K417N': get_mutation_range(22813),
    'E484K': get_mutation_range(23012),
    'K417T': get_mutation_range(22813),
    'L452R': get_mutation_range(22917),
    'T478K': get_mutation_range(23012),
    'G446S': get_mutation_range(22898),
    'F486V': get_mutation_range(23018),
    'L452Q': get_mutation_range(22917),
    'F490S': get_mutation_range(23030),
    'T716I': get_mutation_range(23708),
    'S982A': get_mutation_range(24506),
    'D1118H': get_mutation_range(24914),
    'L18F': get_mutation_range(21614),
    'D80A': get_mutation_range(21802),
    'D215G': get_mutation_range(22205),
    'A701V': get_mutation_range(23664),
    'T20N': get_mutation_range(21620),
    'P26S': get_mutation_range(21638),
    'D138Y': get_mutation_range(21974),
    'R190S': get_mutation_range(22130),
    'H655Y': get_mutation_range(23525),
    'T1027I': get_mutation_range(24641),
    'V1176F': get_mutation_range(25088),
    'T19R': get_mutation_range(21617),
    'G142D': get_mutation_range(21986),
    'R158G': get_mutation_range(22034),
    'R203M': get_mutation_range(22169),
    'D950N': get_mutation_range(24410),
    'G339D': get_mutation_range(22577),
    'S371L': get_mutation_range(22673),
    'S373P': get_mutation_range(22679),
    'S375F': get_mutation_range(22685),
    'Q493R': get_mutation_range(23039),
    'G496S': get_mutation_range(23048),
    'Q498R': get_mutation_range(23054),
    'Y505H': get_mutation_range(23075),
    'N679K': get_mutation_range(23597),
    'N764K': get_mutation_range(23852),
    'S371F': get_mutation_range(22673),
    'T376A': get_mutation_range(22688),
    'D405N': get_mutation_range(22775),
    'R408S': get_mutation_range(22784),
    'N440K': get_mutation_range(22880),
    'Q954H': get_mutation_range(24422),
    'N969K': get_mutation_range(24467),
    'G75V': get_mutation_range(21785),
    'T76I': get_mutation_range(21788),
    'D253N': get_mutation_range(22319),
    'T859N': get_mutation_range(24137),
    'T95I': get_mutation_range(21845),
    'Y144S': get_mutation_range(21992),
    'Y145N': get_mutation_range(21995),
    'R346K': get_mutation_range(22598),
    'N439K': get_mutation_range(22877),
    'K458R': get_mutation_range(22934),
    'S13I': get_mutation_range(21599),
    'W152C': get_mutation_range(22016),
    'Q52R': get_mutation_range(21716),
    'A67V': get_mutation_range(21761),
    'Q677H': get_mutation_range(23591),
    'H66D': get_mutation_range(21758),
    'D253G': get_mutation_range(22319),
    'E154K': get_mutation_range(22020),
    'Q1071H': get_mutation_range(24773),
    '69-70del': (21765, 21770),
    'Y144del': (21990, 21992),
    '242-244del': (22284, 22289),
    '156-157del': (22028, 22033),
    'F157del': (22031, 22033),
    '143-145del': (21989, 21994),
    '246-252del': (22296, 22308),
    '144-145del': (21992, 21996),
    '141-143del': (21983, 21988),
}

def check_mutation(sequence, mutation, positions):
    if mutation.endswith('del'):
        start, end = positions
        return sum(nt == '-' for nt in sequence[start-1:end]) / (end - start + 1)
    else:
        target_aa = mutation[-1]
        for pos in positions:
            if pos < len(sequence):
                if sequence[pos] == target_aa:
                    return 1.0
                elif sequence[pos] in ['A', 'C', 'G', 'T']:
                    codon_start = (pos // 3) * 3
                    codon = sequence[codon_start:codon_start+3]
                    if len(codon) == 3 and genetic_code.get(codon, 'X') == target_aa:
                        return 0.8
    return 0.0


def analyze_sequences(fasta_file):
    variant_counts = Counter()
    total_sequences = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        total_sequences += 1
        sequence = str(record.seq).upper()

        mutation_profile = get_mutation_profile(sequence)
        variant = classify_variant(mutation_profile)
        variant_counts[variant] += 1

        if total_sequences % 10000 == 0:
            print(f"Processed {total_sequences} sequences...")

    return variant_counts, total_sequences

def get_top_variants(variant_counts, total_sequences, top_n=10):
    top_variants = variant_counts.most_common(top_n)
    return [(variant, count, (count / total_sequences) * 100) for variant, count in top_variants]



def print_variant_stats(variant_counts, total_sequences, file_name):
    print(f"\nVariant statistics for {file_name}:")
    print(f"Total sequences processed: {total_sequences}")
    for variant, count in variant_counts.most_common():
        percentage = (count / total_sequences) * 100
        print(f"{variant}: {count} ({percentage:.2f}%)")

def compare_variant_distributions(unfiltered_data, filtered_data):
    all_variants = set([v[0] for v in unfiltered_data + filtered_data])
    
    data = []
    for variant in all_variants:
        unfiltered = next((v for v in unfiltered_data if v[0] == variant), (variant, 0, 0))
        filtered = next((v for v in filtered_data if v[0] == variant), (variant, 0, 0))
        
        percent_change = ((filtered[2] - unfiltered[2]) / unfiltered[2] * 100) if unfiltered[2] > 0 else 0
        
        data.append({
            "Variant": variant,
            "Before %": unfiltered[2],
            "After %": filtered[2],
            "% Change": percent_change
        })
    
    return pd.DataFrame(data).sort_values("Before %", ascending=False).reset_index(drop=True)

def print_variant_table(df):
    print("\nVariant Distribution Comparison:")
    print("{:<15} {:<15} {:<15} {:<15}".format("Variant", "Before %", "After %", "% Change"))
    for _, row in df.iterrows():
        print("{:<15} {:<15.2f} {:<15.2f} {:<15.2f}".format(
            row['Variant'], row['Before %'], row['After %'], row['% Change']))

def calculate_distribution_similarity(df):
    before = df['Before %'].values
    after = df['After %'].values
    
    pearson_corr, _ = pearsonr(before, after)
    spearman_corr, _ = spearmanr(before, after)
    cosine_sim = np.dot(before, after) / (np.linalg.norm(before) * np.linalg.norm(after))
    
    print("\nDistribution Similarity Metrics:")
    print(f"Pearson correlation: {pearson_corr:.4f}")
    print(f"Spearman correlation: {spearman_corr:.4f}")
    print(f"Cosine similarity: {cosine_sim:.4f}")


def plot_mutation_distribution(unfiltered_counts, filtered_counts, output_file):
    mutations = list(unfiltered_counts.keys())
    unfiltered_percentages = [count / sum(unfiltered_counts.values()) * 100 for count in unfiltered_counts.values()]
    filtered_percentages = [filtered_counts[m] / sum(filtered_counts.values()) * 100 for m in mutations]

    plt.figure(figsize=(15, 10))
    plt.bar(range(len(mutations)), unfiltered_percentages, width=0.4, label='Unfiltered', align='center')
    plt.bar(range(len(mutations)), filtered_percentages, width=0.4, label='Filtered', align='edge')
    plt.xlabel('Mutations')
    plt.ylabel('Percentage')
    plt.title('Mutation Distribution Before and After Filtering')
    plt.xticks(range(len(mutations)), mutations, rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def generate_heatmap(df, output_file):
    plt.figure(figsize=(12, 8))
    sns.set(style="whitegrid")
    
    data = df[['Before %', 'After %']].transpose()
    sns.heatmap(data, annot=True, cmap='YlOrRd', fmt='.2f')
    
    plt.title('Variant Distribution Heatmap', fontsize=14)
    plt.xlabel('Variants', fontsize=12)
    plt.ylabel('Before/After Filtering', fontsize=12)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def check_mutation_confidence(sequence, mutation, positions):
    if mutation.endswith('del'):
        start, end = positions
        deletion_count = sum(nt == '-' for nt in sequence[start-1:end])
        return deletion_count / (end - start + 1)
    else:
        target_aa = mutation[-1]
        max_confidence = 0
        for pos in positions:
            if pos < len(sequence):
                if sequence[pos] == target_aa:
                    max_confidence = max(max_confidence, 1.0)
                elif sequence[pos] in ['A', 'C', 'G', 'T']:
                    codon_start = (pos // 3) * 3
                    codon = sequence[codon_start:codon_start+3]
                    if len(codon) == 3:
                        amino_acid = genetic_code.get(codon, 'X')
                        if amino_acid == target_aa:
                            max_confidence = max(max_confidence, 0.8)
        return max_confidence

def plot_variant_distribution(df, output_file):
    plt.figure(figsize=(12, 6))
    sns.set(style="whitegrid")

    x = range(len(df))
    width = 0.35

    plt.bar(x, df['Before %'], width, label='Before Filtering', color='blue', alpha=0.7)
    plt.bar([i + width for i in x], df['After %'], width, label='After Filtering', color='red', alpha=0.7)

    plt.xlabel('Variants', fontsize=12)
    plt.ylabel('Percentage', fontsize=12)
    plt.title('Variant Distribution Before and After Filtering', fontsize=14)
    plt.xticks([i + width/2 for i in x], df['Variant'], rotation=45, ha='right', fontsize=10)
    plt.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def generate_consistency_heatmap(df, output_file):
    plt.figure(figsize=(10, 6))
    sns.set(style="whitegrid")

    data = df[['Before %', 'After %']].transpose()
    sns.heatmap(data, annot=True, cmap='YlOrRd', fmt='.2f', cbar_kws={'label': 'Percentage'})

    plt.title('Variant Distribution Consistency Heatmap', fontsize=14)
    plt.xlabel('Variants', fontsize=12)
    plt.ylabel('Before/After Filtering', fontsize=12)
    plt.xticks(rotation=45, ha='right', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def classify_variant(mutation_profile):
    variant_scores = {}
    for variant, mutations in variant_mutations.items():
        key_score = sum(mutation_profile.get(m, 0) for m in mutations['key']) / len(mutations['key']) if mutations['key'] else 0
        secondary_score = sum(mutation_profile.get(m, 0) for m in mutations['secondary']) / len(mutations['secondary']) if mutations['secondary'] else 0
        variant_scores[variant] = (key_score * 0.7) + (secondary_score * 0.3)
    
    best_variant = max(variant_scores, key=variant_scores.get)
    best_score = variant_scores[best_variant]
    
    if best_score >= 0.5:
        return best_variant
    elif mutation_profile.get('D614G', 0) > 0.8:
        return "D614G"
    elif not any(mutation_profile.values()):
        return "Ancestral"
    else:
        return "Other"

def get_mutation_profile(sequence):
    return {mutation: check_mutation(sequence, mutation, positions) 
            for mutation, positions in mutation_positions.items()}

if __name__ == "__main__":
    unfiltered_file = "SARS-CoV-2.fa"
    filtered_file = "sorted_output.fa"
    output_dir = "analysis_results"

    os.makedirs(output_dir, exist_ok=True)

    print(f"Processing unfiltered file: {unfiltered_file}")
    unfiltered_variant_counts, unfiltered_total = analyze_sequences(unfiltered_file)
    unfiltered_top_variants = get_top_variants(unfiltered_variant_counts, unfiltered_total)

    print(f"\nProcessing filtered file: {filtered_file}")
    filtered_variant_counts, filtered_total = analyze_sequences(filtered_file)
    filtered_top_variants = get_top_variants(filtered_variant_counts, filtered_total)

    df_variants = compare_variant_distributions(unfiltered_top_variants, filtered_top_variants)
    
    print_variant_table(df_variants)
    calculate_distribution_similarity(df_variants)

    # Save results to CSV
    df_variants.to_csv(os.path.join(output_dir, "variant_distribution_comparison.csv"), index=False)

    # Generate visualizations
    plot_variant_distribution(df_variants, os.path.join(output_dir, "variant_distribution_plot.pdf"))
    generate_consistency_heatmap(df_variants, os.path.join(output_dir, "variant_consistency_heatmap.pdf"))

    print(f"\nAnalysis complete. Results saved in {output_dir}")