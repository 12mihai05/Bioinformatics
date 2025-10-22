import matplotlib.pyplot as plt
from collections import Counter
from pathlib import Path
import re
import json

GENETIC_CODE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def clean_dna(seq: str) -> str:
    """Keep only A/C/G/T, uppercase."""
    return "".join(re.findall(r"[ACGTacgt]", seq)).upper()

def read_fasta(filepath: str) -> str:
    """Read FASTA file and return concatenated sequence."""
    text = Path(filepath).read_text()
    if text.startswith(">"):
        lines = [ln.strip() for ln in text.splitlines() if not ln.startswith(">")]
        return "".join(lines)
    return text.strip()

def count_codons(dna_sequence: str) -> Counter:
    """Count all codons in a DNA sequence (reading frame 0)."""
    dna = clean_dna(dna_sequence)
    codons = []

    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        if len(codon) == 3:
            codons.append(codon)

    return Counter(codons)

def get_amino_acid_frequencies(codon_counts: Counter) -> Counter:
    """Convert codon frequencies to amino acid frequencies."""
    aa_counts = Counter()
    for codon, count in codon_counts.items():
        aa = GENETIC_CODE.get(codon, "?")
        if aa != "?":
            aa_counts[aa] += count
    return aa_counts

def plot_top_codons(codon_counts: Counter, title: str, filename: str, top_n: int = 10):
    """Create a bar chart of top N most frequent codons."""
    top_codons = codon_counts.most_common(top_n)
    codons, counts = zip(*top_codons) if top_codons else ([], [])

    plt.figure(figsize=(12, 8))
    bars = plt.bar(codons, counts)
    plt.title(f'{title} - Top {top_n} Most Frequent Codons', fontsize=14, fontweight='bold')
    plt.xlabel('Codons', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.xticks(rotation=45)

    # Add amino acid labels on bars
    for i, (codon, count) in enumerate(top_codons):
        aa = GENETIC_CODE.get(codon, "?")
        plt.text(i, count + max(counts) * 0.01, f'{aa}',
                ha='center', va='bottom', fontweight='bold')

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def compare_genomes(covid_codons: Counter, flu_codons: Counter):
    """Compare codon frequencies between the two genomes."""
    print("\n" + "="*80)
    print("GENOME COMPARISON - MOST FREQUENT CODONS")
    print("="*80)

    covid_top = dict(covid_codons.most_common(10))
    flu_top = dict(flu_codons.most_common(10))

    # Find common codons in top 10 of both
    common_codons = set(covid_top.keys()) & set(flu_top.keys())

    print(f"\nCommon codons in top 10 of both genomes: {len(common_codons)}")
    print("-" * 50)

    if common_codons:
        print(f"{'Codon':<8}{'AA':<4}{'COVID-19':<12}{'Influenza':<12}{'Ratio (C/I)':<12}")
        print("-" * 50)
        for codon in sorted(common_codons, key=lambda x: covid_codons[x], reverse=True):
            aa = GENETIC_CODE.get(codon, "?")
            covid_count = covid_codons[codon]
            flu_count = flu_codons[codon]
            ratio = covid_count / flu_count if flu_count > 0 else float('inf')
            print(f"{codon:<8}{aa:<4}{covid_count:<12}{flu_count:<12}{ratio:<12.2f}")

    # Overall most frequent across both
    print(f"\nOverall most frequent codons (combined):")
    print("-" * 50)
    combined_codons = covid_codons + flu_codons
    print(f"{'Codon':<8}{'AA':<4}{'Combined Count':<15}")
    print("-" * 50)
    for codon, count in combined_codons.most_common(10):
        aa = GENETIC_CODE.get(codon, "?")
        print(f"{codon:<8}{aa:<4}{count:<15}")

def print_top_amino_acids(aa_counts: Counter, genome_name: str):
    """Print top 3 amino acids for a genome."""
    print(f"\nTop 3 amino acids in {genome_name}:")
    print("-" * 40)

    # Full amino acid names
    aa_names = {
        "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic acid",
        "C": "Cysteine", "Q": "Glutamine", "E": "Glutamic acid", "G": "Glycine",
        "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
        "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
        "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
        "*": "Stop codon"
    }

    for i, (aa, count) in enumerate(aa_counts.most_common(3), 1):
        aa_name = aa_names.get(aa, "Unknown")
        percentage = (count / sum(aa_counts.values())) * 100
        print(f"{i}. {aa_name} ({aa}): {count:,} ({percentage:.1f}%)")

def save_results_to_file(covid_aa, flu_aa, covid_codons, flu_codons):
    """Save analysis results to JSON file for use by nutritional analysis."""

    # Get top amino acids for each virus
    covid_top_aa = [(aa, count, (count/sum(covid_aa.values()))*100)
                    for aa, count in covid_aa.most_common(5)]
    flu_top_aa = [(aa, count, (count/sum(flu_aa.values()))*100)
                  for aa, count in flu_aa.most_common(5)]

    # Full amino acid names
    aa_names = {
        "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic acid",
        "C": "Cysteine", "Q": "Glutamine", "E": "Glutamic acid", "G": "Glycine",
        "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
        "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
        "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
        "*": "Stop codon"
    }

    results = {
        "analysis_metadata": {
            "covid_genome_length": len(clean_dna(read_fasta("covid19_genome.fasta"))),
            "influenza_genome_length": len(clean_dna(read_fasta("influenza_genome.fasta"))),
            "covid_total_codons": sum(covid_codons.values()),
            "influenza_total_codons": sum(flu_codons.values())
        },
        "covid_top_amino_acids": [
            {
                "amino_acid_code": aa,
                "amino_acid_name": aa_names.get(aa, "Unknown"),
                "count": count,
                "percentage": round(percentage, 2)
            }
            for aa, count, percentage in covid_top_aa
        ],
        "influenza_top_amino_acids": [
            {
                "amino_acid_code": aa,
                "amino_acid_name": aa_names.get(aa, "Unknown"),
                "count": count,
                "percentage": round(percentage, 2)
            }
            for aa, count, percentage in flu_top_aa
        ],
    }

    covid_top_names = [aa_names.get(aa, "Unknown") for aa, _, _ in covid_top_aa[:3] if aa != "*"]
    flu_top_names = [aa_names.get(aa, "Unknown") for aa, _, _ in flu_top_aa[:3] if aa != "*"]
    results["target_amino_acids"] = list(set(covid_top_names + flu_top_names))

    with open("codon_analysis_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n💾 Results saved to 'codon_analysis_results.json'")
    return results

def main():
    print("Analyzing COVID-19 and Influenza Genome Codon Frequencies")
    print("=" * 60)

    print("\nReading genome files...")
    covid_dna = read_fasta("covid19_genome.fasta")
    flu_dna = read_fasta("influenza_genome.fasta")

    print(f"COVID-19 genome length: {len(clean_dna(covid_dna)):,} bp")
    print(f"Influenza genome length: {len(clean_dna(flu_dna)):,} bp")

    print("\nCounting codons...")
    covid_codons = count_codons(covid_dna)
    flu_codons = count_codons(flu_dna)

    print(f"COVID-19 total codons: {sum(covid_codons.values()):,}")
    print(f"Influenza total codons: {sum(flu_codons.values()):,}")

    print("\nGenerating charts...")
    plot_top_codons(covid_codons, "COVID-19 Genome", "covid19_top10_codons.png")
    plot_top_codons(flu_codons, "Influenza Genome", "influenza_top10_codons.png")

    compare_genomes(covid_codons, flu_codons)

    covid_aa = get_amino_acid_frequencies(covid_codons)
    flu_aa = get_amino_acid_frequencies(flu_codons)

    print_top_amino_acids(covid_aa, "COVID-19")
    print_top_amino_acids(flu_aa, "Influenza")

    save_results_to_file(covid_aa, flu_aa, covid_codons, flu_codons)

    print(f"\n{'='*60}")
    print("Analysis complete! Charts saved as PNG files.")
    print("='*60")

if __name__ == "__main__":
    main()