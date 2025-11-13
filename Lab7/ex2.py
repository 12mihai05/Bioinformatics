# Download 10 influenza genomes. For each genome plot on a chart the most frequent repetitions.

from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

# Set your email for NCBI
Entrez.email = "pasaroiumihai@yahoo.com"

def download_influenza_genomes(count=10):
    """
    Download influenza genome sequences from NCBI.
    
    Args:
        count: Number of genomes to download
    
    Returns:
        List of tuples (accession_id, sequence)
    """
    print(f"Searching for {count} influenza genomes...")
    
    # Search for influenza complete genome sequences
    search_handle = Entrez.esearch(
        db="nucleotide",
        term="(Influenza A virus[Organism] OR Influenza B virus[Organism]) AND (complete genome OR complete sequence) AND 1000:15000[SLEN]",
        retmax=count,
        sort="relevance"
    )
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    id_list = search_results["IdList"]
    print(f"Found {len(id_list)} genome IDs")
    
    genomes = []
    
    for i, genome_id in enumerate(id_list, 1):
        print(f"Downloading genome {i}/{len(id_list)} (ID: {genome_id})...")
        
        try:
            # Fetch the sequence
            fetch_handle = Entrez.efetch(
                db="nucleotide",
                id=genome_id,
                rettype="fasta",
                retmode="text"
            )
            
            record = SeqIO.read(fetch_handle, "fasta")
            fetch_handle.close()
            
            sequence = str(record.seq).upper()
            accession = record.id
            
            # Limit to 3000 bp for consistency
            if len(sequence) > 3000:
                sequence = sequence[:3000]
            
            genomes.append((accession, sequence))
            print(f"  {accession}: {len(sequence)} bp")
            
        except Exception as e:
            print(f"  Error downloading genome {genome_id}: {e}")
            continue
    
    return genomes

def find_most_frequent_repetitions(sequence, pattern_length=4, top_n=10):
    """
    Find the most frequent repetitive patterns in a sequence.
    
    Args:
        sequence: DNA sequence string
        pattern_length: Length of patterns to search for
        top_n: Number of top patterns to return
    
    Returns:
        List of tuples (pattern, count) sorted by frequency
    """
    pattern_counts = {}
    
    # Slide through the sequence
    for i in range(len(sequence) - pattern_length + 1):
        pattern = sequence[i:i + pattern_length]
        
        # Only consider valid DNA bases
        if all(base in 'ACGT' for base in pattern):
            pattern_counts[pattern] = pattern_counts.get(pattern, 0) + 1
    
    # Filter patterns that appear at least 2 times and get top N
    repetitive_patterns = [(pattern, count) for pattern, count in pattern_counts.items() 
                          if count >= 2]
    repetitive_patterns.sort(key=lambda x: x[1], reverse=True)
    
    return repetitive_patterns[:top_n]

def plot_repetitions(genomes_data, pattern_length=4):
    """
    Create visualizations for the most frequent repetitions across genomes.
    
    Args:
        genomes_data: List of tuples (accession, sequence)
        pattern_length: Length of patterns to analyze
    """
    num_genomes = len(genomes_data)
    
    # Create a figure with subplots
    fig, axes = plt.subplots(2, 5, figsize=(20, 10))
    fig.suptitle(f'Most Frequent {pattern_length}-bp Repetitions in Influenza Genomes', 
                 fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    for idx, (accession, sequence) in enumerate(genomes_data):
        ax = axes[idx]
        
        # Find top repetitions
        top_patterns = find_most_frequent_repetitions(sequence, pattern_length, top_n=10)
        
        if top_patterns:
            patterns = [p[0] for p in top_patterns]
            counts = [p[1] for p in top_patterns]
            
            # Create bar chart
            bars = ax.barh(range(len(patterns)), counts, color='steelblue', edgecolor='black')
            ax.set_yticks(range(len(patterns)))
            ax.set_yticklabels(patterns, fontsize=8, fontfamily='monospace')
            ax.set_xlabel('Occurrences', fontsize=9)
            ax.set_title(f'{accession}\n({len(sequence)} bp)', fontsize=9, fontweight='bold')
            ax.invert_yaxis()
            ax.grid(axis='x', alpha=0.3)
            
            # Add count labels on bars
            for i, (bar, count) in enumerate(zip(bars, counts)):
                ax.text(count + 0.5, i, str(count), va='center', fontsize=7)
        else:
            ax.text(0.5, 0.5, 'No repetitions found', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{accession}', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('influenza_repetitions.png', dpi=300, bbox_inches='tight')
    print(f"\nPlot saved as 'influenza_repetitions.png'")
    plt.show()

def plot_comparison_summary(genomes_data, pattern_length=4):
    """
    Create a summary comparison plot showing the top pattern for each genome.
    
    Args:
        genomes_data: List of tuples (accession, sequence)
        pattern_length: Length of patterns to analyze
    """
    genome_labels = []
    top_pattern_names = []
    top_pattern_counts = []
    
    for accession, sequence in genomes_data:
        top_patterns = find_most_frequent_repetitions(sequence, pattern_length, top_n=1)
        
        if top_patterns:
            pattern, count = top_patterns[0]
            genome_labels.append(accession.split('.')[0][:15])  # Truncate for readability
            top_pattern_names.append(pattern)
            top_pattern_counts.append(count)
    
    # Create comparison plot
    fig, ax = plt.subplots(figsize=(14, 6))
    
    bars = ax.bar(range(len(genome_labels)), top_pattern_counts, 
                   color='coral', edgecolor='black', alpha=0.7)
    
    ax.set_xlabel('Genome', fontsize=12, fontweight='bold')
    ax.set_ylabel('Occurrences of Most Frequent Pattern', fontsize=12, fontweight='bold')
    ax.set_title(f'Most Frequent {pattern_length}-bp Pattern in Each Influenza Genome', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(genome_labels)))
    ax.set_xticklabels(genome_labels, rotation=45, ha='right', fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    
    # Add pattern labels on bars
    for i, (bar, pattern, count) in enumerate(zip(bars, top_pattern_names, top_pattern_counts)):
        ax.text(i, count + 1, f'{pattern}\n({count})', 
               ha='center', va='bottom', fontsize=8, fontweight='bold',
               fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig('influenza_comparison.png', dpi=300, bbox_inches='tight')
    print(f"Comparison plot saved as 'influenza_comparison.png'")
    plt.show()

def main():
    print("=" * 70)
    print("INFLUENZA GENOME REPETITION ANALYSIS")
    print("=" * 70)
    
    # Download genomes
    genomes = download_influenza_genomes(count=10)
    
    if not genomes:
        print("Error: No genomes were downloaded successfully.")
        return
    
    print(f"\nSuccessfully downloaded {len(genomes)} genomes")
    print("\nAnalyzing repetitions...")
    
    # Analyze and plot for different pattern lengths
    for pattern_length in [3, 4, 5, 6]:
        print(f"\nAnalyzing {pattern_length}-bp patterns...")
        
        # Show statistics
        for accession, sequence in genomes:
            top_patterns = find_most_frequent_repetitions(sequence, pattern_length, top_n=3)
            if top_patterns:
                print(f"  {accession}: Top pattern: {top_patterns[0][0]} ({top_patterns[0][1]} occurrences)")
    
    # Create visualizations
    print("\n" + "=" * 70)
    print("Generating plots...")
    plot_repetitions(genomes, pattern_length=4)
    plot_comparison_summary(genomes, pattern_length=4)
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()