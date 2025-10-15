# Design an app that uses the sliding window method in order to read the melting temperature (tm) over the seq S.
# Use a sliding window of 8 positions and use the fasta file (my_sequence.fasta) as input.

import matplotlib.pyplot as plt
import numpy as np

def read_fasta(filename):

    sequences = {}
    with open(filename, 'r') as file:
        current_id = None
        current_seq = ""
        
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = current_seq
                current_id = line[1:]
                current_seq = ""
            else:
                current_seq += line.upper()
        
        if current_id is not None:
            sequences[current_id] = current_seq
    
    return sequences

def calculate_tm1(sequence):
    if len(sequence) == 0:
        return 0
    
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    
    tm1 = 4 * (g_count + c_count) + 2 * (a_count + t_count)
    return round(tm1, 2)

def calculate_tm2(sequence, na_concentration=0.01):
    if len(sequence) == 0:
        return 0
    
    import math
    
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    total_length = len(sequence)
    
    gc_percent = ((g_count + c_count) / total_length) * 100
    
    tm2 = 81.5 + 16.6 * math.log10(na_concentration) + 0.41 * gc_percent - 600 / total_length
    
    return round(tm2, 2)

def calculate_both_tm(sequence):
    tm1 = calculate_tm1(sequence)
    tm2 = calculate_tm2(sequence)
    return tm1, tm2

def plot_tm_values(seq_id, results, window_size):
    if not results:
        return
    
    positions = [pos for pos, _, _, _ in results]
    tm1_values = [tm1 for _, _, tm1, _ in results]
    tm2_values = [tm2 for _, _, _, tm2 in results]
    
    plt.figure(figsize=(12, 6))
    
    plt.plot(positions, tm1_values, 'b-o', label='Tm1 (Simple formula)', linewidth=2, markersize=4)
    plt.plot(positions, tm2_values, 'r-o', label='Tm2 (GC content formula)', linewidth=2, markersize=4)
    
    plt.xlabel('Window Position')
    plt.ylabel('Melting Temperature (°C)')
    plt.title(f'Melting Temperature Profile - {seq_id}\n(Window size: {window_size})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.ylim(min(min(tm1_values), min(tm2_values)) - 2, 
             max(max(tm1_values), max(tm2_values)) + 2)
    
    plt.xticks(positions)
    
    plt.tight_layout()
    plt.show()

def sliding_window_tm(sequence, window_size=8):
    results = []
    
    if len(sequence) < window_size:
        print(f"Warning: Sequence length ({len(sequence)}) is shorter than window size ({window_size})")
        return results
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        tm1, tm2 = calculate_both_tm(window)
        results.append((i, window, tm1, tm2))
    
    return results

def analyze_sequences(filename, window_size=8):
    print(f"Reading sequences from: {filename}")
    print(f"Using sliding window size: {window_size}")
    print("=" * 60)
    
    try:
        sequences = read_fasta(filename)
        
        for seq_id, sequence in sequences.items():
            print(f"\nAnalyzing sequence: {seq_id}")
            print(f"Full sequence: {sequence}")
            print(f"Sequence length: {len(sequence)} bases")
            print("-" * 40)
            
            results = sliding_window_tm(sequence, window_size)
            
            if not results:
                print("No windows could be analyzed (sequence too short)")
                continue
            
            print(f"{'Position':<10} {'Window':<12} {'Tm1 (°C)':<10} {'Tm2 (°C)':<10}")
            print("-" * 50)
            
            for position, window, tm1, tm2 in results:
                print(f"{position:<10} {window:<12} {tm1:<10} {tm2:<10}")
            
            tm1_values = [tm1 for _, _, tm1, _ in results]
            tm2_values = [tm2 for _, _, _, tm2 in results]
            
            avg_tm1 = sum(tm1_values) / len(tm1_values)
            min_tm1 = min(tm1_values)
            max_tm1 = max(tm1_values)
            
            avg_tm2 = sum(tm2_values) / len(tm2_values)
            min_tm2 = min(tm2_values)
            max_tm2 = max(tm2_values)
            
            print(f"\nStatistics for {seq_id}:")
            print(f"  Tm1 - Average: {avg_tm1:.2f}°C, Min: {min_tm1:.2f}°C, Max: {max_tm1:.2f}°C")
            print(f"  Tm2 - Average: {avg_tm2:.2f}°C, Min: {min_tm2:.2f}°C, Max: {max_tm2:.2f}°C")
            
            plot_tm_values(seq_id, results, window_size)
            
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found!")
    except Exception as e:
        print(f"Error reading file: {e}")

if __name__ == "__main__":
    fasta_file = "C:/Users/pasar/Desktop/Facultate/Bioinformatics/Lab3/fasta.fasta"
    window_size = 8
    
    analyze_sequences(fasta_file, window_size)
