# Take an arbitrary DNA sequence from NCBI, between 1000 and 3000 nucleotides (letters).
# Implement a software application that detects repetitions between 3b to 6b in this DNA sequence.
# NOTE: repetitive sequences refer to patterns that repeat N times. Minimum number of repetitions is 2.

def read_fasta(filename):
    sequence = ""
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith('>'): 
                sequence += line.strip().upper()
    return sequence

def find_repetitions(sequence, min_length=3, max_length=6, min_repeats=2):
    results = {}
    
    for pattern_length in range(min_length, max_length + 1):
        pattern_counts = {}
        
        for i in range(len(sequence) - pattern_length + 1):
            pattern = sequence[i:i + pattern_length]
            
            if all(base in 'ACGT' for base in pattern):
                pattern_counts[pattern] = pattern_counts.get(pattern, 0) + 1
        
        repetitive_patterns = {pattern: count for pattern, count in pattern_counts.items() 
                              if count >= min_repeats}
        
        if repetitive_patterns:
            results[pattern_length] = dict(sorted(repetitive_patterns.items(), 
                                                 key=lambda x: (-x[1], x[0])))
    
    return results

def display_results(results, sequence_length, top_n=10):
    print(f"DNA Sequence Analysis")
    print(f"Sequence length: {sequence_length} nucleotides")
    print(f"\nRepetitive patterns (3-6 bp) appearing at least 2 times:\n")
    print("=" * 70)
    
    for length in sorted(results.keys()):
        patterns = results[length]
        print(f"\n{length}-base pair patterns: {len(patterns)} unique repetitive patterns found")
        print("-" * 70)
        
        count = 0
        for pattern, occurrences in patterns.items():
            if count >= top_n:
                remaining = len(patterns) - top_n
                if remaining > 0:
                    print(f"... and {remaining} more patterns")
                break
            print(f"  {pattern}: {occurrences} occurrences")
            count += 1
    
    print("\n" + "=" * 70)
    print("Summary:")
    total_patterns = sum(len(patterns) for patterns in results.values())
    print(f"Total unique repetitive patterns found: {total_patterns}")
    
    max_count = 0
    max_pattern = ""
    for patterns in results.values():
        for pattern, count in patterns.items():
            if count > max_count:
                max_count = count
                max_pattern = pattern
    
    if max_pattern:
        print(f"Most frequent pattern: {max_pattern} ({max_count} occurrences)")

def main():
    filename = "sequence.fasta"
    print(f"Reading sequence from {filename}...")
    
    try:
        sequence = read_fasta(filename)
        
        if len(sequence) < 1000:
            print(f"Warning: Sequence is only {len(sequence)} bp (recommended: 1000-3000 bp)")
        elif len(sequence) > 3000:
            print(f"Sequence is {len(sequence)} bp. Analyzing first 3000 bp...")
            sequence = sequence[:3000]
        
        print(f"Analyzing {len(sequence)} nucleotides...\n")
        
        results = find_repetitions(sequence, min_length=3, max_length=6, min_repeats=2)
        
        display_results(results, len(sequence), top_n=10)
        
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found!")
        print("Please ensure the FASTA file is in the same directory as this script.")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()