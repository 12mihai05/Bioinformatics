# The melting point temp (tm) is the temperature at which one half of a particular DNA will disociate and become 
# a single strand of DNA primerland and seq are of critical importance in designing the parameters of a successful
# amplification. The tm of a nucleic acid duplex increases both with its length and with increasing GC content.
# tm = 4(G + C) + 2(A + T)
# tm = 81.5 + 16.6(log 10([Na+])) + 0.41*(%GC) - 600/length
# Implement the application that computes the tm of a DNA seq by using one of these formulas or both of them
# Input = a string of DNA
# Output = temp in Celsius

import math

def calculate_melting_temperature1(dna_sequence):
    dna_sequence = dna_sequence.upper()
    
    g_count = dna_sequence.count('G')
    c_count = dna_sequence.count('C')
    a_count = dna_sequence.count('A')
    t_count = dna_sequence.count('T')
    
    tm = 4 * (g_count + c_count) + 2 * (a_count + t_count)
    
    return tm

def calculate_melting_temperature2(dna_sequence, salt_concentration=0.05):
    dna_sequence = dna_sequence.upper()
    
    length = len(dna_sequence)
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    gc_percentage = (gc_count / length) * 100
    
    tm = 81.5 + 16.6 * math.log10(salt_concentration) + 0.41 * gc_percentage - 600 / length
    
    return tm

def main():
    print("DNA Melting Temperature Calculator")
    print("=" * 40)
    
    dna_sequence = input("Enter DNA sequence (A, T, G, C): ").strip()
    
    valid_nucleotides = set('ATGC')
    if not dna_sequence:
        print("Error: Please enter a valid DNA sequence.")
        return
    
    if not all(nucleotide.upper() in valid_nucleotides for nucleotide in dna_sequence):
        print("Error: DNA sequence should only contain A, T, G, C nucleotides.")
        return
    
    salt_conc = 0.01
    
    tm_simple = calculate_melting_temperature1(dna_sequence)
    tm_complex = calculate_melting_temperature2(dna_sequence, salt_conc)

    dna_upper = dna_sequence.upper()
    length = len(dna_upper)
    gc_count = dna_upper.count('G') + dna_upper.count('C')
    gc_percentage = (gc_count / length) * 100
    
    print(f"\n" + "=" * 50)
    print("RESULTS")
    print("=" * 50)
    print(f"DNA Sequence: {dna_upper}")
    print(f"Length: {length} nucleotides")
    print(f"GC Content: {gc_percentage:.1f}%")
    print(f"Salt Concentration: {salt_conc} M")
    print(f"\nMelting Temperature Calculations:")
    print(f"  Simple Formula [4(G+C) + 2(A+T)]: {tm_simple:.1f}°C")
    print(f"  Complex Formula [81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length]: {tm_complex:.1f}°C")
    print(f"  Difference: {abs(tm_simple - tm_complex):.1f}°C")

if __name__ == "__main__":
    main() 