# A DNA seq. is given: S="ACGGGCATATGCGC".
# Make an app which is able to show the percentage of the components from the alphabet of the seq. S. 
# In other words, the input of the seq. S and the output is the alphabet of the seq. 
# and the percentage of each letter in the alphabet found in seq. S.

def find_alphabet(sequence: str) -> set:

    return set(sequence)

if __name__ == "__main__":
    seq = "ABB"
    alphabet = find_alphabet(seq)
    print("Alphabet of sequence:", alphabet)
    for char in alphabet:
        percentage = seq.count(char)/len(seq) * 100
        print(percentage)
