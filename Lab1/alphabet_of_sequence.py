# Make an application that is able to find the alphabet a sequence of text. 
# This seq. may be an ARN seq. or ADN seq. or protein seq.

def find_alphabet(sequence: str) -> set:

    return set(sequence)

if __name__ == "__main__":
    seq = "ABBBABBBABABABABABC"
    alphabet = find_alphabet(seq)
    print("Alphabet of sequence:", alphabet)

