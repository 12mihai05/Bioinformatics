# Implement an app that converts the coding region of a gene into an amino acid sequence.
# Use the genetic code table from moodle.

import re

DNA_SEQ = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

GENETIC_CODE = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*",
    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",

    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",

    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",

    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
}

AA3 = {
    "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys",
    "Q":"Gln","E":"Glu","G":"Gly","H":"His","I":"Ile",
    "L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro",
    "S":"Ser","T":"Thr","W":"Trp","Y":"Tyr","V":"Val",
    "*":"Stop"
}

def clean_dna(seq: str) -> str:
    return "".join(re.findall(r"[ACGTacgt]", seq)).upper()

def dna_to_rna(dna: str) -> str:
    return dna.replace("T", "U")

def translate(dna: str) -> tuple[str, str]:
    dna = clean_dna(dna)
    rna = dna_to_rna(dna)
    aa = []
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        res = GENETIC_CODE.get(codon, "?")
        if res == "*" or res == "?":
            break
        aa.append(res)
    one = "".join(aa)
    three = "-".join(AA3[a] for a in aa)
    return one, three

if __name__ == "__main__":
    one, three = translate(DNA_SEQ)
    print("DNA :", DNA_SEQ)
    print("AA  (1-letter):", one)
    print("AA (3-letter):", three)
