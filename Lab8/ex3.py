# Download from NCBI 3 bacterial genomes of your choosing try to find in these genomes
# possible transposons. For this one must detect possible inverted repeats without prior
# knowleadge about their existance in the sequence. The inverted repeat should have a 
# minimum length of 4 letters and a maximum of 6 letters.

from Bio import Entrez, SeqIO
import random

Entrez.email = "pasaroiumihai@yahoo.com"  # REQUIRED by NCBI

GENOMES = {
    "E_coli_MG1655": "NC_000913.3",
    "B_subtilis_168": "NC_000964.3",
    "P_aeruginosa_PAO1": "NC_002516.2",
}

def fetch_genome(accession: str) -> str:
    """Download a genome from NCBI (nuccore) in FASTA and return the sequence as a string."""
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq).upper()

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]

def find_inverted_repeats(seq: str, min_len: int = 4, max_len: int = 6, max_spacer: int = 30):
    n = len(seq)
    results = []

    for i in range(n - min_len + 1):
        for L in range(min_len, max_len + 1):
            if i + L > n:
                continue

            left = seq[i:i+L]
            rc = revcomp(left)

            min_j = i + L
            max_j = min(n - L, i + L + max_spacer)

            for j in range(min_j, max_j + 1):
                right = seq[j:j+L]
                if right == rc:
                    spacer = j - (i + L)
                    results.append({
                        "left_start": i,
                        "left_end": i + L - 1,
                        "right_start": j,
                        "right_end": j + L - 1,
                        "length": L,
                        "spacer": spacer,
                        "motif": left,
                    })
    return results

def main():
    genomes_seqs = {}
    for name, acc in GENOMES.items():
        print(f"Downloading genome {name} ({acc}) from NCBI...")
        seq = fetch_genome(acc)
        genomes_seqs[name] = seq
        print(f"  -> length: {len(seq)} bases")

    print("\n=== INVERTED REPEAT SEARCH (length 4–6, spacer ≤ 30) ===\n")

    for name, seq in genomes_seqs.items():
        print(f"--- Genome: {name} ---")
        irs = find_inverted_repeats(seq, min_len=4, max_len=6, max_spacer=30)
        print(f"Total inverted repeats found: {len(irs)}")

        max_show = 10
        for idx, ir in enumerate(irs[:max_show], start=1):
            print(
                f"IR {idx}: motif={ir['motif']}, len={ir['length']}, spacer={ir['spacer']} | "
                f"left=({ir['left_start']}-{ir['left_end']}), right=({ir['right_start']}-{ir['right_end']})"
            )
        if len(irs) > max_show:
            print(f"... ({len(irs) - max_show} more IRs not shown)\n")
        else:
            print()

    print("Note: Positions are 0-based indices. Inverted repeats are candidates for transposon boundaries.")

if __name__ == "__main__":
    main()
