import re
import math
import matplotlib.pyplot as plt

ENZYMES = {
    "EcoRI":   {"site": "GAATTC", "offset": 1},
    "BamHI":   {"site": "GGATCC", "offset": 1},
    "HindIII": {"site": "AAGCTT", "offset": 1},
    "TaqI":    {"site": "TCGA",   "offset": 1},
    "HaeIII":  {"site": "GGCC",   "offset": 2},
}

def get_cut_sites(dna, enzyme_cfg):
    seq = dna.upper()
    pattern = enzyme_cfg["site"].upper()
    offset = enzyme_cfg["offset"]

    hits = [m.start() + offset for m in re.finditer(f"(?={pattern})", seq)]
    return hits

def compute_fragments(dna, enzymes):
    info = {}
    n = len(dna)

    for name, cfg in enzymes.items():
        cuts = get_cut_sites(dna, cfg)
        boundaries = [0] + sorted(cuts) + [n]
        fragments = [
            boundaries[i+1] - boundaries[i] for i in range(len(boundaries) - 1)
        ]

        info[name] = {
            "cuts": sorted(cuts),
            "sizes": sorted(fragments, reverse=True),
        }

    return info

def draw_gel(digest_info):
    plt.figure(figsize=(7, 9))
    ax = plt.gca()
    ax.set_facecolor("black")

    lane_gap = 3.0
    band_half_width = 0.5
    gel_height = 500 

    all_sizes = [s for v in digest_info.values() for s in v["sizes"]]
    max_bp = max(all_sizes)
    min_bp = min(all_sizes)

    def mobility(bp):
        if max_bp == min_bp:
            return gel_height / 2
        frac = (math.log(bp) - math.log(min_bp)) / (math.log(max_bp) - math.log(min_bp))
        return 50 + (1 - frac) * gel_height

    for lane_idx, (name, data) in enumerate(digest_info.items()):
        x_center = 2 + lane_idx * lane_gap

        plt.hlines(gel_height + 60, x_center - 0.9, x_center + 0.9,
                   linewidth=6, color="white")

        for size in data["sizes"]:
            y = mobility(size)
            plt.hlines(y, x_center - band_half_width, x_center + band_half_width,
                       linewidth=4, color="white")

        plt.text(x_center, gel_height + 30, name, color="white",
                 ha="center", va="center", fontsize=9)

    plt.xlim(0, lane_gap * (len(digest_info) + 1))
    plt.ylim(0, gel_height + 90)
    ax.invert_yaxis()
    plt.xticks([])
    plt.yticks([])
    plt.title("Virtual Restriction Digest Gel", color="white", pad=20)
    plt.tight_layout()
    plt.show()


dna_sequence = "TCAATTATATTCAGCATGGAAAGAATAAAAGAACTACGGAATCTAATGTCGCAGTCTCGCACCCGCGAGATACTAACAAAAACCACAGTGGACCATATGGCCATAATTAAGAAGTACACATCGGGGAGACAGGAAAAGAACCCGTCACTTAGAATGAAATGGATGATGGCAATGAAATATCCAATTACTGCTGACAAAAGGATAACAGAAATGGTTCCAGAGAGAAATGAACAAGGACAAACCCTATGGAGTAAAATGAGTGATGCTGGGTCAGATAGAGTGATGGTATCACCTTTGGCTGTAACATGGTGGAATAGAAATGGGCCCGTGACAAATACGGTCCATTACCCAAAAGTGTACAAAACTTATTTTGACAAAGTCGAAAGGTTGAAACATGGAACCTTCGGCCCTGTCCATTTTAGAAACCAAGTCAAAATACGTAGAAGAGTAGACACAAACCCTGGTCATGCAGACCTCAGTGCCAAAGAGGCACAAGATGTAATTATGGAAGTTGTTTTTCCCAATGAAGTGGGGGCCAGAATACTAACATCAGAATCACAGCTAACAATAACCAAAGAGAAAAAAGAAGAACTCCGAGATTGCAAAATTTCCCCCTTGATGGTCGCATACATGCTAGAGAGAGAACTTGTGCGGAAAACAAGATTTCTCCCAGTTGCTGGCGGAACAAGCAGTATATACATTGAAGTTTTACATTTGACTCAAGGAACGTGTTGGGAACAAATGTACACTCCAGGTGGAGGAGTGAGGAATGACGATGTTGACCAAAGCCTAATTATTGCGGCCAGGAACATAGTGAGAAGAGCCGCAGTGTCAGCAGATCCACTCGCATCTTTATTGGAGATGTGCCACAGCACGCAAATTGGCGGAACAAGGATGGTGGACATTCTTAGGCAGAACCCGACTGAAGAACAAGCTGTGGATATATGCAAAGCTGCAATGGGATTGAGAATCAGCTCATCTTTCAGCTTTGGTGGCTTTACATTTAAAAGAACGAGCGGGTCGTCAGTCAAAAGAGATGAAGAGGTTCTTACAGGTAATCTCCAAACATTGAGAATAAGAGTACATGAGGGGTATGAGGAATTCACAATGGTGGGGAAAAGAGCAACAGCTATACTAAGAAAAGCAACCAGAAGACTGGTTCAACTCATAGTGAGTGGAAGAGACGAACAGTCAGTAGCCGAGGCAATAATCGTGGCCATGGTTTTTTCCCAAGAAGATTGCATGATAAAAGCAGTTAGAGGTGACCTGAATTTTGTCAACAGAGCAAATCAGCGGTTGAACCCCATGCATCAGCTTTTAAGGCATTTTCAGAAAGATGCGAAAGTACTCTTTCAAAATTGGGGAGTTGAACACATCGACAGTGTGATGGGAATGGTTGGAGTATTACCAGATATGACTCCAAGCACAGAGATGTCAATGAGAGGAATAAGAGTCAGCAAAATGGGCGTGGATGAATACTCCAGTACAGAGAGGGTGGTGGTTAGCATTGATAGGTTTTTGAGAGTTCGAGACCAACGGGGGAATGTATTGTTATCTCCTGAGGAAGTCAGTGAAACACAAGGAACTGAAAGACTGACCATAACTTATTCATCATCGATGATGTGGGAAATTAATGGGCCTGAGTCGGTTTTGGTCAATACCTATCAATGGATCATCAGGAATTGGGAAGCTATCAAAATTCAGTGGTCTCAGAACCCTGCAATGTTGTACAACAAAATGGAATTTGAACCATTTCAATCTTTAGTCCCCAAGGCCACTAGAAGCCAATACAGTGGGTTTGTCAGAACTCTATTCCAACAAATGAGAGACGTACTTGGGACATTTGACACTGCCCAGATAATAAAGCTTCTCCCTTTTGCAGCTGCTCCACCAAAGCAAAGCAGAATGCAGTTCTCTTCACTGACTGTGAATGTGAGGGGATCAGGGATGAGAATACTTGTAAGGGGCAATTCTCCTGTATTCAACTACAACAAGACCACTAAAAGGCTAACAATTCTTGGAAAAGATGCCGGCACTTTAATTGAAGACCCAGATGAAAGCACATCCGGAGTGGAGTCCGCCGTCTTGAGAGGGTTCCTCATTATAGGTAAAGAAGACAGAAGATACGGACCAGCATTAAGCATCAATGAACTGAGTAACCTTGCAAAAGGGGAAAAGGCTAATGTGTTAATTGGGCAAGGAGACGTGGTGTTGGTAATGAAACGGAAACGGGACTCTAGTATACTTACTGACAGCCAGACAGCGACCAAACGAATTCGGATGGCCATCAATTAATATTGAATAGTTTAAAAACGA"

digest_data = compute_fragments(dna_sequence, ENZYMES)


print(f"DNA Sequence Length: {len(dna_sequence)} bp")
print(f"Source: Influenza A virus segment 4 (NCBI: NC_002016.1)")

for enzyme, result in digest_data.items():
    num_cuts = len(result["cuts"])
    num_fragments = len(result["sizes"])
    
    print(f"\n{enzyme} ({ENZYMES[enzyme]['site']}):")
    print(f"  Number of cleavages: {num_cuts}")
    
    if num_cuts > 0:
        print(f"  Cleavage positions: {result['cuts']}")
    else:
        print(f"  Cleavage positions: None found")
    
    print(f"  Number of fragments: {num_fragments}")
    print(f"  Fragment sizes (bp): {result['sizes']}")


draw_gel(digest_data)
