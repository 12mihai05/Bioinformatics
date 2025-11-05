# Gel electrophoresis is an analysis method implemented in all disciplines of life sciences. The results of gel electrophoresis indicate the relative sizes of fragments, which is useful for restriction mapping and analyzing PCR fragments. 

# 1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).

# 2. Take 10 random samples from this sequence, between 100-3000 bases.

# 3. Store these samples in an array.

# 4. Simulate the migration of these DNA segments on the electrophoresis gel, based on their molecular weights - however, their length should be sufficient for this exercise (show a visual representation).

# Note: Short DNA fragments meet small friction forces and travel faster through the electrophoresis gel. Long DNA fragments exhibit a high friction force and travel slowly through the electrophoresis gel.

# Gel electrophoresis simulation using a real NCBI DNA sequence
# Shorter fragments migrate farther (lower on the gel)

import math
import random
import matplotlib.pyplot as plt
import pandas as pd

# --------------------- Config ---------------------
SEED = 42
N_SAMPLES = 10
MIN_LEN, MAX_LEN = 100, 3000
LADDER_BP = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000,
             1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]
random.seed(SEED)

# --------------------- Sequence -------------------
dna = (
"TCAATTATATTCAGCATGGAAAGAATAAAAGAACTACGGAATCTAATGTCGCAGTCTCGCACCCGCGAGATACTAACAAAAACCACAGTGGACCATATGGCCATAATTAAGAAGTACACATCGGGGAGACAGGAAAAGAACCCGTCACTTAGAATGAAATGGATGATGGCAATGAAATATCCAATTACTGCTGACAAAAGGATAACAGAAATGGTTCCAGAGAGAAATGAACAAGGACAAACCCTATGGAGTAAAATGAGTGATGCTGGGTCAGATAGAGTGATGGTATCACCTTTGGCTGTAACATGGTGGAATAGAAATGGGCCCGTGACAAATACGGTCCATTACCCAAAAGTGTACAAAACTTATTTTGACAAAGTCGAAAGGTTGAAACATGGAACCTTCGGCCCTGTCCATTTTAGAAACCAAGTCAAAATACGTAGAAGAGTAGACACAAACCCTGGTCATGCAGACCTCAGTGCCAAAGAGGCACAAGATGTAATTATGGAAGTTGTTTTTCCCAATGAAGTGGGGGCCAGAATACTAACATCAGAATCACAGCTAACAATAACCAAAGAGAAAAAAGAAGAACTCCGAGATTGCAAAATTTCCCCCTTGATGGTCGCATACATGCTAGAGAGAGAACTTGTGCGGAAAACAAGATTTCTCCCAGTTGCTGGCGGAACAAGCAGTATATACATTGAAGTTTTACATTTGACTCAAGGAACGTGTTGGGAACAAATGTACACTCCAGGTGGAGGAGTGAGGAATGACGATGTTGACCAAAGCCTAATTATTGCGGCCAGGAACATAGTGAGAAGAGCCGCAGTGTCAGCAGATCCACTCGCATCTTTATTGGAGATGTGCCACAGCACGCAAATTGGCGGAACAAGGATGGTGGACATTCTTAGGCAGAACCCGACTGAAGAACAAGCTGTGGATATATGCAAAGCTGCAATGGGATTGAGAATCAGCTCATCTTTCAGCTTTGGTGGCTTTACATTTAAAAGAACGAGCGGGTCGTCAGTCAAAAGAGATGAAGAGGTTCTTACAGGTAATCTCCAAACATTGAGAATAAGAGTACATGAGGGGTATGAGGAATTCACAATGGTGGGGAAAAGAGCAACAGCTATACTAAGAAAAGCAACCAGAAGACTGGTTCAACTCATAGTGAGTGGAAGAGACGAACAGTCAGTAGCCGAGGCAATAATCGTGGCCATGGTTTTTTCCCAAGAAGATTGCATGATAAAAGCAGTTAGAGGTGACCTGAATTTTGTCAACAGAGCAAATCAGCGGTTGAACCCCATGCATCAGCTTTTAAGGCATTTTCAGAAAGATGCGAAAGTACTCTTTCAAAATTGGGGAGTTGAACACATCGACAGTGTGATGGGAATGGTTGGAGTATTACCAGATATGACTCCAAGCACAGAGATGTCAATGAGAGGAATAAGAGTCAGCAAAATGGGCGTGGATGAATACTCCAGTACAGAGAGGGTGGTGGTTAGCATTGATAGGTTTTTGAGAGTTCGAGACCAACGGGGGAATGTATTGTTATCTCCTGAGGAAGTCAGTGAAACACAAGGAACTGAAAGACTGACCATAACTTATTCATCATCGATGATGTGGGAAATTAATGGGCCTGAGTCGGTTTTGGTCAATACCTATCAATGGATCATCAGGAATTGGGAAGCTATCAAAATTCAGTGGTCTCAGAACCCTGCAATGTTGTACAACAAAATGGAATTTGAACCATTTCAATCTTTAGTCCCCAAGGCCACTAGAAGCCAATACAGTGGGTTTGTCAGAACTCTATTCCAACAAATGAGAGACGTACTTGGGACATTTGACACTGCCCAGATAATAAAGCTTCTCCCTTTTGCAGCTGCTCCACCAAAGCAAAGCAGAATGCAGTTCTCTTCACTGACTGTGAATGTGAGGGGATCAGGGATGAGAATACTTGTAAGGGGCAATTCTCCTGTATTCAACTACAACAAGACCACTAAAAGGCTAACAATTCTTGGAAAAGATGCCGGCACTTTAATTGAAGACCCAGATGAAAGCACATCCGGAGTGGAGTCCGCCGTCTTGAGAGGGTTCCTCATTATAGGTAAAGAAGACAGAAGATACGGACCAGCATTAAGCATCAATGAACTGAGTAACCTTGCAAAAGGGGAAAAGGCTAATGTGTTAATTGGGCAAGGAGACGTGGTGTTGGTAATGAAACGGAAACGGGACTCTAGTATACTTACTGACAGCCAGACAGCGACCAAACGAATTCGGATGGCCATCAATTAATATTGAATAGTTTAAAAACGA"
)

# --------------------- Sampling -------------------
def random_subseq(seq: str, min_len: int, max_len: int) -> dict:
    max_len = min(max_len, len(seq))
    k = random.randint(min_len, max_len)
    start = random.randint(0, len(seq) - k)
    return {"start": start, "length": k, "seq": seq[start:start + k]}

samples = [random_subseq(dna, MIN_LEN, MAX_LEN) for _ in range(N_SAMPLES)]

# --------------------- Checks & table --------------
N = len(dna)
for s in samples:
    assert 0 <= s["start"] <= N - s["length"]
    assert MIN_LEN <= s["length"] <= min(MAX_LEN, N)

df = pd.DataFrame([{"Sample": i + 1, "Start": s["start"], "Length (bp)": s["length"]}
                   for i, s in enumerate(samples)])
print(df.to_string(index=False))

# --------------------- Gel mapping -----------------
all_bp = LADDER_BP + [s["length"] for s in samples]
min_bp, max_bp = min(all_bp), max(all_bp)

GEL_H, GEL_W = 700, 400
LANE_W, LANE_GAP = 40, 120

def bp_to_y(bp: int, top=60, bottom=40) -> float:
    log_min, log_max = math.log10(min_bp), math.log10(max_bp)
    frac = (math.log10(bp) - log_min) / (log_max - log_min)
    return top + (GEL_H - top - bottom) * (1 - frac)

# --------------------- Plot ------------------------
fig = plt.figure(figsize=(5, 9), dpi=120)
ax = plt.gca()
ax.set_xlim(0, GEL_W); ax.set_ylim(GEL_H, 0); ax.axis("off")
ax.add_patch(plt.Rectangle((0, 0), GEL_W, GEL_H, facecolor="black"))

# Ladder lane
lane_x = 60
for bp in LADDER_BP:
    ax.add_patch(plt.Rectangle((lane_x, bp_to_y(bp) - 2), LANE_W, 4, facecolor="white"))
for ref in (3000, 1500, 500):
    ax.text(10, bp_to_y(ref) + 4, f"{ref} bp", color="white", fontsize=8, va="top")
ax.text(lane_x + LANE_W / 2, 25, "Ladder", color="white", fontsize=10, ha="center")

# Sample lane + inline labels (with light de-overlap)
lane2_x = lane_x + LANE_W + LANE_GAP
bands = sorted([(bp_to_y(s["length"]), s["length"]) for s in samples], key=lambda x: x[0])

# draw bands
for y, bp in bands:
    ax.add_patch(plt.Rectangle((lane2_x, y - 2.5), LANE_W, 5, facecolor="white"))

# label bands next to each line
label_x = lane2_x + LANE_W + 8
min_gap = 14  # pixels; push labels so they don't collide
y_text_prev = -1e9
for y, bp in bands:
    y_text = max(y, y_text_prev + min_gap)
    ax.text(label_x, y_text, f"{bp} bp", color="white", fontsize=9, va="center", ha="left")
    y_text_prev = y_text

ax.text(lane2_x + LANE_W / 2, 25, "Sample mix", color="white", fontsize=10, ha="center")

plt.tight_layout()
plt.savefig("gel_simulation.png", bbox_inches="tight", facecolor="black")
plt.show()
print("\nSaved: gel_simulation.png")
