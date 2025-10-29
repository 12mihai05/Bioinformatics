# Take an arbitrary DNA seq. from the NCBI betweeen 1000 - 3000 nucleotides (letters)
# a) Take 2000 random samples from this seq of about 100-150 bases
# b) Store thus samples inside an array or list
# c) Rebuild the original DNA seq. by using only these random samples
# d) Make a text file called answer.txt in which you explain what will be the main issue with the algorithm approach

import random
import textwrap
from typing import List, Tuple, Optional, Dict, Set

DNA_SEQ = ("TCAATTATATTCAGCATGGAAAGAATAAAAGAACTACGGAATCTAATGTCGCAGTCTCGCACCCGCGAGATACTAACAAAAACCACAGTGGACCATATGGCCATAATTAAGAAGTACACATCGGGGAGACAGGAAAAGAACCCGTCACTTAGAATGAAATGGATGATGGCAATGAAATATCCAATTACTGCTGACAAAAGGATAACAGAAATGGTTCCAGAGAGAAATGAACAAGGACAAACCCTATGGAGTAAAATGAGTGATGCTGGGTCAGATAGAGTGATGGTATCACCTTTGGCTGTAACATGGTGGAATAGAAATGGGCCCGTGACAAATACGGTCCATTACCCAAAAGTGTACAAAACTTATTTTGACAAAGTCGAAAGGTTGAAACATGGAACCTTCGGCCCTGTCCATTTTAGAAACCAAGTCAAAATACGTAGAAGAGTAGACACAAACCCTGGTCATGCAGACCTCAGTGCCAAAGAGGCACAAGATGTAATTATGGAAGTTGTTTTTCCCAATGAAGTGGGGGCCAGAATACTAACATCAGAATCACAGCTAACAATAACCAAAGAGAAAAAAGAAGAACTCCGAGATTGCAAAATTTCCCCCTTGATGGTCGCATACATGCTAGAGAGAGAACTTGTGCGGAAAACAAGATTTCTCCCAGTTGCTGGCGGAACAAGCAGTATATACATTGAAGTTTTACATTTGACTCAAGGAACGTGTTGGGAACAAATGTACACTCCAGGTGGAGGAGTGAGGAATGACGATGTTGACCAAAGCCTAATTATTGCGGCCAGGAACATAGTGAGAAGAGCCGCAGTGTCAGCAGATCCACTCGCATCTTTATTGGAGATGTGCCACAGCACGCAAATTGGCGGAACAAGGATGGTGGACATTCTTAGGCAGAACCCGACTGAAGAACAAGCTGTGGATATATGCAAAGCTGCAATGGGATTGAGAATCAGCTCATCTTTCAGCTTTGGTGGCTTTACATTTAAAAGAACGAGCGGGTCGTCAGTCAAAAGAGATGAAGAGGTTCTTACAGGTAATCTCCAAACATTGAGAATAAGAGTACATGAGGGGTATGAGGAATTCACAATGGTGGGGAAAAGAGCAACAGCTATACTAAGAAAAGCAACCAGAAGACTGGTTCAACTCATAGTGAGTGGAAGAGACGAACAGTCAGTAGCCGAGGCAATAATCGTGGCCATGGTTTTTTCCCAAGAAGATTGCATGATAAAAGCAGTTAGAGGTGACCTGAATTTTGTCAACAGAGCAAATCAGCGGTTGAACCCCATGCATCAGCTTTTAAGGCATTTTCAGAAAGATGCGAAAGTACTCTTTCAAAATTGGGGAGTTGAACACATCGACAGTGTGATGGGAATGGTTGGAGTATTACCAGATATGACTCCAAGCACAGAGATGTCAATGAGAGGAATAAGAGTCAGCAAAATGGGCGTGGATGAATACTCCAGTACAGAGAGGGTGGTGGTTAGCATTGATAGGTTTTTGAGAGTTCGAGACCAACGGGGGAATGTATTGTTATCTCCTGAGGAAGTCAGTGAAACACAAGGAACTGAAAGACTGACCATAACTTATTCATCATCGATGATGTGGGAAATTAATGGGCCTGAGTCGGTTTTGGTCAATACCTATCAATGGATCATCAGGAATTGGGAAGCTATCAAAATTCAGTGGTCTCAGAACCCTGCAATGTTGTACAACAAAATGGAATTTGAACCATTTCAATCTTTAGTCCCCAAGGCCACTAGAAGCCAATACAGTGGGTTTGTCAGAACTCTATTCCAACAAATGAGAGACGTACTTGGGACATTTGACACTGCCCAGATAATAAAGCTTCTCCCTTTTGCAGCTGCTCCACCAAAGCAAAGCAGAATGCAGTTCTCTTCACTGACTGTGAATGTGAGGGGATCAGGGATGAGAATACTTGTAAGGGGCAATTCTCCTGTATTCAACTACAACAAGACCACTAAAAGGCTAACAATTCTTGGAAAAGATGCCGGCACTTTAATTGAAGACCCAGATGAAAGCACATCCGGAGTGGAGTCCGCCGTCTTGAGAGGGTTCCTCATTATAGGTAAAGAAGACAGAAGATACGGACCAGCATTAAGCATCAATGAACTGAGTAACCTTGCAAAAGGGGAAAAGGCTAATGTGTTAATTGGGCAAGGAGACGTGGTGTTGGTAATGAAACGGAAACGGGACTCTAGTATACTTACTGACAGCCAGACAGCGACCAAACGAATTCGGATGGCCATCAATTAATATTGAATAGTTTAAAAACGA")

N_READS = 2000
READ_LEN_RANGE: Tuple[int, int] = (100, 150)
MIN_OVERLAP = 20
SEED = 42
random.seed(SEED)

def sample_reads(seq: str, n_reads: int, read_len_range: Tuple[int, int]) -> List[str]:
    reads: List[str] = []
    lo, hi = read_len_range
    L = len(seq)
    if hi > L:
        raise ValueError(f"Read upper bound {hi} exceeds sequence length {L}.")
    while len(reads) < n_reads:
        k = random.randint(lo, hi)
        start = random.randint(0, L - k)
        r = seq[start:start + k]
        if len(r) == k:
            reads.append(r)
    return reads

def validate_and_report_reads(reads: List[str], lo: int, hi: int, label: str = "reads") -> None:
    assert len(reads) == 2000, f"{label} must have 2000 items, got {len(reads)}."
    bad = [i for i, r in enumerate(reads) if not (lo <= len(r) <= hi)]
    assert not bad, f"{label} contains items out of length range [{lo},{hi}] at indices {bad[:5]}..."
    from collections import Counter
    lengths = Counter(map(len, reads))
    smallest = min(lengths)
    largest = max(lengths)
    print(f"{label}: OK (n={len(reads)}). Lengths range {smallest}..{largest}.")

def longest_overlap(a: str, b: str, min_olap: int) -> int:
    max_olap = min(len(a), len(b))
    if a[-min_olap:] != b[:min_olap]:
        return 0
    for olen in range(max_olap, min_olap - 1, -1):
        if a[-olen:] == b[:olen]:
            return olen
    return 0

def build_prefix_index(reads: List[Optional[str]], k: int) -> Dict[str, Set[int]]:
    idx: Dict[str, Set[int]] = {}
    for i, r in enumerate(reads):
        if not r or len(r) < k:
            continue
        p = r[:k]
        idx.setdefault(p, set()).add(i)
    return idx

def update_prefix_index(idx: Dict[str, Set[int]], i: int, old: Optional[str], new: Optional[str], k: int) -> None:
    if old and len(old) >= k:
        p_old = old[:k]
        s = idx.get(p_old)
        if s and i in s:
            s.remove(i)
            if not s:
                idx.pop(p_old, None)
    if new and len(new) >= k:
        p_new = new[:k]
        idx.setdefault(p_new, set()).add(i)

def greedy_assemble_indexed(reads_in: List[str], k: int, min_olap: int) -> str:
    reads: List[Optional[str]] = list(reads_in)
    idx = build_prefix_index(reads, k)
    made_progress = True
    while made_progress:
        made_progress = False
        for i, r in enumerate(reads):
            if not r or len(r) < k:
                continue
            suf = r[-k:]
            cand_idxs = idx.get(suf, set())
            if not cand_idxs:
                continue
            best_j = None
            best_ol = 0
            for j in list(cand_idxs):
                if j == i:
                    continue
                rj = reads[j]
                if not rj:
                    continue
                ol = longest_overlap(r, rj, min_olap)
                if ol > best_ol:
                    best_ol = ol
                    best_j = j
            if best_j is None or best_ol < min_olap:
                continue
            old_i = reads[i]
            old_j = reads[best_j]
            merged = r + reads[best_j][best_ol:]
            reads[i] = merged
            reads[best_j] = None
            update_prefix_index(idx, i, old_i, merged, k)
            update_prefix_index(idx, best_j, old_j, None, k)
            made_progress = True
    contigs = [r for r in reads if r]
    contigs.sort(key=len, reverse=True)
    return "".join(contigs)

def write_answer_file(content: str, filename: str = "answer.txt") -> None:
    with open(filename, "w", encoding="utf-8") as f:
        f.write(content)

def main():
    seq = "".join(c for c in DNA_SEQ.upper() if c in "ACGTN")
    L = len(seq)
    assert 1000 <= L <= 3000, f"Provided sequence length {L} not in [1000, 3000]."
    print(f"Using provided sequence of length {L} nt.")
    reads = sample_reads(seq, N_READS, READ_LEN_RANGE)
    validate_and_report_reads(reads, READ_LEN_RANGE[0], READ_LEN_RANGE[1])
    print("Generated 2000 reads of 100–150 bases from within the sequence.")
    assert len(reads) == 2000, "Did not generate 2000 reads."
    assert all(READ_LEN_RANGE[0] <= len(r) <= READ_LEN_RANGE[1] for r in reads), "Read out of 100–150 bp range."
    assert all(r in seq for r in reads), "A sampled read is not a substring of the original sequence."
    print("Generated 2000 reads of 100–150 bases from within the sequence.")
    assembled = greedy_assemble_indexed(reads, k=MIN_OVERLAP, min_olap=MIN_OVERLAP)
    print(f"Assembled sequence length: {len(assembled)} nt")
    exact_match = (assembled == seq)
    contains_original = seq in assembled
    assembled_contains = assembled in seq
    explanation = f"""
    DNA Assembly from Random Samples — Report

    Provided sequence length: {L} nt
    Random reads: {N_READS}
    Read length range: {READ_LEN_RANGE[0]}–{READ_LEN_RANGE[1]} bases
    Minimum overlap (k): {MIN_OVERLAP} bases
    Assembled length: {len(assembled)} nt

    Reconstruction check:
    - Exact match with original: {exact_match}
    - Assembled contains original: {contains_original}
    - Original contains assembled: {assembled_contains}

    Main issue with this algorithmic approach:
    Even with 2000 reads of 100–150 bp from the original sequence, overlaps are frequently
    ambiguous in repetitive regions and some positions may lack sufficient overlapping coverage.
    A greedy, local-overlap strategy makes irreversible choices, which can misjoin or fragment contigs.
    Robust assembly typically relies on coverage-aware, graph-based methods (de Bruijn/overlap graphs)
    with error models; this simple greedy approach is not reliable for exact reconstruction.
    """
    write_answer_file(textwrap.dedent(explanation).strip())
    print('Wrote "answer.txt" explaining the main issue.')

if __name__ == "__main__":
    main()
