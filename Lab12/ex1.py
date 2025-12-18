import math
from collections import Counter

motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTA ACT".replace(" ", ""),  # "ATGGTAACT"
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT",
]

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

ALPHABET = "ACGT"

def make_count_matrix(seqs):
    L = len(seqs[0])
    for s in seqs:
        if len(s) != L:
            raise ValueError("All motif sequences must have the same length.")
        if any(ch not in ALPHABET for ch in s):
            raise ValueError(f"Invalid character in motif: {s}")

    counts = {b: [0]*L for b in ALPHABET}
    for i in range(L):
        col = [s[i] for s in seqs]
        c = Counter(col)
        for b in ALPHABET:
            counts[b][i] = c.get(b, 0)
    return counts

def add_pseudocounts(counts, alpha=1.0):
    L = len(next(iter(counts.values())))
    weights = {b: [counts[b][i] + alpha for i in range(L)] for b in ALPHABET}
    return weights

def relative_frequencies(weights):
    L = len(next(iter(weights.values())))
    col_sums = [sum(weights[b][i] for b in ALPHABET) for i in range(L)]
    freqs = {b: [weights[b][i] / col_sums[i] for i in range(L)] for b in ALPHABET}
    return freqs

def log_likelihood_matrix(freqs, null_prob=0.25, log_base="e"):
    L = len(next(iter(freqs.values())))
    loglik = {b: [0.0]*L for b in ALPHABET}

    def log_fn(x):
        if log_base == "e":
            return math.log(x)
        elif log_base == 2:
            return math.log(x, 2)
        elif log_base == 10:
            return math.log10(x)
        else:
            raise ValueError("log_base must be 'e', 2, or 10")

    for b in ALPHABET:
        for i in range(L):
            p = freqs[b][i]
            loglik[b][i] = log_fn(p / null_prob)
    return loglik

def score_window(window, loglik):
    return sum(loglik[window[i]][i] for i in range(len(window)))

def scan_sequence(S, L, loglik):
    scores = []
    for start in range(len(S) - L + 1):
        w = S[start:start+L]
        if any(ch not in ALPHABET for ch in w):
            continue
        scores.append((start, w, score_window(w, loglik)))
    return scores

def print_matrix(name, mat, digits=3):
    print(f"\n{name}")
    header = "pos  " + "  ".join(f"{i+1:>5}" for i in range(len(next(iter(mat.values())))))
    print(header)
    for b in ALPHABET:
        row = " ".join(f"{v:>5.{digits}f}" if isinstance(v, float) else f"{v:>5d}" for v in mat[b])
        print(f"{b:>3}  {row}")

if __name__ == "__main__":
    L = len(motifs[0])
    N = len(motifs)

    counts = make_count_matrix(motifs)

    alpha = 1.0
    weights = add_pseudocounts(counts, alpha=alpha)

    freqs = relative_frequencies(weights)

    loglik = log_likelihood_matrix(freqs, null_prob=0.25, log_base="e")

    print(f"Motif length L={L}, number of sequences N={N}, pseudocount alpha={alpha}")
    print_matrix("1) Count matrix", counts, digits=0)
    print_matrix("2) Weight matrix (counts + alpha)", weights, digits=0)
    print_matrix("3) Relative frequencies matrix", freqs, digits=3)
    print_matrix("4) Log-likelihood matrix (ln(P/0.25))", loglik, digits=3)

    scores = scan_sequence(S, L, loglik)
    scores_sorted = sorted(scores, key=lambda x: x[2], reverse=True)

    print("\n5) Sliding window scores on S")
    for start, w, sc in scores:
        print(f"start={start:>2}  window={w}  score={sc:.3f}")

    top_k = 5
    print(f"\nTop {top_k} windows (highest scores):")
    for start, w, sc in scores_sorted[:top_k]:
        print(f"start={start:>2}  window={w}  score={sc:.3f}")
