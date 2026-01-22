import math
from collections import defaultdict

ALPHABET = ["A", "C", "G", "T"]
IDX = {b: i for i, b in enumerate(ALPHABET)}

def count_transitions(seq: str):
    seq = seq.strip().upper()
    counts = {f: {t: 0 for t in ALPHABET} for f in ALPHABET}
    totals = {f: 0 for f in ALPHABET}

    for a, b in zip(seq, seq[1:]):
        if a in IDX and b in IDX:
            counts[a][b] += 1
            totals[a] += 1

    return counts, totals

def transition_probabilities(counts, totals):
    Tr = {f: {t: 0.0 for t in ALPHABET} for f in ALPHABET}
    for f in ALPHABET:
        if totals[f] == 0:
            continue
        for t in ALPHABET:
            Tr[f][t] = counts[f][t] / totals[f]
    return Tr

def log2(x: float) -> float:
    return math.log(x) / math.log(2)

def log_likelihood_matrix(Tr_plus, Tr_minus):
    beta = {f: {t: 0.0 for t in ALPHABET} for f in ALPHABET}
    for f in ALPHABET:
        for t in ALPHABET:
            p = Tr_plus[f][t]
            q = Tr_minus[f][t]
            if p == 0.0 and q == 0.0:
                beta[f][t] = 0.0
            elif q == 0.0 and p > 0.0:
                beta[f][t] = float("inf")
            elif p == 0.0 and q > 0.0:
                beta[f][t] = float("-inf")
            else:
                beta[f][t] = log2(p / q)
    return beta

def score_sequence(seq: str, beta):
    seq = seq.strip().upper()
    score = 0.0
    for a, b in zip(seq, seq[1:]):
        if a in IDX and b in IDX:
            score += beta[a][b]
    return score

def print_matrix(title, M, digits=3):
    print(title)
    header = "     " + "  ".join(f"{c:>7}" for c in ALPHABET)
    print(header)
    for r in ALPHABET:
        row = []
        for c in ALPHABET:
            v = M[r][c]
            if v == float("inf"):
                s = "  +inf "
            elif v == float("-inf"):
                s = "  -inf "
            else:
                s = f"{v:7.{digits}f}"
            row.append(s)
        print(f"{r:>3}  " + "  ".join(row))
    print()

S1 = "ATCGATTCGATATCATACACGTAT"   
S2 = "CTCGACTAGTATGAAGTCCACGCTTG" 
S  = "CAGGTTGGAAACGTAA"           

counts_p, totals_p = count_transitions(S1)
Tr_plus = transition_probabilities(counts_p, totals_p)

counts_m, totals_m = count_transitions(S2)
Tr_minus = transition_probabilities(counts_m, totals_m)

beta = log_likelihood_matrix(Tr_plus, Tr_minus)

print_matrix("Transition matrix Tr+ (CpG island, +):", Tr_plus, digits=2)
print_matrix("Transition matrix Tr- (Non-island, -):", Tr_minus, digits=2)
print_matrix("Log-likelihood matrix beta = log2(Tr+/Tr-):", beta, digits=3)

ll_score = score_sequence(S, beta)
print(f"Sequence S = {S}")
print(f"log-likelihood score (sum beta over transitions) = {ll_score:.4f}")

if ll_score > 0:
    print("Decision: CpG ISLAND (+) is more likely.")
elif ll_score < 0:
    print("Decision: NON-ISLAND (-) is more likely.")
else:
    print("Decision: Tie / no preference (score = 0).")
