import math
import time
import urllib.parse
import urllib.request
from collections import Counter

import matplotlib.pyplot as plt

ALPHABET = "ACGT"

motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT",
]

def make_count_matrix(seqs):
    L = len(seqs[0])
    for s in seqs:
        if len(s) != L:
            raise ValueError("All motif sequences must have the same length.")
        if any(ch not in ALPHABET for ch in s):
            raise ValueError(f"Invalid character in motif: {s}")

    counts = {b: [0] * L for b in ALPHABET}
    for i in range(L):
        col = [s[i] for s in seqs]
        c = Counter(col)
        for b in ALPHABET:
            counts[b][i] = c.get(b, 0)
    return counts


def add_pseudocounts(counts, alpha=1.0):
    L = len(next(iter(counts.values())))
    return {b: [counts[b][i] + alpha for i in range(L)] for b in ALPHABET}


def relative_frequencies(weights):
    L = len(next(iter(weights.values())))
    col_sums = [sum(weights[b][i] for b in ALPHABET) for i in range(L)]
    return {b: [weights[b][i] / col_sums[i] for i in range(L)] for b in ALPHABET}


def log_likelihood_matrix(freqs, null_prob=0.25):
    L = len(next(iter(freqs.values())))
    loglik = {b: [0.0] * L for b in ALPHABET}
    for b in ALPHABET:
        for i in range(L):
            p = freqs[b][i]
            loglik[b][i] = math.log(p / null_prob)
    return loglik


def score_window(window, loglik):
    return sum(loglik[window[i]][i] for i in range(len(window)))


def scan_sequence(seq, L, loglik):
    seq = seq.upper()
    scores = []
    for start in range(len(seq) - L + 1):
        w = seq[start : start + L]
        if any(ch not in ALPHABET for ch in w):
            scores.append(float("nan"))
        else:
            scores.append(score_window(w, loglik))
    return scores


def find_top_peaks(scores, k=5, min_distance=200):
    indexed = [(i, s) for i, s in enumerate(scores) if not (isinstance(s, float) and math.isnan(s))]
    indexed.sort(key=lambda x: x[1], reverse=True)

    peaks = []
    for idx, s in indexed:
        if all(abs(idx - chosen) >= min_distance for chosen, _ in peaks):
            peaks.append((idx, s))
            if len(peaks) == k:
                break
    return peaks


def ncbi_esearch(db, term, retmax=10):
    params = {"db": db, "term": term, "retmax": str(retmax), "retmode": "json"}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url) as r:
        import json
        data = json.loads(r.read().decode("utf-8"))
    return data["esearchresult"]["idlist"]


def ncbi_efetch_fasta(db, ids):
    params = {"db": db, "id": ",".join(ids), "rettype": "fasta", "retmode": "text"}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url) as r:
        return r.read().decode("utf-8")


def parse_fasta(text):
    records = []
    header = None
    seq_parts = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def main():
    counts = make_count_matrix(motifs)
    weights = add_pseudocounts(counts, alpha=1.0)
    freqs = relative_frequencies(weights)
    loglik = log_likelihood_matrix(freqs, null_prob=0.25)
    motif_len = len(motifs[0])

    query = 'Influenza A virus[Organism] AND "complete genome"[Title]'
    ids = ncbi_esearch(db="nuccore", term=query, retmax=10)
    if not ids:
        raise RuntimeError("No genomes found. Try changing the query.")

    time.sleep(0.4)
    fasta_text = ncbi_efetch_fasta(db="nuccore", ids=ids)
    records = parse_fasta(fasta_text)

    if len(records) < 10:
        print(f"Warning: only got {len(records)} FASTA records.")

    records = records[:10]

    fig, axes = plt.subplots(5, 2, figsize=(14, 18))
    axes = axes.flatten()

    for i, (hdr, seq) in enumerate(records, start=1):
        ax = axes[i - 1]

        scores = scan_sequence(seq, motif_len, loglik)
        x = list(range(len(scores)))

        peaks = find_top_peaks(scores, k=5, min_distance=200)

        ax.plot(x, scores, linewidth=0.8, alpha=0.7, label="Motif Score Landscape")

        if peaks:
            px = [p[0] for p in peaks]
            py = [p[1] for p in peaks]
            ax.scatter(px, py, s=35, label="Top 5 Candidate Motifs")

        ax.set_title(f"Influenza Genome {i} - Exon-Intron Boundary Signals", fontsize=11)
        ax.set_xlabel("Genomic Position (bp)", fontsize=9)
        ax.set_ylabel("Log-likelihood Score", fontsize=9)

        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right", fontsize=8, frameon=True)

    for j in range(len(records), 10):
        axes[j].axis("off")

    plt.tight_layout()
    out = "influenza_10_genomes_motif_signals.png"
    plt.savefig(out, dpi=220)
    plt.close()

    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
