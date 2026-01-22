# We are insde a courtroom, Mihai is accused of plagiarism. The lawywer must prove his innocence, and asks you for your help as a specialist in informations. You take 2 poetries: one of Mihai Eminescu and the other of Nichita Stanescu. You use these 2 texts to create 2 transition matrices, which is able to captiure the probability of transition between  words. You combine the matrices into a log likelihood matrix. We use the LLM to scan, by using a sliding window, the text of Mihai, which has been acused of plagiarism. Negative scores -> one model, positive scores -> another model, zero values -> neither. Emulate the text of Mihai using ChatGPT and ask it to make it a combination between 2 poetrues of Eminescu and Stanescu. At the end of the process you should be able to tell which parts of the text is written by Eminescu or by Stanescu or by neither.

import math
import re
from collections import defaultdict, Counter
from typing import List, Tuple, Dict

EMINESCU_TEXT = r"""
Somnoroase păsărele
Pe la cuiburi se adună,
Se ascund în rămurele —
Noapte bună!
Lunea plină-n zarea lin
Varsă razele-n păduri,
Iar vântul șoapte-n codri plouă
Iar flori'n codrii se visează
Iar codrul lin se-nclină sub vânt
Iar pasul lin pe frunzele cade
Iar vântul șoapte-n codri cântă
Iar noaptea lin se-ntinde'n vale
Iar linul vântului s-apropie
"""

STANESCU_TEXT = r"""
Leoaică tânără, iubirea
mi-a sărit în față.
Mă pândise-n încordare
mai demult.
Mi-a cerut o inimă
și o trecere prin foc,
iar eu nu aveam
decât cuvinte și tăcere.
Plouă în codrul tău
plouă-n cuvintele mele
plouă-n ochii mei de om
plouă peste tot ce-i sfânt
"""

MIHAI_TEXT = r"""
Somnoroase păsărele în codrii se adună,
Leoaică tânără, iubirea, pe cărări se-ntinde
Se ascund în rămurele și noaptea cântă lin
Iubirea mi-a sărit în grădini de plouă
Lunea plină varsă razele în păduri
Mă pândise-n încordare cu cuvinte și tăcere
Pasul lin pe frunzele noaptei se-nclină
Iar vântul șoapte prin codrul se-nchină
"""

def tokenize(text: str) -> List[str]:
    text = text.lower().strip()
    tokens = re.findall(r"[a-zăâîșţț]+", text, flags=re.IGNORECASE)
    return [t for t in tokens if t]

def bigrams(tokens: List[str]) -> List[Tuple[str, str]]:
    if len(tokens) < 2:
        return []
    return list(zip(tokens[:-1], tokens[1:]))

def build_transition_counts(tokens: List[str]) -> Dict[str, Counter]:
    trans = defaultdict(Counter)
    for a, b in bigrams(tokens):
        trans[a][b] += 1
    return trans

def build_vocabulary(*token_lists: List[str]) -> List[str]:
    vocab = set()
    for toks in token_lists:
        vocab.update(toks)
    vocab.add("<START>")
    vocab.add("<END>")
    vocab.add("<UNK>")
    return sorted(vocab)

def cond_prob(trans_counts: Dict[str, Counter], prev: str, nxt: str,
              vocab_size: int, alpha: float = 1.0) -> float:

    out_counts = trans_counts.get(prev)
    if out_counts is None:
        return alpha / (alpha * vocab_size)
    
    numerator = out_counts.get(nxt, 0) + alpha
    total_out = sum(out_counts.values())
    denom = total_out + alpha * vocab_size
    
    if denom == 0:
        return alpha / (alpha * vocab_size)
    
    return numerator / denom

def build_llr_lookup(em_tokens: List[str], st_tokens: List[str], alpha: float = 1.0):
    em_counts = build_transition_counts(em_tokens)
    st_counts = build_transition_counts(st_tokens)
    vocab = build_vocabulary(em_tokens, st_tokens)
    V = len(vocab)

    def llr(prev: str, nxt: str) -> float:
        prev_safe = prev if prev in vocab else "<UNK>"
        nxt_safe = nxt if nxt in vocab else "<UNK>"
        
        pe = cond_prob(em_counts, prev_safe, nxt_safe, V, alpha=alpha)
        ps = cond_prob(st_counts, prev_safe, nxt_safe, V, alpha=alpha)
        
        pe = max(pe, 1e-10)
        ps = max(ps, 1e-10)
        
        return math.log(pe) - math.log(ps)

    return llr, vocab

def score_tokens_llr(tokens: List[str], llr_func) -> float:
    if len(tokens) < 2:
        return 0.0
    s = 0.0
    count = 0
    for a, b in bigrams(tokens):
        s += llr_func(a, b)
        count += 1
    
    return s / count if count > 0 else 0.0

def sliding_windows(tokens: List[str], window_size: int, step: int):
    if not tokens:
        return
    
    if len(tokens) <= window_size:
        yield 0, tokens
    else:
        for i in range(0, len(tokens) - window_size + 1, step):
            yield i, tokens[i:i+window_size]

def label_score(score: float, eps: float = 0.25) -> str:
    if score > eps:
        return "EMINESCU-like (positive)"
    if score < -eps:
        return "STANESCU-like (negative)"
    return "NEITHER / ambiguous (near zero)"

def token_level_attribution(tokens: List[str], window_size: int, step: int, llr_func, eps: float = 0.25):
    if not tokens:
        return [], []
    
    votes = [0.0] * len(tokens)
    counts = [0] * len(tokens)

    for start, w in sliding_windows(tokens, window_size, step):
        sc = score_tokens_llr(w, llr_func)
        win_len = len(w)
        for j in range(start, min(start + win_len, len(tokens))):
            votes[j] += sc
            counts[j] += 1

    avg = []
    for i in range(len(tokens)):
        if counts[i] > 0:
            a = votes[i] / counts[i]
        else:
            a = 0.0
        avg.append(a)

    labels = []
    for a in avg:
        if a > eps:
            labels.append("E")
        elif a < -eps:
            labels.append("S")
        else:
            labels.append("N")
    return avg, labels

def compress_spans(tokens: List[str], labels: List[str], min_span_len: int = 1) -> List[Tuple[str, str]]:
    if not tokens:
        return []
    
    spans = []
    cur_lab = labels[0]
    cur_words = [tokens[0]]

    for t, lab in zip(tokens[1:], labels[1:]):
        if lab == cur_lab:
            cur_words.append(t)
        else:
            if len(cur_words) >= min_span_len:
                spans.append((cur_lab, " ".join(cur_words)))
            cur_lab = lab
            cur_words = [t]
    
    if len(cur_words) >= min_span_len:
        spans.append((cur_lab, " ".join(cur_words)))
    
    return spans

if __name__ == "__main__":
    em_tokens = tokenize(EMINESCU_TEXT)
    st_tokens = tokenize(STANESCU_TEXT)
    mi_tokens = tokenize(MIHAI_TEXT)

    llr_func, vocab = build_llr_lookup(em_tokens, st_tokens, alpha=0.5)

    mihai_len = len(mi_tokens)
    if mihai_len < 15:
        WINDOW_SIZE = max(3, mihai_len // 3)
        STEP = 1
    elif mihai_len < 30:
        WINDOW_SIZE = 5
        STEP = 1
    elif mihai_len < 60:
        WINDOW_SIZE = 7
        STEP = 2
    else:
        WINDOW_SIZE = 10
        STEP = 3

    EPS = 0.25

    print("=" * 80)
    print("PLAGIARISM DETECTION ANALYSIS")
    print("=" * 80)


    avg_scores, tok_labels = token_level_attribution(
        mi_tokens, WINDOW_SIZE, STEP, llr_func, eps=EPS
    )

    legend = {"E": "EMINESCU", "S": "STANESCU", "N": "NEITHER/Mixed"}
    spans = compress_spans(mi_tokens, tok_labels, min_span_len=1)

    if spans:
        for lab, text in spans:
            print(f"[{legend[lab]}]: {text}\n")
    else:
        print("No significant spans found.")

    total_tokens = len(tok_labels) if tok_labels else 0
    if total_tokens > 0:
        e_pct = tok_labels.count("E") * 100.0 / total_tokens
        s_pct = tok_labels.count("S") * 100.0 / total_tokens
        n_pct = tok_labels.count("N") * 100.0 / total_tokens
        print("=" * 80)
        print("OVERALL STYLE PERCENTAGES")
        print("=" * 80)
        print(f"Eminescu : {e_pct:5.1f}%")
        print(f"Stanescu  : {s_pct:5.1f}%")
        print(f"Neither   : {n_pct:5.1f}%")
