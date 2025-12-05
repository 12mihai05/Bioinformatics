import math
from collections import Counter
import matplotlib.pyplot as plt

S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
WINDOW = 30



def cg_total(seq: str) -> float:

    seq = seq.upper()
    cg = sum(1 for b in seq if b in "CG")
    return round(cg / len(seq) * 100, 2)


def cg_pattern(seq: str, w: int):

    seq = seq.upper()
    n = len(seq)

    total_cg = sum(1 for b in seq if b in "CG")
    cgtot = 100 * total_cg / n 

    xs = []
    for i in range(n - w + 1):
        win = seq[i:i + w]
        cg_sw = sum(1 for b in win if b in "CG")
        cgsw = cgtot * cg_sw / len(win)
        xs.append(cgsw)

    return round(cgtot, 2), xs



def kappa_ic_window(win: str) -> float:
    
    A = win.upper()
    m = len(A)
    N = m - 1
    T = 0.0

    for u in range(1, N + 1):
        B = A[u:N] 
        if not B:
            continue
        C = sum(1 for i in range(len(B)) if A[i] == B[i])
        T += (C / len(B)) * 100.0

    return T / N 


def kappa_ic_total(seq: str) -> float:
    return round(kappa_ic_window(seq), 2)


def kappa_pattern(seq: str, w: int):
    seq = seq.upper()
    n = len(seq)
    ys = []
    for i in range(n - w + 1):
        win = seq[i:i + w]
        ys.append(kappa_ic_window(win))
    return ys



def build_pattern(seq: str, w: int):
    cgtot, xs = cg_pattern(seq, w)
    ys = kappa_pattern(seq, w)

    if len(xs) == 0:
        raise ValueError(
            f"Sequence too short (len={len(seq)}). It must be at least {w} bases long."
        )

    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    kic_tot = kappa_ic_total(seq)

    return xs, ys, cx, cy, cgtot, kic_tot




def plot_pattern(xs, ys, cx, cy, title_prefix="Promoter"):
    plt.figure(figsize=(6, 5))
    plt.scatter(xs, ys, s=10)
    plt.xlabel("(C+G)% per window (CGSW)")
    plt.ylabel("Kappa IC per window")
    plt.title(f"{title_prefix} pattern")
    plt.grid(True)

    plt.figure(figsize=(5, 4))
    plt.scatter([cx], [cy], s=60, marker="x")
    plt.xlabel("Center (C+G)%")
    plt.ylabel("Center Kappa IC")
    plt.title(f"{title_prefix} center of weight")
    plt.grid(True)

    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    xs, ys, cx, cy, cgtot, kic_tot = build_pattern(S, WINDOW)

    print(f"Sequence length: {len(S)}")
    print(f"Global C+G% (CGTOT)    = {cgtot}")     
    print(f"Global Kappa IC        = {kic_tot}")  
    print(f"Center of weight X (CG) = {round(cx, 2)}")
    print(f"Center of weight Y (IC) = {round(cy, 2)}")

    plot_pattern(xs, ys, cx, cy, title_prefix="Test sequence")



