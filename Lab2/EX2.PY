# find in sequence S only the dinucleodites and trinucleodites that exist, without the use of the brute force engine.
# In order to achive the result one must verify each combination starting from the beggining of sequence S 
# Example: ABAA we have AB, BA, AA, ABA

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
from collections import Counter

def scan_kmers(sequence: str, k: int):
    n = len(sequence)
    total = max(n - k + 1, 0)
    counts = Counter()
    for i in range(total):
        kmer = sequence[i:i+k]
        counts[kmer] += 1
    return counts, total

def rel_freq(count: int, total_windows: int) -> float:
    return (count / total_windows * 100.0) if total_windows > 0 else 0.0


class KmerScannerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Observed Di/Trinucleotides (No Brute Force)")
        self.root.geometry("800x600")
        self.root.resizable(True, True)

        self.default_sequence = "ATTGTCCCAATCTGTTG"
        self._build_ui()

    def _build_ui(self):
        bar = ttk.Frame(self.root)
        bar.pack(fill=tk.X, padx=8, pady=6)

        ttk.Label(bar, text="DNA sequence:").pack(side=tk.LEFT)
        self.sequence_entry = ttk.Entry(bar, width=80)
        self.sequence_entry.pack(side=tk.LEFT, padx=6)
        self.sequence_entry.insert(0, self.default_sequence)

        ttk.Button(bar, text="Analyze", command=self.analyze_sequence).pack(side=tk.LEFT, padx=4)
        ttk.Button(bar, text="Clear", command=self.clear_results).pack(side=tk.LEFT, padx=4)

        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill=tk.BOTH, expand=True, padx=8, pady=6)

        self.tab2 = ttk.Frame(self.nb)
        self.nb.add(self.tab2, text="Observed dinucleotides (k=2)")
        self.k2_text = scrolledtext.ScrolledText(self.tab2, wrap=tk.NONE, font=("Consolas", 10))
        self.k2_text.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.tab3 = ttk.Frame(self.nb)
        self.nb.add(self.tab3, text="Observed trinucleotides (k=3)")
        self.k3_text = scrolledtext.ScrolledText(self.tab3, wrap=tk.NONE, font=("Consolas", 10))
        self.k3_text.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.tabS = ttk.Frame(self.nb)
        self.nb.add(self.tabS, text="Summary")
        self.sum_text = scrolledtext.ScrolledText(self.tabS, wrap=tk.WORD, font=("Consolas", 10))
        self.sum_text.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(self.root, textvariable=self.status_var, anchor="w").pack(fill=tk.X, padx=8, pady=(0,6))

    def clear_results(self):
        self.k2_text.delete("1.0", tk.END)
        self.k3_text.delete("1.0", tk.END)
        self.sum_text.delete("1.0", tk.END)
        self.status_var.set("Cleared")

    def analyze_sequence(self):
        s = self.sequence_entry.get().strip().upper()

        if not s:
            messagebox.showwarning("Input", "Please enter a DNA sequence.")
            return

        valid = set("ACGT")
        if not all(ch in valid for ch in s):
            messagebox.showerror("Invalid", "Sequence must contain only A, C, G, T.")
            return

        self.status_var.set("Scanning...")
        self.root.update_idletasks()

        self.clear_results()

        k2_counts, total2 = scan_kmers(s, 2)
        k3_counts, total3 = scan_kmers(s, 3)

        lines2 = []
        lines2.append(f"Sequence length: {len(s)}")
        lines2.append(f"Total windows (k=2): {total2}")
        lines2.append("")
        lines2.append(f"{'Dinuc':<6} {'Count':>5} {'%':>7}")
        lines2.append("-" * 22)

        k2_items = [(kmer, cnt) for kmer, cnt in k2_counts.items() if cnt > 0]
        k2_items.sort(key=lambda x: (-x[1], x[0]))

        for kmer, cnt in k2_items:
            pct = rel_freq(cnt, total2)
            lines2.append(f"{kmer:<6} {cnt:>5} {pct:>6.2f}")

        self.k2_text.insert("1.0", "\n".join(lines2))

        lines3 = []
        lines3.append(f"Sequence length: {len(s)}")
        lines3.append(f"Total windows (k=3): {total3}")
        lines3.append("")
        lines3.append(f"{'Trinuc':<7} {'Count':>5} {'%':>7}")
        lines3.append("-" * 24)

        k3_items = [(kmer, cnt) for kmer, cnt in k3_counts.items() if cnt > 0]
        k3_items.sort(key=lambda x: (-x[1], x[0]))

        for kmer, cnt in k3_items:
            pct = rel_freq(cnt, total3)
            lines3.append(f"{kmer:<7} {cnt:>5} {pct:>6.2f}")

        self.k3_text.insert("1.0", "\n".join(lines3))

        sum_lines = []
        sum_lines.append(f"Sequence: {s}")
        sum_lines.append(f"Length: {len(s)}")
        sum_lines.append("")
        sum_lines.append(f"Observed dinucleotides: {len(k2_items)} (out of {max(len(s)-1,0)} windows)")
        if k2_items:
            sum_lines.append("  First few: " + ", ".join(k for k, _ in k2_items[:8]))
        sum_lines.append("")
        sum_lines.append(f"Observed trinucleotides: {len(k3_items)} (out of {max(len(s)-2,0)} windows)")
        if k3_items:
            sum_lines.append("  First few: " + ", ".join(k for k, _ in k3_items[:8]))
        sum_lines.append("")
        sum_lines.append("Note: We only scan existing windows (no brute-force of all possibilities).")
        self.sum_text.insert("1.0", "\n".join(sum_lines))

        self.status_var.set("Done")


def main():
    root = tk.Tk()
    app = KmerScannerGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()