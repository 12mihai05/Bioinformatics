# Make a brute force engine which is able to generate all the dinucleotides and trinucleotides combinations.
# Search for each combination inside sequence S and calculate their relative frequencies.
# S = "ATTGTCCCAATCTGTTG"

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox


def generate_dinucleotides():
    nucleotides = ['A', 'C', 'G', 'T']
    dinucleotides = []
    for n1 in nucleotides:
        for n2 in nucleotides:
            dinucleotides.append(n1 + n2)
    return dinucleotides

def generate_trinucleotides():
    nucleotides = ['A', 'C', 'G', 'T']
    trinucleotides = []
    for n1 in nucleotides:
        for n2 in nucleotides:
            for n3 in nucleotides:
                trinucleotides.append(n1 + n2 + n3)
    return trinucleotides

def count_occurrences(sequence, pattern):
    count = 0
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i:i+len(pattern)] == pattern:
            count += 1
    return count

def calculate_percentage(sequence, pattern):
    pattern_length = len(pattern)
    total_possible = len(sequence) - pattern_length + 1
    if total_possible <= 0:
        return 0.0
    count = count_occurrences(sequence, pattern)
    return (count / total_possible) * 100


class NucleotideAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("K-mer Analysis")
        self.root.geometry("800x600")
        self.root.resizable(True, True)

        self.default_sequence = "ATTGTCCCAATCTGTTG"
        self._build_ui()

    def _build_ui(self):
        top = ttk.Frame(self.root)
        top.pack(fill=tk.X, padx=8, pady=6)

        ttk.Label(top, text="DNA sequence:").pack(side=tk.LEFT)
        self.sequence_entry = ttk.Entry(top, width=80)
        self.sequence_entry.pack(side=tk.LEFT, padx=6)
        self.sequence_entry.insert(0, self.default_sequence)

        ttk.Button(top, text="Analyze", command=self.analyze_sequence).pack(side=tk.LEFT, padx=4)
        ttk.Button(top, text="Clear", command=self.clear_results).pack(side=tk.LEFT, padx=4)

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=8, pady=6)

        self.dinuc_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.dinuc_frame, text="Dinucleotides (k=2)")
        self.dinuc_text = scrolledtext.ScrolledText(self.dinuc_frame, wrap=tk.NONE, font=("Consolas", 10))
        self.dinuc_text.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.trinuc_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.trinuc_frame, text="Trinucleotides (k=3)")
        self.trinuc_text = scrolledtext.ScrolledText(self.trinuc_frame, wrap=tk.NONE, font=("Consolas", 10))
        self.trinuc_text.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.summary_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.summary_frame, text="Summary")
        self.summary_text = scrolledtext.ScrolledText(self.summary_frame, wrap=tk.WORD, font=("Consolas", 10))
        self.summary_text.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(self.root, textvariable=self.status_var, anchor="w").pack(fill=tk.X, padx=8, pady=(0,6))

    def clear_results(self):
        self.dinuc_text.delete(1.0, tk.END)
        self.trinuc_text.delete(1.0, tk.END)
        self.summary_text.delete(1.0, tk.END)
        self.status_var.set("Cleared")

    def analyze_sequence(self):
        sequence = self.sequence_entry.get().strip().upper()

        if not sequence:
            messagebox.showwarning("Input", "Please enter a DNA sequence.")
            return

        valid = set("ACGT")
        if not all(c in valid for c in sequence):
            messagebox.showerror("Invalid", "Sequence must contain only A, C, G, T.")
            return

        self.status_var.set("Analyzing...")
        self.root.update_idletasks()

        self.clear_results()
        self.analyze_dinucleotides(sequence)
        self.analyze_trinucleotides(sequence)
        self.generate_summary(sequence)
        self.status_var.set(f"Done (length={len(sequence)})")


    def analyze_dinucleotides(self, sequence):
        dinucleotides = generate_dinucleotides()

        dinuc_results = []
        for dinuc in dinucleotides:
            cnt = count_occurrences(sequence, dinuc)
            pct = calculate_percentage(sequence, dinuc)
            dinuc_results.append((dinuc, cnt, pct))

        dinuc_results.sort(key=lambda x: (-x[2], x[0]))

        lines = []
        lines.append(f"Sequence length: {len(sequence)}")
        lines.append(f"Windows (k=2): {max(len(sequence) - 1, 0)}")
        lines.append("")
        lines.append(f"{'Dinuc':<6} {'Count':>5} {'%':>7}")
        lines.append("-" * 22)

        for dinuc, cnt, pct in dinuc_results:
            lines.append(f"{dinuc:<6} {cnt:>5} {pct:>6.2f}")

        self.dinuc_text.insert("1.0", "\n".join(lines))
        self.dinuc_results = dinuc_results

    def analyze_trinucleotides(self, sequence):
        trinucleotides = generate_trinucleotides()

        trinuc_results = []
        for trinuc in trinucleotides:
            cnt = count_occurrences(sequence, trinuc)
            pct = calculate_percentage(sequence, trinuc)
            trinuc_results.append((trinuc, cnt, pct))

        trinuc_results.sort(key=lambda x: (-x[2], x[0]))

        lines = []
        lines.append(f"Sequence length: {len(sequence)}")
        lines.append(f"Windows (k=3): {max(len(sequence) - 2, 0)}")
        lines.append("")
        lines.append(f"{'Trinuc':<7} {'Count':>5} {'%':>7}")
        lines.append("-" * 24)

        for trinuc, cnt, pct in trinuc_results:
            lines.append(f"{trinuc:<7} {cnt:>5} {pct:>6.2f}")

        self.trinuc_text.insert("1.0", "\n".join(lines))
        self.trinuc_results = trinuc_results

    def generate_summary(self, sequence):
        lines = []
        lines.append(f"Sequence: {sequence}")
        lines.append(f"Length: {len(sequence)}")
        lines.append("")

        lines.append("Base composition:")
        for nuc in ['A', 'C', 'G', 'T']:
            c = sequence.count(nuc)
            p = (c / len(sequence)) * 100 if sequence else 0
            lines.append(f"  {nuc}: {c} ({p:.2f}%)")
        lines.append("")

        lines.append("Top 5 dinucleotides:")
        for i, (dinuc, cnt, pct) in enumerate(self.dinuc_results[:5], 1):
            lines.append(f"  {i}. {dinuc}  {cnt}  ({pct:.2f}%)")
        lines.append("")

        lines.append("Top 5 trinucleotides:")
        for i, (trinuc, cnt, pct) in enumerate(self.trinuc_results[:5], 1):
            lines.append(f"  {i}. {trinuc}  {cnt}  ({pct:.2f}%)")
        lines.append("")

        self.summary_text.insert("1.0", "\n".join(lines))


def main():
    root = tk.Tk()
    app = NucleotideAnalyzerGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
