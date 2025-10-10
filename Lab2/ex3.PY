# Design an application by using AI, which contains a GUI that allows the user to select a fasta file. The content
# of the fasta file should be analyzed by using a sliding window of 30 positions. The content for each sliding window 
# should be used in order to extract/compute the relative frequencies of the nucleoids formed in the aphabet of the 
# sequence . The output of the app should be a chart containing 4 signals, one signal for each symbol for the alphabet of the sequence.
# In order to plot this chart first, the vector that contains the relative frequencies for each symbol must be ccalculated
# And only then plotted

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


ALPHABET = ("A", "C", "G", "T")

def read_fasta(path: str) -> str:
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts).upper()

def sliding_window(seq: str, k: int):
    n = len(seq)
    for i in range(0, n - k + 1):
        yield seq[i:i+k]

def relative_freqs_for_window(win: str):
    L = len(win)
    if L == 0:
        return {b: 0.0 for b in ALPHABET}
    counts = {b: 0 for b in ALPHABET}
    for ch in win:
        if ch in counts:
            counts[ch] += 1
    return {b: counts[b] / L for b in ALPHABET}

def frequency_vectors(seq: str, k: int):
    vecs = {b: [] for b in ALPHABET}
    for win in sliding_window(seq, k):
        freqs = relative_freqs_for_window(win)
        for b in ALPHABET:
            vecs[b].append(freqs[b])
    return vecs


class App(ttk.Frame):
    def __init__(self, master):
        super().__init__(master)
        master.title("FASTA Sliding-Window Nucleotide Frequencies")
        master.geometry("900x600")
        self.pack(fill=tk.BOTH, expand=True)

        self.seq = ""
        self.freq_vecs = None

        self._build_controls()
        self._build_plot()

    def _build_controls(self):
        bar = ttk.Frame(self)
        bar.pack(fill=tk.X, padx=8, pady=6)

        self.path_var = tk.StringVar()
        ttk.Label(bar, text="FASTA:").pack(side=tk.LEFT)
        self.path_entry = ttk.Entry(bar, textvariable=self.path_var, width=60, state="readonly")
        self.path_entry.pack(side=tk.LEFT, padx=6)

        ttk.Button(bar, text="Browse…", command=self.on_browse).pack(side=tk.LEFT, padx=4)

        ttk.Label(bar, text="Window (k):").pack(side=tk.LEFT, padx=(12,4))
        self.k_var = tk.StringVar(value="30")
        self.k_entry = ttk.Entry(bar, textvariable=self.k_var, width=6)
        self.k_entry.pack(side=tk.LEFT)

        ttk.Button(bar, text="Analyze & Plot", command=self.on_analyze).pack(side=tk.LEFT, padx=8)

        self.status_var = tk.StringVar(value="Pick a FASTA file to begin.")
        ttk.Label(self, textvariable=self.status_var, anchor="w").pack(fill=tk.X, padx=8, pady=(0, 6))

    def _build_plot(self):
        self.fig = Figure(figsize=(6.5, 4.5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("Window index")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_title("A/C/G/T frequencies over sliding windows")

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)


    def on_browse(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa *.frn"), ("All files", "*.*")]
        )
        if path:
            self.path_var.set(path)
            self.status_var.set(f"Selected: {os.path.basename(path)}")

    def on_analyze(self):
        path = self.path_var.get().strip()
        if not path:
            messagebox.showwarning("No file", "Please select a FASTA file.")
            return
        try:
            k = int(self.k_var.get().strip())
            if k <= 0:
                raise ValueError
        except Exception:
            messagebox.showerror("Invalid window", "Window size (k) must be a positive integer.")
            return

        try:
            self.seq = read_fasta(path)
        except Exception as e:
            messagebox.showerror("Read error", f"Could not read FASTA:\n{e}")
            return

        if len(self.seq) < k:
            messagebox.showwarning(
                "Sequence too short",
                f"Sequence length ({len(self.seq)}) is smaller than window size {k}."
            )
            self._plot_empty()
            return

        self.status_var.set("Computing frequency vectors…")
        self.update_idletasks()
        self.freq_vecs = frequency_vectors(self.seq, k)

        self._plot_vectors(self.freq_vecs)
        self.status_var.set(f"Done: length={len(self.seq)} | windows={len(self.freq_vecs['A'])} | k={k}")


    def _plot_empty(self):
        self.ax.clear()
        self.ax.set_xlabel("Window index")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_title("A/C/G/T frequencies over sliding windows")
        self.canvas.draw_idle()

    def _plot_vectors(self, vecs):
        self.ax.clear()
        n = len(next(iter(vecs.values()))) if vecs else 0
        x = list(range(n))

        self.ax.plot(x, vecs["A"], label="A")
        self.ax.plot(x, vecs["C"], label="C")
        self.ax.plot(x, vecs["G"], label="G")
        self.ax.plot(x, vecs["T"], label="T")

        self.ax.set_xlabel("Window index")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_title("A/C/G/T frequencies (sliding window)")
        self.ax.set_ylim(0.0, 1.0) 
        self.ax.legend(loc="upper right")
        self.ax.grid(True, linestyle=":", linewidth=0.5)

        self.canvas.draw_idle()


def main():
    root = tk.Tk()
    app = App(root)
    root.mainloop()

if __name__ == "__main__":
    main()
