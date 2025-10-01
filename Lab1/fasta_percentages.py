#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from collections import Counter
from pathlib import Path

# --- CONFIG: hardcoded default path (Windows-safe raw string) ---
DEFAULT_FASTA_PATH = r"C:\Users\pasar\Desktop\Facultate\Bioinformatics\Lab1\my_sequence.fasta"

# ---------- FASTA utilities ----------
def parse_fasta(path: str):
    """Yield (header, sequence_string) from a FASTA file."""
    header, seq_chunks = None, []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header, seq_chunks = line[1:].strip(), []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)

def clean_seq(s: str) -> str:
    return "".join(c for c in s.upper() if not c.isspace())

def percentages(seq: str):
    """Return list of (symbol, count, percent) sorted by symbol."""
    seq = clean_seq(seq)
    if not seq:
        return [], 0
    total = len(seq)
    counts = Counter(seq)
    items = [(sym, cnt, (cnt / total) * 100.0) for sym, cnt in sorted(counts.items())]
    return items, total

# ---------- GUI ----------
class FastaApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Alphabet Percentages")
        self.geometry("720x520")
        self.minsize(680, 480)

        # Data holders
        self.records = []           # list[(header, sequence)]
        self.concat_seq = ""        # combined sequence

        # Top bar: path entry + buttons
        top = ttk.Frame(self, padding=8)
        top.pack(fill="x")

        ttk.Label(top, text="FASTA path:").pack(side="left")
        self.path_var = tk.StringVar(value=DEFAULT_FASTA_PATH)
        self.path_entry = ttk.Entry(top, textvariable=self.path_var)
        self.path_entry.pack(side="left", fill="x", expand=True, padx=6)

        ttk.Button(top, text="Browse…", command=self.browse).pack(side="left", padx=2)
        ttk.Button(top, text="Load", command=self.load_fasta).pack(side="left")

        # Middle: record list (left) + table (right)
        mid = ttk.Frame(self, padding=8)
        mid.pack(fill="both", expand=True)

        # Record list
        left = ttk.Frame(mid)
        left.pack(side="left", fill="y")
        ttk.Label(left, text="Records").pack(anchor="w", pady=(0, 4))
        self.record_list = tk.Listbox(left, height=12, exportselection=False)
        self.record_list.pack(fill="y")
        self.record_list.bind("<<ListboxSelect>>", self.on_record_select)

        # Table
        right = ttk.Frame(mid)
        right.pack(side="left", fill="both", expand=True, padx=(10, 0))

        self.info_var = tk.StringVar(value="Length: -    Alphabet: -")
        ttk.Label(right, textvariable=self.info_var).pack(anchor="w", pady=(0, 6))

        cols = ("symbol", "count", "percent")
        self.tree = ttk.Treeview(right, columns=cols, show="headings", height=14)
        self.tree.heading("symbol", text="Symbol")
        self.tree.heading("count", text="Count")
        self.tree.heading("percent", text="Percent (%)")
        self.tree.column("symbol", width=90, anchor="center")
        self.tree.column("count", width=100, anchor="center")
        self.tree.column("percent", width=120, anchor="center")
        self.tree.pack(fill="both", expand=True)

        # Status bar
        self.status = tk.StringVar(value="Ready")
        ttk.Label(self, textvariable=self.status, anchor="w", padding=6).pack(fill="x")

        # Try autoload on startup (if file exists)
        self.after(100, self.autoload_default)

    # --- actions ---
    def autoload_default(self):
        p = Path(self.path_var.get())
        if p.exists() and p.is_file():
            self.load_fasta()

    def browse(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.faa"), ("All files", "*.*")]
        )
        if path:
            self.path_var.set(path)

    def load_fasta(self):
        path = self.path_var.get().strip()
        try:
            path_obj = Path(path)
            if not path_obj.exists():
                raise FileNotFoundError(f"File not found:\n{path}")
            self.records = list(parse_fasta(path))
            if not self.records:
                messagebox.showinfo("FASTA", "No FASTA records found.")
                self.clear_views()
                return

            # Build combined sequence
            self.concat_seq = "".join(seq for _, seq in self.records)

            # Populate record list
            self.record_list.delete(0, tk.END)
            self.record_list.insert(tk.END, "Overall (all records)")
            for i, (hdr, _) in enumerate(self.records, 1):
                # Trim very long headers for display
                display = hdr if len(hdr) <= 60 else hdr[:57] + "..."
                self.record_list.insert(tk.END, f"{i}. {display}")
            self.record_list.selection_clear(0, tk.END)
            self.record_list.selection_set(0)  # select Overall by default

            # Render table for Overall
            self.render_table(self.concat_seq, label="Overall")
            self.status.set(f"Loaded {len(self.records)} record(s) from: {path_obj.name}")

        except Exception as e:
            self.status.set("Error")
            messagebox.showerror("Error", str(e))

    def on_record_select(self, _evt):
        if not self.records:
            return
        idx = self.record_list.curselection()
        if not idx:
            return
        i = idx[0]
        if i == 0:
            self.render_table(self.concat_seq, label="Overall")
        else:
            hdr, seq = self.records[i - 1]
            self.render_table(seq, label=hdr)

    def render_table(self, seq: str, label: str):
        # Clear table
        for row in self.tree.get_children():
            self.tree.delete(row)

        rows, total = percentages(seq)
        # Update info label
        alphabet = "".join(sym for sym, _, _ in rows)
        self.info_var.set(f"Length: {total}    Alphabet: {{{', '.join(alphabet)}}}")

        # Fill rows
        for sym, cnt, pct in rows:
            self.tree.insert("", "end", values=(sym, cnt, f"{pct:.2f}"))

        self.status.set(f"Showing: {label}")

    def clear_views(self):
        self.record_list.delete(0, tk.END)
        for row in self.tree.get_children():
            self.tree.delete(row)
        self.info_var.set("Length: -    Alphabet: -")

# ---------- run ----------
if __name__ == "__main__":
    # nice default ttk theme if available
    try:
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)  # crisp on Windows (optional)
    except Exception:
        pass

    app = FastaApp()
    app.mainloop()
