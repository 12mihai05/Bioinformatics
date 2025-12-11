import tkinter as tk
from tkinter import ttk, messagebox

class AlignmentApp:
    def __init__(self, root):
        self.root = root
        self.root.title("DNA Sequence Aligner (Needleman-Wunsch)")
        self.root.geometry("1100x700")
        
        self.default_s1 = "ACCGTGAAGCCAATAC"
        self.default_s2 = "AGCGTGCAGCCAATAC"
        
        self.top_frame = tk.Frame(root)
        self.top_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.left_panel = tk.LabelFrame(self.top_frame, text="Controls & Inputs", width=250)
        self.left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=5)
        self.left_panel.pack_propagate(False)

        self.mid_panel = tk.LabelFrame(self.top_frame, text="Matrix Heatmap")
        self.mid_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)
        
        self.right_panel = tk.LabelFrame(self.top_frame, text="Traceback Path")
        self.right_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)

        self.bottom_frame = tk.LabelFrame(root, text="Show Alignment:", height=200)
        self.bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=10)
        self.bottom_frame.pack_propagate(False)

        tk.Label(self.left_panel, text="Sequence 1:").pack(anchor="w", padx=5)
        self.entry_s1 = tk.Entry(self.left_panel)
        self.entry_s1.insert(0, self.default_s1)
        self.entry_s1.pack(fill=tk.X, padx=5, pady=(0, 10))

        tk.Label(self.left_panel, text="Sequence 2:").pack(anchor="w", padx=5)
        self.entry_s2 = tk.Entry(self.left_panel)
        self.entry_s2.insert(0, self.default_s2)
        self.entry_s2.pack(fill=tk.X, padx=5, pady=(0, 10))

        param_frame = tk.Frame(self.left_panel)
        param_frame.pack(fill=tk.X, padx=5)

        tk.Label(param_frame, text="Gap").grid(row=0, column=0, sticky="w")
        self.spin_gap = tk.Spinbox(param_frame, from_=-10, to=10, width=5)
        self.spin_gap.delete(0, "end")
        self.spin_gap.insert(0, 0)
        self.spin_gap.grid(row=0, column=1, pady=2)

        tk.Label(param_frame, text="Match").grid(row=1, column=0, sticky="w")
        self.spin_match = tk.Spinbox(param_frame, from_=-10, to=10, width=5)
        self.spin_match.delete(0, "end")
        self.spin_match.insert(0, 1) 
        self.spin_match.grid(row=1, column=1, pady=2)

        tk.Label(param_frame, text="Mismatch").grid(row=2, column=0, sticky="w")
        self.spin_mismatch = tk.Spinbox(param_frame, from_=-10, to=10, width=5)
        self.spin_mismatch.delete(0, "end")
        self.spin_mismatch.insert(0, -1)
        self.spin_mismatch.grid(row=2, column=1, pady=2)

        self.btn_align = tk.Button(self.left_panel, text="Align", command=self.run_alignment, height=2, bg="#ddd")
        self.btn_align.pack(fill=tk.X, padx=10, pady=20)

        tk.Label(self.left_panel, text="Presets").pack(anchor="w", padx=5)
        preset_frame = tk.Frame(self.left_panel)
        preset_frame.pack(fill=tk.BOTH, expand=True, padx=5)
        for i in range(16):
            btn = tk.Button(preset_frame, text=f"Setting {i+1}", font=("Arial", 7))
            btn.grid(row=i//2, column=i%2, sticky="ew", padx=1, pady=1)

        self.heatmap_canvas = tk.Canvas(self.mid_panel, bg="black")
        self.heatmap_canvas.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.grid_canvas = tk.Canvas(self.right_panel, bg="#FFFFE0") 
        self.grid_canvas.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.txt_output = tk.Text(self.bottom_frame, font=("Courier New", 10))
        self.txt_output.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.root.after(100, self.run_alignment)

    def get_color_gradient(self, value, min_val, max_val):
        if max_val == min_val:
            ratio = 0.5
        else:
            ratio = (value - min_val) / (max_val - min_val)
        red = int(255 * ratio)
        blue = int(100 * (1 - ratio))
        return f'#{red:02x}00{blue:02x}'

    def run_alignment(self):
        seq1 = self.entry_s1.get().upper()
        seq2 = self.entry_s2.get().upper()
        try:
            match_score = int(self.spin_match.get())
            mismatch_score = int(self.spin_mismatch.get())
            gap_penalty = int(self.spin_gap.get())
        except ValueError:
            messagebox.showerror("Error", "Parameters must be integers.")
            return

        n = len(seq1)
        m = len(seq2)

        score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
        
        for i in range(n + 1): score_matrix[i][0] = i * gap_penalty
        for j in range(m + 1): score_matrix[0][j] = j * gap_penalty

        min_score = 0
        max_score = 0
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                if seq1[i - 1] == seq2[j - 1]:
                    diag = score_matrix[i - 1][j - 1] + match_score
                else:
                    diag = score_matrix[i - 1][j - 1] + mismatch_score
                
                up = score_matrix[i - 1][j] + gap_penalty
                left = score_matrix[i][j - 1] + gap_penalty
                
                val = max(diag, up, left)
                score_matrix[i][j] = val
                
                if val > max_score: max_score = val
                if val < min_score: min_score = val

        align1 = ""
        align2 = ""
        i, j = n, m
        matches_count = 0
        path_coords = set()
        path_coords.add((i, j))

        while i > 0 and j > 0:
            current = score_matrix[i][j]
            diag_score = score_matrix[i - 1][j - 1]
            up_score = score_matrix[i - 1][j]
            left_score = score_matrix[i][j - 1]
            
            is_match = seq1[i - 1] == seq2[j - 1]
            diag_step = match_score if is_match else mismatch_score
            
            check_diag = (current == diag_score + diag_step)
            check_up = (current == up_score + gap_penalty)
            check_left = (current == left_score + gap_penalty)
            
            if is_match and check_diag:
                align1 += seq1[i - 1]
                align2 += seq2[j - 1]
                matches_count += 1
                i -= 1
                j -= 1
            elif check_up:
                align1 += seq1[i - 1]
                align2 += "-"
                i -= 1
            elif check_left:
                align1 += "-"
                align2 += seq2[j - 1]
                j -= 1
            else:
                align1 += seq1[i - 1]
                align2 += seq2[j - 1]
                i -= 1
                j -= 1
            
            path_coords.add((i, j))

        while i > 0:
            align1 += seq1[i - 1]
            align2 += "-"
            i -= 1
            path_coords.add((i, j))
        while j > 0:
            align1 += "-"
            align2 += seq2[j - 1]
            j -= 1
            path_coords.add((i, j))

        align1 = align1[::-1]
        align2 = align2[::-1]

        self.draw_heatmap(score_matrix, n, m, min_score, max_score)
        self.draw_grid(path_coords, n, m)
        self.update_text_output(align1, align2, matches_count)

    def draw_heatmap(self, matrix, rows, cols, min_val, max_val):
        self.heatmap_canvas.delete("all")
        w = self.heatmap_canvas.winfo_width()
        h = self.heatmap_canvas.winfo_height()
        if w < 10: w, h = 300, 300
        
        cell_w = w / (cols + 1)
        cell_h = h / (rows + 1)

        for i in range(rows + 1):
            for j in range(cols + 1):
                color = self.get_color_gradient(matrix[i][j], min_val, max_val)
                x0, y0 = j * cell_w, i * cell_h
                x1, y1 = x0 + cell_w, y0 + cell_h
                self.heatmap_canvas.create_rectangle(x0, y0, x1, y1, fill=color, outline="")

    def draw_grid(self, path_set, rows, cols):
        self.grid_canvas.delete("all")
        w = self.grid_canvas.winfo_width()
        h = self.grid_canvas.winfo_height()
        if w < 10: w, h = 300, 300
        
        cell_w = w / (cols + 1)
        cell_h = h / (rows + 1)

        for r in range(rows + 1):
            for c in range(cols + 1):
                x0 = c * cell_w
                y0 = r * cell_h
                x1 = x0 + cell_w
                y1 = y0 + cell_h
                
                if (r, c) in path_set:
                    fill_color = "#CD2626"
                else:
                    fill_color = "#FFFACD"
                
                self.grid_canvas.create_rectangle(x0, y0, x1, y1, fill=fill_color, outline="black", width=1)

    def update_text_output(self, s1, s2, matches):
        connections = ""
        for c1, c2 in zip(s1, s2):
            connections += "|" if (c1 == c2 and c1 != "-") else " "
        
        length = len(s1)
        similarity = int((matches / length) * 100) if length > 0 else 0

        text = f"{s1}\n{connections}\n{s2}\n\n"
        text += f"Matches = {matches}\n"
        text += f"Length  = {length}\n"
        text += f"Similarity = {similarity} %\n"
        text += f"Tracing back: M[{len(s1)},{len(s2)}]"
        
        self.txt_output.delete("1.0", tk.END)
        self.txt_output.insert(tk.END, text)

if __name__ == "__main__":
    root = tk.Tk()
    app = AlignmentApp(root)
    root.mainloop()