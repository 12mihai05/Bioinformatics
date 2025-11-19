# Make an artificial DNA sequence of 200-400 bases in length in which to determine 3-4 transposable elements
import random

BASES = "ACGT"

TRANSPOSABLE_ELEMENTS = {
    "TE1": "ACGTTGCA",   
    "TE2": "GTTG",      
    "TE3": "TGCATTT",    
    "TE4": "CTAGCTAG"    
}

def random_dna(length: int) -> str:
    return "".join(random.choice(BASES) for _ in range(length))

def print_sequence_with_indices(seq: str, line_length: int = 60):
    for i in range(0, len(seq), line_length):
        chunk = seq[i:i+line_length]
        print(f"{i:4d}: {chunk}")

def build_dna() -> tuple[str, int]:
    total_len = random.randint(200, 400)

    prefix_len = random.randint(50, 80)
    prefix = random_dna(prefix_len)

    core_region = TRANSPOSABLE_ELEMENTS["TE1"] + "TTT"

   
    middle = core_region + TRANSPOSABLE_ELEMENTS["TE4"]

    remaining = total_len - (len(prefix) + len(middle))
    if remaining < 0:
        remaining = 0  
    suffix = random_dna(remaining)

    dna = prefix + middle + suffix
    te1_start = len(prefix) 
    return dna, te1_start

def detect_elements(dna: str, elements: dict) -> list[dict]:
    detections = []
    for name, seq in elements.items():
        start_pos = 0
        while True:
            pos = dna.find(seq, start_pos)
            if pos == -1:
                break
            detections.append({
                "name": name,
                "sequence": seq,
                "start": pos,
                "end": pos + len(seq)
            })
            start_pos = pos + 1
    return detections

def main():
    dna, te1_start = build_dna()

    detections = detect_elements(dna, TRANSPOSABLE_ELEMENTS)

    print("=== DNA SEQUENCE ===")
    print(f"Length: {len(dna)} bases\n")
    print_sequence_with_indices(dna)

    print("\n=== ALL DETECTED TRANSPOSABLE ELEMENTS ===")
    for d in detections:
        print(f"{d['name']} at {d['start']}–{d['end']-1} : {d['sequence']}")

    te1_list = [d for d in detections if d["name"] == "TE1"]
    if not te1_list:
        print("\n[Unexpected] No TE1 found.")
        return

    te1 = te1_list[0]
    s1, e1 = te1["start"], te1["end"]
    print(f"\nWe focus on TE1 at {s1}–{e1-1} : {te1['sequence']}")

    print("\nCase 1: TE2 completely inside TE1")
    found_nested = False
    for d in detections:
        if d["name"] == "TE2":
            if d["start"] >= s1 and d["end"] <= e1:
                print(f"  -> TE2 at {d['start']}–{d['end']-1} is inside TE1")
                found_nested = True
    if not found_nested:
        print("  (No TE2 found inside TE1 this time, which would be unexpected.)")

    print("\nCase 2: TE3 starts inside TE1 and ends outside it")
    found_overlap = False
    for d in detections:
        if d["name"] == "TE3":
            if d["start"] >= s1 and d["start"] < e1 and d["end"] > e1:
                print(f"  -> TE3 at {d['start']}–{d['end']-1} starts in TE1 and ends after TE1")
                found_overlap = True
    if not found_overlap:
        print("  (No TE3 with partial overlap found, which would be unexpected.)")

    print("\nCase 3: TE4 separate (non-overlapping with TE1)")
    for d in detections:
        if d["name"] == "TE4":
            if d["end"] <= s1 or d["start"] >= e1:
                print(f"  -> TE4 at {d['start']}–{d['end']-1} does not overlap TE1")

if __name__ == "__main__":
    main()
