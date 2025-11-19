# Implement a software app which is able to detect the problem of this trasnposable element (start, end)within the created DNA sequence.

import random

BASES = "ACGT"

TRANSPOSABLE_ELEMENTS = {
    "TE1": "ACGTTGCA",
    "TE2": "TTAAGGCC",
    "TE3": "GGGTTTAA",
    "TE4": "CTAGCTAG"
}

def random_dna(length: int) -> str:
    return "".join(random.choice(BASES) for _ in range(length))

def insert_elements(dna: str, elements_to_insert: dict):
    dna_list = list(dna)
    length = len(dna_list)
    inserts = []

    def overlaps(a_start, a_end, b_start, b_end):
        return not (a_end <= b_start or b_end <= a_start)

    used_intervals = []

    for name, seq in elements_to_insert.items():
        te_len = len(seq)

        while True:
            start = random.randint(0, length - te_len)
            end = start + te_len
            if all(not overlaps(start, end, u_start, u_end) for u_start, u_end in used_intervals):
                break

        dna_list[start:end] = list(seq)
        used_intervals.append((start, end))
        inserts.append({
            "name": name,
            "sequence": seq,
            "start": start,
            "end": end
        })

    return "".join(dna_list), inserts

def detect_elements(dna: str, elements: dict):
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

def analyze_transposable_problems(detections):
    problems = []
    n = len(detections)

    for i in range(n):
        for j in range(i + 1, n):
            a = detections[i]
            b = detections[j]
            s1, e1 = a["start"], a["end"]
            s2, e2 = b["start"], b["end"]

            if e1 <= s2 or e2 <= s1:
                continue

            if s1 <= s2 and e2 <= e1:
                problems.append(
                    f"Problem: {b['name']} ({s2}-{e2-1}) is completely inside "
                    f"{a['name']} ({s1}-{e1-1})."
                )
            elif s2 <= s1 and e1 <= e2:
                problems.append(
                    f"Problem: {a['name']} ({s1}-{e1-1}) is completely inside "
                    f"{b['name']} ({s2}-{e2-1})."
                )
            else:
                if s1 < s2 < e1 < e2:
                    problems.append(
                        f"Problem: {b['name']} ({s2}-{e2-1}) starts inside "
                        f"{a['name']} ({s1}-{e1-1}) and finishes outside it."
                    )
                elif s2 < s1 < e2 < e1:
                    problems.append(
                        f"Problem: {a['name']} ({s1}-{e1-1}) starts inside "
                        f"{b['name']} ({s2}-{e2-1}) and finishes outside it."
                    )
                else:
                    problems.append(
                        f"Problem: {a['name']} ({s1}-{e1-1}) and "
                        f"{b['name']} ({s2}-{e2-1}) overlap in a complex way."
                    )

    return problems

def print_sequence_with_indices(seq: str, line_length: int = 50):
    for i in range(0, len(seq), line_length):
        chunk = seq[i:i+line_length]
        print(f"{i:4d}: {chunk}")

def main():
    length = random.randint(200, 400)
    dna = random_dna(length)

    k = random.randint(3, 4)
    chosen_names = random.sample(list(TRANSPOSABLE_ELEMENTS.keys()), k=k)
    chosen_elements = {name: TRANSPOSABLE_ELEMENTS[name] for name in chosen_names}

    final_dna, true_inserts = insert_elements(dna, chosen_elements)


    outer_start = 50
    outer_end = 80
    outer_seq = final_dna[outer_start:outer_end]

    inner_start = outer_start + 5
    inner_end = outer_start + 15
    inner_seq = final_dna[inner_start:inner_end]

    overlap_start = outer_start + 10
    overlap_end = outer_end + 5
    overlap_seq = final_dna[overlap_start:overlap_end]

    extra_elements = {
        "TE_OUTER": outer_seq,
        "TE_INNER": inner_seq,
        "TE_OVERLAP": overlap_seq,
    }

    print("\n=== EXTRA TEs FOR SPECIAL USE CASES ===")
    print(f"TE_OUTER    start={outer_start}, end={outer_end-1}")
    print(f"TE_INNER    start={inner_start}, end={inner_end-1}   (inside TE_OUTER)")
    print(f"TE_OVERLAP  start={overlap_start}, end={overlap_end-1} (starts inside, ends outside TE_OUTER)")

    all_elements = {**chosen_elements, **extra_elements}
    detections = detect_elements(final_dna, all_elements)

    print("\n=== CREATED DNA SEQUENCE ===")
    print(f"Length: {len(final_dna)} bases\n")
    print_sequence_with_indices(final_dna)

    print("\n=== TRUE INSERT POSITIONS (GROUND TRUTH, NON-OVERLAPPING INSERTS) ===")
    for ins in true_inserts:
        print(
            f"  {ins['name']}  start={ins['start']}, end={ins['end']-1}, "
            f"seq={ins['sequence']}"
        )

    print("\n=== DETECTED TRANSPOSABLE ELEMENTS (INCLUDING SPECIAL TEs) ===")
    for det in detections:
        print(
            f"  {det['name']}  start={det['start']}, end={det['end']-1}, "
            f"seq={det['sequence']}"
        )

    print("\n=== CHECK: DID WE FIND ALL TRUE INSERTS? ===")
    true_set = {(i["name"], i["start"], i["end"]) for i in true_inserts}
    det_set_original = {(d["name"], d["start"], d["end"]) for d in detections if d["name"] in chosen_elements}

    missing = true_set - det_set_original
    extra   = det_set_original - true_set

    if not missing and not extra:
        print("OK ✅  All inserted original TEs were correctly detected.")
    else:
        if missing:
            print("Missing detections (inserts not found):")
            for m in missing:
                print(" ", m)
        if extra:
            print("Extra detections (found but not in ground truth):")
            for e in extra:
                print(" ", e)

    print("\n=== PROBLEM ANALYSIS BETWEEN TRANSPOSABLE ELEMENTS (start, end) ===")
    problems = analyze_transposable_problems(detections)
    if not problems:
        print("No problems detected: no TE is inside or partially overlapping another TE.")
    else:
        for p in problems:
            print(p)

if __name__ == "__main__":
    main()
