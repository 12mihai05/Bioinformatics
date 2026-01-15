# 1. Use a random DNA sequence of about 50 letters. Use the sequence to compute the transtion probabilities between letters. 
#    The output should be the transition matrix stored as a json file.

import random
import json
from collections import defaultdict

DNA_BASES = ['A', 'T', 'G', 'C']

random_sequence = ''.join(random.choices(DNA_BASES, k=50))
print(f"Generated DNA sequence ({len(random_sequence)} letters):")
print(random_sequence)
print()

transition_counts = defaultdict(lambda: defaultdict(int))

for i in range(len(random_sequence) - 1):
    current_base = random_sequence[i]
    next_base = random_sequence[i + 1]
    transition_counts[current_base][next_base] += 1

transition_matrix = {}
for from_base in DNA_BASES:
    total_transitions = sum(transition_counts[from_base].values())
    transition_matrix[from_base] = {}
    
    for to_base in DNA_BASES:
        if total_transitions > 0:
            probability = transition_counts[from_base][to_base] / total_transitions
        else:
            probability = 0.0
        transition_matrix[from_base][to_base] = probability

print("Transition Probability Matrix:")
print("From \\ To", end="")
for base in DNA_BASES:
    print(f"\t{base}", end="")
print()

for from_base in DNA_BASES:
    print(f"{from_base}", end="")
    for to_base in DNA_BASES:
        prob = transition_matrix[from_base][to_base]
        print(f"\t{prob:.4f}", end="")
    print()

print()

output_file = "transition_matrix.json"
with open(output_file, 'w') as f:
    json.dump(transition_matrix, f, indent=2)

print(f"Transition matrix saved to '{output_file}'")
