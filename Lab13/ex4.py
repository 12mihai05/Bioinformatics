# 4. Use the transition matrix from the JSON matrix in order to synthesize new sequences of text, 
#    based on the transition matrix.

import json
import random

with open("word_transition_matrix.json", 'r') as f:
    transition_matrix = json.load(f)

print("Loaded transition matrix with", len(transition_matrix), "unique words")
print()

def synthesize_text(transition_matrix, start_word=None, length=20):

    if start_word is None:
        start_word = random.choice(list(transition_matrix.keys()))
    
    sequence = [start_word]
    current_word = start_word
    
    for _ in range(length - 1):
        if current_word in transition_matrix and transition_matrix[current_word]:
            next_words = list(transition_matrix[current_word].keys())
            probabilities = list(transition_matrix[current_word].values())
            
            current_word = random.choices(next_words, weights=probabilities, k=1)[0]
            sequence.append(current_word)
        else:
            current_word = random.choice(list(transition_matrix.keys()))
            sequence.append(current_word)
    
    return ' '.join(sequence)

print("=" * 70)
print("SYNTHESIZED TEXT SEQUENCES")
print("=" * 70)
print()

print("Sequence 1 (15 words):")
text1 = synthesize_text(transition_matrix, length=15)
print(text1)
print()

print("Sequence 2 (25 words):")
text2 = synthesize_text(transition_matrix, length=25)
print(text2)
print()

print("Sequence 3 (40 words):")
text3 = synthesize_text(transition_matrix, length=40)
print(text3)
print()

print("=" * 70)
print("SEQUENCES WITH SPECIFIC STARTING WORDS")
print("=" * 70)
print()

common_words = ['the', 'and', 'cat', 'dog', 'quick', 'fox', 'world']
available_start_words = [w for w in common_words if w in transition_matrix]

if available_start_words:
    for start_word in available_start_words[:3]: 
        print(f"Starting with '{start_word}':")
        text = synthesize_text(transition_matrix, start_word=start_word, length=20)
        print(text)
        print()
else:
    for start_word in list(transition_matrix.keys())[:3]:
        print(f"Starting with '{start_word}':")
        text = synthesize_text(transition_matrix, start_word=start_word, length=20)
        print(text)
        print()

print("=" * 70)
print("SAVING GENERATED SEQUENCES")
print("=" * 70)
print()

output_sequences = []
for i in range(5):
    seq = synthesize_text(transition_matrix, length=30)
    output_sequences.append(f"Sequence {i+1}:\n{seq}\n")

with open("synthesized_sequences.txt", 'w') as f:
    f.write('\n'.join(output_sequences))

print("5 synthesized sequences (30 words each) saved to 'synthesized_sequences.txt'")
