# 2. Use a random english text of about 300 letters (that implies spaces, commas, dots). Compute the transition probabilities between words.
#    Store the transition matrix as a json file. For ease of implementation, you can represent each new word by using a symbol of your coice (ASCII).

import random
import json
from collections import defaultdict
import re

word_list = [
    "the", "quick", "brown", "fox", "jumps", "over", "lazy", "dog", "and",
    "cat", "runs", "through", "forest", "mountains", "river", "sky", "sun",
    "moon", "stars", "night", "day", "morning", "evening", "birds", "tree",
    "water", "fire", "wind", "earth", "world", "life", "time", "people",
    "place", "house", "road", "city", "country", "sea", "beach", "island"
]

text = []
char_count = 0
sentence_length = 0

while char_count < 300:
    if sentence_length == 0:
        sentence_length = random.randint(5, 12)
    
    word = random.choice(word_list)
    text.append(word)
    char_count += len(word) + 1
    sentence_length -= 1
    
    if sentence_length == 0:
        punctuation = random.choice(['.', '.', '.', '!', '?'])  
        text.append(punctuation)
        char_count += 1
        
        if random.random() < 0.3 and len(text) > 3:
            text.insert(-1, ',')

random_text = ' '.join(text)
random_text = re.sub(r'\s+([.,!?])', r'\1', random_text)

print("Generated English Text:")
print(random_text)
print(f"\nText length: {len(random_text)} characters")
print()

words = re.findall(r'\b[a-z]+\b', random_text.lower())
print(f"Total words: {len(words)}")
print()

word_to_symbol = {}
symbol_to_word = {}
symbols = []

for i, word in enumerate(set(words)):
    if i < 256:  
        symbol = chr(65 + i)  
        word_to_symbol[word] = symbol
        symbol_to_word[symbol] = word

print("Word-to-Symbol Mapping:")
for word, symbol in sorted(word_to_symbol.items()):
    print(f"  {word:15} -> {symbol}")
print()

symbol_sequence = [word_to_symbol[word] for word in words]
print(f"Symbol sequence: {''.join(symbol_sequence)}")
print()

reconstructed_words = [symbol_to_word[symbol] for symbol in symbol_sequence]
reconstructed_text = ' '.join(reconstructed_words)
print("Reconstructed text from symbols:")
print(reconstructed_text)
print()

transition_counts = defaultdict(lambda: defaultdict(int))

for i in range(len(words) - 1):
    current_word = words[i]
    next_word = words[i + 1]
    transition_counts[current_word][next_word] += 1

transition_matrix = {}
for word in sorted(set(words)):
    total_transitions = sum(transition_counts[word].values())
    transition_matrix[word] = {}
    
    for next_word in sorted(set(words)):
        if total_transitions > 0:
            probability = transition_counts[word][next_word] / total_transitions
        else:
            probability = 0.0
        
        if probability > 0:
            transition_matrix[word][next_word] = round(probability, 6)

print("Sample Transition Probabilities (words with transitions):")
for word in sorted(transition_matrix.keys())[:5]:
    print(f"  From '{word}': ", end="")
    if transition_matrix[word]:
        for next_word, prob in sorted(transition_matrix[word].items()):
            print(f"'{next_word}'={prob:.4f} ", end="")
    print()

print()

output_file = "word_transition_matrix.json"
with open(output_file, 'w') as f:
    json.dump(transition_matrix, f, indent=2)

mapping_file = "word_symbol_mapping.json"
with open(mapping_file, 'w') as f:
    json.dump(word_to_symbol, f, indent=2)

print(f"Word transition matrix saved to '{output_file}'")
print(f"Word-to-symbol mapping saved to '{mapping_file}'")
