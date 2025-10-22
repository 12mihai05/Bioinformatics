#!/usr/bin/env python3

import matplotlib.pyplot as plt
import json
from pathlib import Path

def load_codon_analysis_results():
    """Load results from codon analysis."""
    results_file = "codon_analysis_results.json"

    if not Path(results_file).exists():
        print(f"❌ Could not find {results_file}")
        print("Please run codon_analysis.py first to generate the required data.")
        return None

    try:
        with open(results_file, 'r') as f:
            results = json.load(f)
        return results
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"❌ Error reading {results_file}: {e}")
        return None

def analyze_viral_amino_acids():
    """Analyze amino acids most frequent in viruses for dietary recommendations."""

    print("🧬 VIRAL AMINO ACID ANALYSIS & NUTRITIONAL RECOMMENDATIONS")
    print("=" * 70)

    # Load results from codon analysis
    results = load_codon_analysis_results()
    if not results:
        return None, None

    print(f"\n📈 ANALYSIS METADATA:")
    print("-" * 50)
    metadata = results['analysis_metadata']
    print(f"• COVID-19 genome length: {metadata['covid_genome_length']:,} bp")
    print(f"• Influenza genome length: {metadata['influenza_genome_length']:,} bp")
    print(f"• COVID-19 total codons: {metadata['covid_total_codons']:,}")
    print(f"• Influenza total codons: {metadata['influenza_total_codons']:,}")

    print(f"\n🦠 TOP AMINO ACIDS FROM VIRAL GENOME ANALYSIS:")
    print("-" * 50)
    print("COVID-19 Top Amino Acids:")
    for aa in results['covid_top_amino_acids'][:3]:
        if aa['amino_acid_code'] != '*':  # Skip stop codons
            print(f"  • {aa['amino_acid_name']} ({aa['amino_acid_code']}): {aa['percentage']:.1f}%")

    print("\nInfluenza Top Amino Acids:")
    for aa in results['influenza_top_amino_acids'][:3]:
        if aa['amino_acid_code'] != '*':  # Skip stop codons
            print(f"  • {aa['amino_acid_name']} ({aa['amino_acid_code']}): {aa['percentage']:.1f}%")

    target_amino_acids = results['target_amino_acids']

    print(f"\n📊 TARGET AMINO ACIDS FOR DIETARY RESTRICTION:")
    print("-" * 50)
    print("Based on viral genome analysis, these amino acids are highly")
    print("utilized by COVID-19 and Influenza viruses:")
    print()
    for aa in target_amino_acids:
        print(f"• {aa}")

    return results, target_amino_acids

def calculate_virus_specific_scores(food_data, virus_amino_acids):
    """Calculate food scores based on specific virus amino acid usage."""
    virus_scores = []

    # Create a mapping for virus-specific amino acids
    virus_aa_map = {
        'Leucine': 'Leucine_mg_per_100g',
        'Serine': 'Serine_mg_per_100g',
        'Threonine': 'Threonine_mg_per_100g',
        'Glutamine': 'Glutamine_mg_per_100g'
    }

    for i, food in enumerate(food_data['Food']):
        score = 0
        for aa_info in virus_amino_acids:
            aa_name = aa_info['amino_acid_name']
            if aa_name in virus_aa_map and aa_info['amino_acid_code'] != '*':
                aa_key = virus_aa_map[aa_name]
                score += food_data[aa_key][i] * (aa_info['percentage'] / 100)
        virus_scores.append(score)

    return virus_scores

def get_food_amino_acid_data():
    """Get amino acid content data for various foods."""
    # Amino acid content in foods (mg per 100g)
    food_amino_acid_data = {
        'Food': [
            # Low protein/amino acid foods
            'White Rice', 'Potato', 'Sweet Potato', 'Banana', 'Apple',
            'Cucumber', 'Lettuce', 'Watermelon', 'Orange', 'Carrot',
            'Celery', 'Spinach', 'Broccoli', 'Cauliflower', 'Bell Pepper',

            # Medium protein foods
            'Quinoa', 'Oats', 'Whole Wheat Bread', 'Pasta', 'Brown Rice',
            'Almonds', 'Walnuts', 'Sunflower Seeds', 'Avocado', 'Olive Oil',

            # High protein foods (to avoid/limit)
            'Chicken Breast', 'Beef', 'Salmon', 'Eggs', 'Greek Yogurt',
            'Cottage Cheese', 'Tofu', 'Lentils', 'Black Beans', 'Chickpeas'
        ],
        'Leucine_mg_per_100g': [
            # Low protein foods (very low leucine)
            150, 100, 120, 68, 19,
            60, 80, 18, 30, 72,
            40, 220, 190, 160, 80,

            # Medium protein foods
            810, 1200, 800, 400, 180,
            1100, 900, 750, 160, 5,

            # High protein foods
            1800, 1700, 1600, 1100, 950,
            1200, 800, 1850, 1400, 1500
        ],
        'Serine_mg_per_100g': [
            # Low protein foods
            180, 150, 140, 40, 17,
            80, 60, 9, 40, 80,
            50, 250, 200, 180, 90,

            # Medium protein foods
            460, 800, 450, 350, 200,
            600, 500, 400, 130, 3,

            # High protein foods
            900, 850, 800, 750, 600,
            700, 500, 1100, 900, 950
        ],
        'Threonine_mg_per_100g': [
            # Low protein foods
            140, 120, 110, 28, 12,
            70, 50, 7, 20, 60,
            35, 180, 160, 140, 70,

            # Medium protein foods
            380, 600, 350, 280, 150,
            450, 400, 350, 110, 2,

            # High protein foods
            1000, 950, 900, 600, 500,
            600, 400, 800, 700, 750
        ],
        'Glutamine_mg_per_100g': [
            # Low protein foods
            300, 250, 280, 150, 40,
            140, 120, 30, 100, 150,
            80, 600, 500, 400, 180,

            # Medium protein foods
            900, 1400, 800, 600, 350,
            1200, 1000, 900, 250, 8,

            # High protein foods
            3200, 3000, 2800, 1400, 1800,
            2000, 1600, 4000, 3500, 3800
        ],
        'Category': [
            # Low protein foods
            'Low Protein', 'Low Protein', 'Low Protein', 'Low Protein', 'Low Protein',
            'Low Protein', 'Low Protein', 'Low Protein', 'Low Protein', 'Low Protein',
            'Low Protein', 'Low Protein', 'Low Protein', 'Low Protein', 'Low Protein',

            # Medium protein foods
            'Medium Protein', 'Medium Protein', 'Medium Protein', 'Medium Protein', 'Medium Protein',
            'Medium Protein', 'Medium Protein', 'Medium Protein', 'Medium Protein', 'Medium Protein',

            # High protein foods
            'High Protein', 'High Protein', 'High Protein', 'High Protein', 'High Protein',
            'High Protein', 'High Protein', 'High Protein', 'High Protein', 'High Protein'
        ]
    }

    # Calculate total target amino acids for each food
    total_aa = []
    for i in range(len(food_amino_acid_data['Food'])):
        total = (food_amino_acid_data['Leucine_mg_per_100g'][i] +
                food_amino_acid_data['Serine_mg_per_100g'][i] +
                food_amino_acid_data['Threonine_mg_per_100g'][i] +
                food_amino_acid_data['Glutamine_mg_per_100g'][i])
        total_aa.append(total)

    food_amino_acid_data['Total_Target_AA'] = total_aa

    return food_amino_acid_data

def generate_virus_specific_recommendations(data, analysis_results):
    """Generate separate recommendations for each virus."""

    print(f"\n🦠 VIRUS-SPECIFIC FOOD RECOMMENDATIONS")
    print("=" * 70)

    # Calculate virus-specific scores
    covid_scores = calculate_virus_specific_scores(data, analysis_results['covid_top_amino_acids'][:3])
    flu_scores = calculate_virus_specific_scores(data, analysis_results['influenza_top_amino_acids'][:3])

    # Create sorted lists for each virus
    covid_food_data = list(zip(data['Food'], covid_scores, data['Category']))
    flu_food_data = list(zip(data['Food'], flu_scores, data['Category']))

    covid_food_data.sort(key=lambda x: x[1])  # Sort by score (ascending - lowest first)
    flu_food_data.sort(key=lambda x: x[1])

    # COVID-19 recommendations
    print(f"\n🧬 COVID-19 SPECIFIC RECOMMENDATIONS")
    print("-" * 50)
    print("Foods with lowest amino acids used by COVID-19:")
    print()
    print(f"{'Rank':<4} {'Food':<20} {'Score':<12} {'Category':<15}")
    print("-" * 55)

    for i, (food, score, category) in enumerate(covid_food_data[:15], 1):
        print(f"{i:<4} {food:<20} {score:<12.1f} {category:<15}")

    # Influenza recommendations
    print(f"\n🦠 INFLUENZA SPECIFIC RECOMMENDATIONS")
    print("-" * 50)
    print("Foods with lowest amino acids used by Influenza:")
    print()
    print(f"{'Rank':<4} {'Food':<20} {'Score':<12} {'Category':<15}")
    print("-" * 55)

    for i, (food, score, category) in enumerate(flu_food_data[:15], 1):
        print(f"{i:<4} {food:<20} {score:<12.1f} {category:<15}")

    return covid_food_data, flu_food_data

def generate_recommendations(data, target_amino_acids):
    """Generate dietary recommendations based on amino acid content."""

    print(f"\n🥗 RECOMMENDED FOODS (Low in Target Amino Acids)")
    print("=" * 70)

    # Create list of tuples and sort by total amino acids (ascending)
    food_data = list(zip(data['Food'], data['Total_Target_AA'], data['Category']))
    food_data.sort(key=lambda x: x[1])  # Sort by total AA content

    low_aa_foods = food_data[:15]  # Get 15 lowest

    print("These foods are naturally low in the amino acids most used by")
    print("COVID-19 and Influenza viruses:\n")

    print(f"{'Rank':<4} {'Food':<20} {'Total AA (mg/100g)':<18} {'Category':<15}")
    print("-" * 65)

    for i, (food, total_aa, category) in enumerate(low_aa_foods, 1):
        print(f"{i:<4} {food:<20} {total_aa:<18.0f} {category:<15}")

    print(f"\n🚫 FOODS TO LIMIT (High in Target Amino Acids)")
    print("=" * 70)

    # Sort by total amino acids (descending) and get top 10
    food_data.sort(key=lambda x: x[1], reverse=True)
    high_aa_foods = food_data[:10]

    print("These foods are high in amino acids frequently used by viruses")
    print("and might be limited during active infection:\n")

    print(f"{'Rank':<4} {'Food':<20} {'Total AA (mg/100g)':<18} {'Category':<15}")
    print("-" * 65)

    for i, (food, total_aa, category) in enumerate(high_aa_foods, 1):
        print(f"{i:<4} {food:<20} {total_aa:<18.0f} {category:<15}")

def create_virus_specific_charts(covid_food_data, flu_food_data):
    """Create separate charts for COVID-19 and Influenza food recommendations."""

    # COVID-19 Chart
    plt.figure(figsize=(14, 8))

    covid_low_foods = covid_food_data[:15]  # Top 15 lowest scoring foods
    foods, scores, categories = zip(*covid_low_foods)

    colors = ['green' if cat == 'Low Protein' else 'orange' if cat == 'Medium Protein' else 'red'
              for cat in categories]

    bars = plt.bar(range(len(foods)), scores, color=colors, alpha=0.7)
    plt.title('COVID-19 Specific: 15 Best Foods (Lowest Amino Acid Scores)',
              fontsize=16, fontweight='bold', pad=20)
    plt.ylabel('Weighted Amino Acid Score', fontsize=12)
    plt.xlabel('Foods (Ranked by COVID-19 Specific Score)', fontsize=12)
    plt.xticks(range(len(foods)), foods, rotation=45, ha='right')

    # Add value labels on bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + max(scores) * 0.01,
                f'{height:.1f}', ha='center', va='bottom', fontweight='bold', fontsize=10)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='green', alpha=0.7, label='Low Protein'),
                      Patch(facecolor='orange', alpha=0.7, label='Medium Protein'),
                      Patch(facecolor='red', alpha=0.7, label='High Protein')]
    plt.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig('covid19_food_recommendations.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Influenza Chart
    plt.figure(figsize=(14, 8))

    flu_low_foods = flu_food_data[:15]  # Top 15 lowest scoring foods
    foods, scores, categories = zip(*flu_low_foods)

    colors = ['green' if cat == 'Low Protein' else 'orange' if cat == 'Medium Protein' else 'red'
              for cat in categories]

    bars = plt.bar(range(len(foods)), scores, color=colors, alpha=0.7)
    plt.title('Influenza Specific: 15 Best Foods (Lowest Amino Acid Scores)',
              fontsize=16, fontweight='bold', pad=20)
    plt.ylabel('Weighted Amino Acid Score', fontsize=12)
    plt.xlabel('Foods (Ranked by Influenza Specific Score)', fontsize=12)
    plt.xticks(range(len(foods)), foods, rotation=45, ha='right')

    # Add value labels on bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + max(scores) * 0.01,
                f'{height:.1f}', ha='center', va='bottom', fontweight='bold', fontsize=10)

    # Add legend
    plt.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig('influenza_food_recommendations.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_visualization(data):
    """Create a visualization of amino acid content in foods."""

    # Create a bar chart showing total target amino acids by food category
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Chart 1: Average amino acid content by category
    categories = ['Low Protein', 'Medium Protein', 'High Protein']
    category_totals = {cat: [] for cat in categories}

    for i, category in enumerate(data['Category']):
        category_totals[category].append(data['Total_Target_AA'][i])

    category_averages = [(cat, sum(values)/len(values)) for cat, values in category_totals.items()]
    category_averages.sort(key=lambda x: x[1])

    cat_names, cat_values = zip(*category_averages)
    colors_map = {'Low Protein': 'green', 'Medium Protein': 'orange', 'High Protein': 'red'}
    colors1 = [colors_map[cat] for cat in cat_names]

    bars1 = ax1.bar(cat_names, cat_values, color=colors1, alpha=0.7)
    ax1.set_title('Average Target Amino Acids by Food Category', fontweight='bold')
    ax1.set_ylabel('Total Target Amino Acids (mg/100g)')
    ax1.tick_params(axis='x', rotation=45)

    # Add value labels on bars
    for i, bar in enumerate(bars1):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 50,
                f'{height:.0f}', ha='center', va='bottom', fontweight='bold')

    # Chart 2: Top 15 lowest amino acid foods
    food_data = list(zip(data['Food'], data['Total_Target_AA'], data['Category']))
    food_data.sort(key=lambda x: x[1])  # Sort by total AA
    low_foods = food_data[:15]

    foods, aa_values, categories = zip(*low_foods)
    colors2 = ['green' if cat == 'Low Protein' else 'orange' if cat == 'Medium Protein' else 'red'
              for cat in categories]

    bars2 = ax2.bar(range(len(foods)), aa_values, color=colors2, alpha=0.7)
    ax2.set_title('15 Foods Lowest in Target Amino Acids', fontweight='bold')
    ax2.set_ylabel('Total Target Amino Acids (mg/100g)')
    ax2.set_xlabel('Foods (Ranked by Total AA Content)')
    ax2.set_xticks(range(len(foods)))
    ax2.set_xticklabels(foods, rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig('amino_acid_food_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def immune_supporting_recommendations():
    """Provide additional immune-supporting dietary recommendations."""

    print(f"\n🛡️ IMMUNE SYSTEM SUPPORT RECOMMENDATIONS")
    print("=" * 70)

    print("While limiting certain amino acids, it's crucial to support")
    print("immune function with these nutrients and foods:\n")

    immune_foods = {
        "Vitamin C Rich": ["Citrus fruits", "Bell peppers", "Strawberries", "Kiwi", "Broccoli"],
        "Vitamin D Sources": ["Sunlight exposure", "Fatty fish (small amounts)", "Egg yolks", "Fortified foods"],
        "Zinc Rich (moderate)": ["Pumpkin seeds", "Cashews", "Spinach", "Dark chocolate"],
        "Antioxidants": ["Berries", "Green tea", "Turmeric", "Ginger", "Garlic"],
        "Prebiotics": ["Garlic", "Onions", "Bananas", "Asparagus", "Jerusalem artichokes"]
    }

    for category, foods in immune_foods.items():
        print(f"🔹 {category}:")
        for food in foods:
            print(f"   • {food}")
        print()

    print("⚠️  IMPORTANT DISCLAIMERS:")
    print("-" * 30)
    print("• This analysis is theoretical and based on genomic data")
    print("• Consult healthcare providers before making dietary changes")
    print("• Maintain balanced nutrition for overall health")
    print("• This is NOT medical advice for treating viral infections")
    print("• Severe protein restriction can harm immune function")

def main():
    """Main function to run the nutritional analysis."""

    # Analyze viral amino acids
    analysis_results, target_amino_acids = analyze_viral_amino_acids()

    if not analysis_results:
        print("❌ Cannot proceed without codon analysis results.")
        return 1

    # Get food amino acid data
    food_data = get_food_amino_acid_data()

    # Generate virus-specific recommendations
    covid_food_data, flu_food_data = generate_virus_specific_recommendations(food_data, analysis_results)

    # Generate general recommendations
    generate_recommendations(food_data, target_amino_acids)

    # Create virus-specific charts
    print(f"\n📊 GENERATING VIRUS-SPECIFIC CHARTS...")
    print("-" * 50)
    create_virus_specific_charts(covid_food_data, flu_food_data)

    # Create general visualization
    create_visualization(food_data)

    # Provide immune support recommendations
    immune_supporting_recommendations()

    print(f"\n📈 Analysis complete! Charts saved:")
    print("  • covid19_food_recommendations.png")
    print("  • influenza_food_recommendations.png")
    print("  • amino_acid_food_analysis.png")
    print("=" * 70)

    return 0

if __name__ == "__main__":
    main()