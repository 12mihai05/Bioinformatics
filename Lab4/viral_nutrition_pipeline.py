#!/usr/bin/env python3
"""
Viral Genome Analysis & Nutritional Recommendations Pipeline

This script orchestrates the complete analysis workflow:
1. Downloads viral genome FASTA files (if needed)
2. Runs codon frequency analysis
3. Generates nutritional recommendations based on results
4. Creates comprehensive reports and visualizations
"""

import subprocess
import sys
import os
from pathlib import Path
import json

def check_dependencies():
    """Check if required files and dependencies are available."""
    print("🔍 Checking dependencies and files...")

    required_files = [
        "covid19_genome.fasta",
        "influenza_genome.fasta",
        "codon_analysis.py",
        "nutritional_analysis.py"
    ]

    missing_files = []
    for file in required_files:
        if not Path(file).exists():
            missing_files.append(file)

    if missing_files:
        print(f"❌ Missing files: {', '.join(missing_files)}")
        if any('fasta' in f for f in missing_files):
            print("📥 Downloading missing genome files...")
            download_genomes()
        return False

    print("✅ All required files found!")
    return True

def download_genomes():
    """Download genome FASTA files if missing."""
    genome_urls = {
        "covid19_genome.fasta": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_045512.2&rettype=fasta&retmode=text",
        "influenza_genome.fasta": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=CY121680.1&rettype=fasta&retmode=text"
    }

    for filename, url in genome_urls.items():
        if not Path(filename).exists():
            print(f"⬇️  Downloading {filename}...")
            try:
                result = subprocess.run(
                    ["curl", "-o", filename, url],
                    capture_output=True,
                    text=True,
                    check=True
                )
                print(f"✅ {filename} downloaded successfully")
            except subprocess.CalledProcessError as e:
                print(f"❌ Failed to download {filename}: {e}")
                return False
    return True

def run_codon_analysis():
    """Run the codon frequency analysis."""
    print("\n🧬 STEP 1: Running Codon Frequency Analysis")
    print("=" * 50)

    try:
        result = subprocess.run(
            [sys.executable, "codon_analysis.py"],
            capture_output=True,
            text=True,
            check=True
        )
        print(result.stdout)
        if result.stderr:
            print(f"⚠️  Warnings: {result.stderr}")

        # Check if results file was created
        if Path("codon_analysis_results.json").exists():
            print("✅ Codon analysis completed successfully!")
            return True
        else:
            print("❌ Results file not created")
            return False

    except subprocess.CalledProcessError as e:
        print(f"❌ Codon analysis failed: {e}")
        print(f"Error output: {e.stderr}")
        return False

def run_nutritional_analysis():
    """Run the nutritional recommendations analysis."""
    print("\n🥗 STEP 2: Generating Nutritional Recommendations")
    print("=" * 50)

    try:
        result = subprocess.run(
            [sys.executable, "nutritional_analysis.py"],
            capture_output=True,
            text=True,
            check=True
        )
        print(result.stdout)
        if result.stderr:
            print(f"⚠️  Warnings: {result.stderr}")

        print("✅ Nutritional analysis completed successfully!")
        return True

    except subprocess.CalledProcessError as e:
        print(f"❌ Nutritional analysis failed: {e}")
        print(f"Error output: {e.stderr}")
        return False

def display_summary():
    """Display a summary of generated files and key results."""
    print("\n📊 PIPELINE SUMMARY")
    print("=" * 50)

    # Check for generated files
    generated_files = []
    expected_files = [
        "codon_analysis_results.json",
        "covid19_top10_codons.png",
        "influenza_top10_codons.png",
        "amino_acid_food_analysis.png"
    ]

    for file in expected_files:
        if Path(file).exists():
            generated_files.append(file)

    print(f"📁 Generated Files ({len(generated_files)}/{len(expected_files)}):")
    for file in generated_files:
        size = Path(file).stat().st_size
        print(f"   ✅ {file} ({size:,} bytes)")

    # Display key results from JSON file
    if Path("codon_analysis_results.json").exists():
        try:
            with open("codon_analysis_results.json", "r") as f:
                results = json.load(f)

            print(f"\n🎯 Key Findings:")
            print(f"   • COVID-19 genome: {results['analysis_metadata']['covid_genome_length']:,} bp")
            print(f"   • Influenza genome: {results['analysis_metadata']['influenza_genome_length']:,} bp")
            print(f"   • Target amino acids for dietary focus:")

            for aa in results['target_amino_acids']:
                print(f"     - {aa}")

        except (json.JSONDecodeError, KeyError) as e:
            print(f"⚠️  Could not read results file: {e}")

    print(f"\n🎉 Pipeline completed successfully!")
    print(f"📖 Review the generated charts and nutritional recommendations.")

def main():
    """Main pipeline orchestrator."""
    print("🦠 VIRAL GENOME ANALYSIS & NUTRITIONAL PIPELINE")
    print("=" * 60)
    print("This pipeline analyzes viral genomes and provides nutritional")
    print("recommendations based on amino acid frequencies.")
    print("=" * 60)

    # Step 0: Check dependencies
    if not check_dependencies():
        print("❌ Dependency check failed. Please ensure all required files are present.")
        return 1

    # Step 1: Run codon analysis
    if not run_codon_analysis():
        print("❌ Pipeline failed at codon analysis step.")
        return 1

    # Step 2: Run nutritional analysis
    if not run_nutritional_analysis():
        print("❌ Pipeline failed at nutritional analysis step.")
        return 1

    # Step 3: Display summary
    display_summary()

    return 0

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)