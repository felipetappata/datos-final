#!/usr/bin/env python3
"""
Table 3 Generator for Al Sadoon et al. (2019) Replication
=========================================================

This script reads Monte Carlo simulation results from CSV files and generates
a LaTeX tabular environment (not full table) replicating Table 3 from Al Sadoon et al. (2019).
This is exactly like Table 2 but for N = 5000.

Author: Felipe I. Tappata
Date: January 2025
"""

import pandas as pd
import os
import glob
import re
from pathlib import Path

# Configuration - Modified for Table 3 (N = 5000)
PARTIAL_OUTPUT_DIR = "output/partial_tab3"  # Updated to use partial_tab3 for N = 5000
TABLE_OUTPUT_DIR = "output/tables"
TABLE_OUTPUT_FILE = "table3.tex"  # Changed from table2.tex to table3.tex

def find_first_significant_decimal(value):
    """
    Find the position of the first significant (non-zero) digit after the decimal point.
    
    Args:
        value: Numeric value
    
    Returns:
        Position of first significant digit (1-indexed), or None if no significant decimals
    """
    if pd.isna(value) or value == 0:
        return None
    
    # Get the decimal part
    decimal_part = abs(value) % 1
    if decimal_part == 0:
        return None
    
    # Convert to string and find first non-zero digit after decimal
    decimal_str = f"{decimal_part:.10f}"  # Use more precision to be safe
    decimal_str = decimal_str[2:]  # Remove "0."
    
    for i, digit in enumerate(decimal_str):
        if digit != '0':
            return i + 1  # 1-indexed position
    
    return None

def format_number_enhanced(value, underline_significant=True, use_decimal_alignment=True):
    """
    Enhanced number formatting with optional underlining of first significant decimal
    and support for decimal alignment with proper centering.
    
    Args:
        value: Numeric value to format
        underline_significant: Whether to underline the first significant decimal
        use_decimal_alignment: Whether to format for decimal alignment (plain text vs math mode)
    
    Returns:
        Formatted string for LaTeX
    """
    if pd.isna(value):
        return "---"
    
    # Handle zero
    if value == 0:
        if use_decimal_alignment:
            return "$\\phantom{-}0.00000$"
        else:
            return "$0.00000$"
    
    # Check if we need scientific notation
    if abs(value) < 1e-5 and abs(value) > 0:
        # Use scientific notation
        formatted = f"{value:.2e}"
        if "e" in formatted:
            base, exp = formatted.split("e")
            exp = int(exp)
            if base.startswith('-'):
                base = f"-{base[1:]}"
            if use_decimal_alignment:
                return f"${base} \\times 10^{{{exp}}}$"  # Fixed: restored \times in math mode
            else:
                return f"${base} \\times 10^{{{exp}}}$"
    
    # Regular formatting with 5 decimal places
    formatted = f"{value:.5f}"
    
    if underline_significant:
        # Find the first significant decimal
        sig_pos = find_first_significant_decimal(value)
        if sig_pos is not None:
            # Split the formatted number
            if '.' in formatted:
                integer_part, decimal_part = formatted.split('.')
                if sig_pos <= len(decimal_part):
                    # Underline the significant digit
                    sig_digit = decimal_part[sig_pos-1]
                    new_decimal = (decimal_part[:sig_pos-1] + 
                                 f"\\underline{{{sig_digit}}}" + 
                                 decimal_part[sig_pos:])
                    formatted = f"{integer_part}.{new_decimal}"
    
    if use_decimal_alignment:
        # For decimal alignment with proper centering, pad positive numbers with phantom minus
        # Always wrap in math mode for proper formatting
        if not formatted.startswith('-'):
            # Add phantom minus sign for alignment
            formatted = f"$\\phantom{{-}}{formatted}$"
        else:
            formatted = f"${formatted}$"
        return formatted
    else:
        # Handle negative signs properly for LaTeX math mode  
        if formatted.startswith('-'):
            # For negative numbers, remove leading zero after minus sign
            if formatted.startswith('-0.'):
                formatted = f"-{formatted[2:]}"  # Remove the '0' after '-'
            return f"${formatted}$"
        else:
            # For positive numbers less than 1, remove leading zero
            if formatted.startswith('0.'):
                formatted = formatted[1:]  # Remove leading '0'
            return f"${formatted}$"

def format_number(value):
    """
    Original format function for backward compatibility.
    Shows 5 decimal places. Uses format similar to stata output:
    - Negative: -.00351 (no leading zero)
    - Positive < 1: .13394 (no leading zero) 
    - Positive >= 1: 1.13394 (with leading digit)
    
    Args:
        value: Numeric value to format
    
    Returns:
        Formatted string for LaTeX
    """
    return format_number_enhanced(value, underline_significant=False, use_decimal_alignment=False)

def read_simulation_results():
    """
    Read all experiment CSV files and organize them into a structured format.
    Filters for N = 5000 results only.
    
    Returns:
        Dictionary with results organized by (experiment, model, rho)
    """
    results = {}
    
    # Pattern to match our CSV files
    pattern = os.path.join(PARTIAL_OUTPUT_DIR, "*.csv")
    csv_files = glob.glob(pattern)
    
    print(f"Found {len(csv_files)} CSV files in {PARTIAL_OUTPUT_DIR}")
    
    for file_path in csv_files:
        filename = os.path.basename(file_path)
        print(f"Processing: {filename}")
        
        # Parse filename to extract parameters
        # Expected format: exp{1-5}_model{A|B}_N{N}_rho{rho}.csv or endo_model{A}_N{N}_rho{rho}.csv
        try:
            parts = filename.replace('.csv', '').split('_')
            
            # Handle both naming conventions
            if filename.startswith('exp'):
                experiment = int(parts[0].replace('exp', ''))  # 1-5
                model = parts[1].replace('model', '')  # A or B
                n_value = int(parts[2].replace('N', ''))
                rho_value = float(parts[3].replace('rho', ''))
            elif filename.startswith('endo_model'):
                # This is the basic endogenous experiment (experiment 1 equivalent)
                experiment = 1
                model = parts[1].replace('model', '')  # A or B
                n_value = int(parts[2].replace('N', ''))
                rho_value = float(parts[3].replace('rho', ''))
            else:
                print(f"  -> Skipping unknown format: {filename}")
                continue
            
            # Only process N = 5000 files for Table 3
            if n_value != 5000:
                print(f"  -> Skipping N={n_value} (Table 3 is for N=5000)")
                continue
            
            # Read the CSV file
            df = pd.read_csv(file_path)
            
            # Extract the results we need
            if len(df) > 0:
                # Get the first (and likely only) row of results
                row = df.iloc[0]
                
                # Store results - key is (experiment, model, rho)
                key = (experiment, model, rho_value)
                results[key] = {
                    'ab_bias': row.get('AB_bias', pd.NA),
                    'ab_se': row.get('AB_se', pd.NA),
                    'sys_bias': row.get('SYS_bias', pd.NA),
                    'sys_se': row.get('SYS_se', pd.NA),
                    'ab_valid': row.get('AB_valid', pd.NA),
                    'sys_valid': row.get('SYS_valid', pd.NA),
                    'n_value': n_value
                }
                
                print(f"  -> Experiment {experiment}, Model {model}, ρ={rho_value}, N={n_value}")
                if not pd.isna(row.get('AB_bias')):
                    print(f"     AB bias: {row.get('AB_bias'):.6f}, SYS bias: {row.get('SYS_bias'):.6f}")
            
        except Exception as e:
            print(f"  -> Error parsing {filename}: {e}")
            continue
    
    return results

def generate_latex_table(results, output_path="output/tables/table3.tex"):
    """
    Generate the LaTeX tabular environment replicating Table 3 from Al Sadoon et al. (2019).
    Uses decimal alignment and underlining of first significant decimal.
    Only outputs the tabular environment, not the full table wrapper.
    
    Args:
        results: Dictionary with simulation results
        output_path: Path where to save the LaTeX table
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Experiment names for headers
    experiment_names = {
        1: "Very short $T$ ($T = 4$)",
        2: "More sample selection (25\\%)",
        3: "Increasing the ratio of variances: $\\sigma_\eta/\\sigma_\\varepsilon = 2$",
        4: "Reducing the correlation of the errors: $\\rho = 0.25$",
        5: "Non-stationary time-varying error components"
    }
    
    with open(output_path, 'w') as f:
        # Table with 8 columns: Model (l), Type (l), and 6 centered data columns (c)
        # Using centered columns since we handle alignment manually with \phantom
        f.write("\\begin{tabular}{@{}ll*{6}{c}@{}}\n")
        f.write("\\toprule\n")
        
        # Column headers - first two columns are Model and blank for bias/s.e.
        f.write("& & \\multicolumn{2}{c}{$\\rho = 0.25$} & \\multicolumn{2}{c}{$\\rho = 0.5$} & \\multicolumn{2}{c}{$\\rho = 0.75$} \\\\\n")
        f.write("\\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}\n")
        f.write("Model & & AB & SYS & AB & SYS & AB & SYS \\\\\n")
        f.write("\\midrule\n")
        
        # Generate rows for each experiment
        for exp_num in [1, 2, 3, 4, 5]:
            # Experiment header spanning all 8 columns, centered
            roman_num = ['I', 'II', 'III', 'IV', 'V'][exp_num-1]
            f.write(f"\\multicolumn{{8}}{{c}}{{Experiment {roman_num}: {experiment_names[exp_num]}}} \\\\\n")
            
            # Model A and B rows for this experiment
            for model in ['A', 'B']:
                # Bias row
                f.write(f"{model} & bias")
                
                # Get bias results for each rho value
                for rho in [0.25, 0.50, 0.75]:
                    key = (exp_num, model, rho)
                    
                    if key in results:
                        result = results[key]
                        ab_bias = format_number_enhanced(result['ab_bias'], underline_significant=True, use_decimal_alignment=True)
                        sys_bias = format_number_enhanced(result['sys_bias'], underline_significant=True, use_decimal_alignment=True)
                    else:
                        ab_bias = "---"
                        sys_bias = "---"
                    
                    f.write(f" & {ab_bias} & {sys_bias}")
                
                f.write(" \\\\\n")
                
                # Standard error row for this model
                f.write(" & s.e.")
                for rho in [0.25, 0.50, 0.75]:
                    key = (exp_num, model, rho)
                    
                    if key in results:
                        result = results[key]
                        ab_se = format_number_enhanced(result['ab_se'], underline_significant=True, use_decimal_alignment=True)
                        sys_se = format_number_enhanced(result['sys_se'], underline_significant=True, use_decimal_alignment=True)
                    else:
                        ab_se = "---"
                        sys_se = "---"
                    
                    f.write(f" & {ab_se} & {sys_se}")
                
                f.write(" \\\\\n")
            
            # Add space between experiments (except after the last one)
            if exp_num < 5:
                f.write("\\addlinespace\n")
        
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")

def main():
    """
    Main function to coordinate the table generation process.
    """
    print("=" * 70)
    print("Table 3 Generator for Al Sadoon et al. (2019) Replication")
    print("=" * 70)
    print()
    
    # Read simulation results
    print("Reading simulation results for N = 5000...")
    results = read_simulation_results()
    
    print(f"\nLoaded {len(results)} parameter combinations")
    print()
    
    # Generate LaTeX table
    output_path = os.path.join(TABLE_OUTPUT_DIR, TABLE_OUTPUT_FILE)
    print(f"Generating LaTeX table: {output_path}")
    generate_latex_table(results, output_path)
    
    # Check for missing combinations
    print("\nChecking for missing combinations...")
    print("-" * 40)
    
    experiments = [1, 2, 3, 4, 5]
    models = ["A", "B"]
    rho_values = [0.25, 0.50, 0.75]
    
    total_expected = len(experiments) * len(models) * len(rho_values)
    available = 0
    
    for exp in experiments:
        for model in models:
            for rho in rho_values:
                key = (exp, model, rho)
                if key in results and not pd.isna(results[key].get('ab_bias', pd.NA)):
                    available += 1
                else:
                    print(f"Missing: Experiment {exp} model {model}, rho={rho}")
    
    print(f"Available: {available}/{total_expected} parameter combinations")
    
    if available < total_expected:
        print("\nNote: Missing data will appear as em-dashes (—) in the table")
    
    print()
    print("=" * 70)
    print("Table 3 generation complete!")
    print("=" * 70)
    print()
    print("To include in your LaTeX document:")
    print(f"\\input{{{output_path}}}")
    print()
    print("Required LaTeX packages:")
    print("\\usepackage{booktabs}  % for \\toprule, \\midrule, \\bottomrule")
    print("\\usepackage{array}     % for enhanced tabular environments")
    print("\\usepackage{ulem}      % for \\underline (if not already loaded)")
    print()
    print("Note: Numbers are formatted with manual decimal alignment using \\phantom{-}")
    print("      and first significant decimal digits are underlined for better readability.")

if __name__ == "__main__":
    main()
