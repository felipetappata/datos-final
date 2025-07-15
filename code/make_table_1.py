#!/usr/bin/env python3
"""
Table 1 Generator for Al Sadoon et al. (2019) Replication
=========================================================

Author: Felipe I. Tappata
Date: July 2025
"""

import pandas as pd
import os
import glob
import re
from pathlib import Path

# Configuration
PARTIAL_OUTPUT_DIR = "output/partial"
TABLE_OUTPUT_DIR = "output/tables"
TABLE_OUTPUT_FILE = "table1.tex"

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

def format_number_enhanced(value, underline_significant=True, use_centered_alignment=True):
    """
    Enhanced number formatting with optional underlining of first significant decimal
    and centered alignment that preserves decimal alignment.
    
    Args:
        value: Numeric value to format
        underline_significant: Whether to underline the first significant decimal
        use_centered_alignment: Whether to use centered alignment (vs math mode)
    
    Returns:
        Formatted string for LaTeX
    """
    if pd.isna(value):
        return "---"
    
    # Handle zero
    if value == 0:
        if use_centered_alignment:
            return "$\\phantom{-}0.00000$"  # Add phantom minus for alignment in math mode
        else:
            return "$0.00000$"
    
    # Check if we need scientific notation
    if abs(value) < 1e-5 and abs(value) > 0:
        # Use scientific notation
        formatted = f"{value:.2e}"
        if "e" in formatted:
            base, exp = formatted.split("e")
            exp = int(exp)
            if use_centered_alignment:
                if base.startswith('-'):
                    return f"${base} \\times 10^{{{exp}}}$"
                else:
                    return f"$\\phantom{{-}}{base} \\times 10^{{{exp}}}$"
            else:
                if base.startswith('-'):
                    return f"${base} \\times 10^{{{exp}}}$"
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
    
    if use_centered_alignment:
        # For centered alignment, add phantom minus to positive numbers for consistent spacing
        # Always wrap in math mode for proper formatting
        if formatted.startswith('-'):
            return f"${formatted}$"
        else:
            return f"$\\phantom{{-}}{formatted}$"
    else:
        # Handle negative signs properly for LaTeX math mode
        if formatted.startswith('-'):
            return f"$-{formatted[1:]}$"
        else:
            return f"${formatted}$"

def format_number(value):
    """
    Original format_number function for backward compatibility.
    Uses the enhanced version with default settings.
    """
    return format_number_enhanced(value, underline_significant=False, use_centered_alignment=False)

def read_simulation_results():
    """
    Read all simulation result CSV files and organize them into a structured format.
    
    Returns:
        Dictionary with results organized by scenario
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
        # Expected format: {endo|nonendo}_model{A|B}_N{500|5000}_rho{0.25|0.50|0.75}.csv
        try:
            parts = filename.replace('.csv', '').split('_')
            selection_type = 'endogenous' if parts[0] == 'endo' else 'nonendogenous'
            model = parts[1].replace('model', '')  # A or B
            n_value = int(parts[2].replace('N', ''))  # 500 or 5000
            rho_value = float(parts[3].replace('rho', ''))  # 0.25, 0.50, 0.75
            
            # Read the CSV file
            df = pd.read_csv(file_path)
            
            # Extract the results we need
            if len(df) > 0:
                # Get the first (and likely only) row of results
                row = df.iloc[0]
                
                # Store results
                key = (selection_type, model, n_value, rho_value)
                results[key] = {
                    'ab_bias': row.get('AB_bias', pd.NA),
                    'ab_se': row.get('AB_se', pd.NA),
                    'sys_bias': row.get('SYS_bias', pd.NA),
                    'sys_se': row.get('SYS_se', pd.NA),
                    'ab_valid': row.get('AB_valid', pd.NA),
                    'sys_valid': row.get('SYS_valid', pd.NA)
                }
                
                print(f"  -> {selection_type}, Model {model}, N={n_value}, ρ={rho_value}")
                if not pd.isna(row.get('AB_bias')):
                    print(f"     AB bias: {row.get('AB_bias'):.6f}, SYS bias: {row.get('SYS_bias'):.6f}")
            
        except Exception as e:
            print(f"  -> Error parsing {filename}: {e}")
            continue
    
    return results

def check_missing_files():
    """
    Check for missing parameter combinations and warn about them.
    
    Returns:
        Set of missing parameter combinations
    """
    models = ['A', 'B']
    n_values = [500, 5000]
    rho_values = [0.25, 0.50, 0.75]
    selection_types = ['nonendogenous', 'endogenous']
    
    missing = set()
    
    for selection in selection_types:
        selection_prefix = 'endo' if selection == 'endogenous' else 'nonendo'
        for model in models:
            for n in n_values:
                for rho in rho_values:
                    # Try different rho formats
                    rho_formats = [f"{rho:.2f}", f"{rho:.1f}"]
                    found = False
                    
                    for rho_str in rho_formats:
                        filename = f"{selection_prefix}_model{model}_N{n}_rho{rho_str}.csv"
                        filepath = os.path.join(PARTIAL_OUTPUT_DIR, filename)
                        if os.path.exists(filepath):
                            found = True
                            break
                    
                    if not found:
                        missing.add((selection, model, n, rho))
                        print(f"Warning: {filename} not found; treating as NA")
    
    return missing

def generate_latex_table(results, output_path="output/tables/table1.tex"):
    """
    Generate the LaTeX table replicating Table 1 from Al Sadoon et al. (2019).
    Uses decimal alignment and underlining of first significant decimal.
    
    Args:
        results: Dictionary with simulation results
        output_path: Path where to save the LaTeX table
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    latex_content = []
    
    # Table header with centered columns for numeric data
    # Use regular centered columns (c) instead of decimal alignment (D)
    latex_content.append("\\begin{tabular}{@{}ccc*{4}{c}@{}}")
    latex_content.append("\\toprule")
    latex_content.append("& & & \\multicolumn{2}{c}{No endogenous} & \\multicolumn{2}{c}{Endogenous} \\\\")
    latex_content.append("& & & \\multicolumn{2}{c}{selection} & \\multicolumn{2}{c}{selection} \\\\")
    latex_content.append("\\cmidrule(lr){4-5} \\cmidrule(lr){6-7}")
    latex_content.append("Select. & & & (1) & (2) & (3) & (4) \\\\")
    latex_content.append("Model & $\\rho$ & & AB & SYS & AB & SYS \\\\")
    latex_content.append("\\midrule")
    
    # Data rows - reorganized order: A models (N=500, N=5000), then B models (N=500, N=5000)
    models = ["A", "B"]
    sample_sizes = [500, 5000]
    rho_values = [0.25, 0.50, 0.75]
    
    def format_enhanced_number(value):
        """Format number with centered alignment and significant digit underlining"""
        return format_number_enhanced(value, underline_significant=True, use_centered_alignment=True)
    
    for model in models:
        for N in sample_sizes:
            # Sample size header with rules above and below
            latex_content.append("\\midrule")
            latex_content.append(f"\\multicolumn{{7}}{{c}}{{$N = {N}$}} \\\\")
            latex_content.append("\\midrule")
            
            for rho in rho_values:
                # Get results for non-endogenous and endogenous
                nonendo_key = ('nonendogenous', model, N, rho)
                endo_key = ('endogenous', model, N, rho)
                
                nonendo_results = results.get(nonendo_key, {})
                endo_results = results.get(endo_key, {})
                
                # Format the values with enhanced formatting
                ab_bias_nonendo = format_enhanced_number(nonendo_results.get('ab_bias', pd.NA))
                ab_se_nonendo = format_enhanced_number(nonendo_results.get('ab_se', pd.NA))
                sys_bias_nonendo = format_enhanced_number(nonendo_results.get('sys_bias', pd.NA))
                sys_se_nonendo = format_enhanced_number(nonendo_results.get('sys_se', pd.NA))
                
                ab_bias_endo = format_enhanced_number(endo_results.get('ab_bias', pd.NA))
                ab_se_endo = format_enhanced_number(endo_results.get('ab_se', pd.NA))
                sys_bias_endo = format_enhanced_number(endo_results.get('sys_bias', pd.NA))
                sys_se_endo = format_enhanced_number(endo_results.get('sys_se', pd.NA))
                
                # Format rho value in math mode
                rho_str = f"$.{int(rho*100):02d}$" if rho != 0.75 else "$.75$"
                
                # Bias row - model appears on this row
                latex_content.append(f"{model} & {rho_str} & bias & {ab_bias_nonendo} & {sys_bias_nonendo} & {ab_bias_endo} & {sys_bias_endo} \\\\")
                
                # Standard error row - no model letter
                latex_content.append(f"& & s.e. & {ab_se_nonendo} & {sys_se_nonendo} & {ab_se_endo} & {sys_se_endo} \\\\")
    
    # Table footer
    latex_content.append("\\bottomrule")
    latex_content.append("\\end{tabular}")
    
    # Write to file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(latex_content))
    
    return output_path

def main():
    """Main function to generate Table 1."""
    
    print("=" * 70)
    print("TABLE 1 GENERATOR - Al Sadoon et al. (2019) Replication")
    print("=" * 70)
    
    # Configuration
    partial_dir = "output/partial"
    output_path = "output/tables/table1.tex"
    
    print(f"Reading results from: {partial_dir}")
    print(f"Output will be saved to: {output_path}")
    print()
    
    # Load results data
    print("Loading simulation results...")
    
    # Check for missing files first
    missing = check_missing_files()
    if missing:
        print(f"\nFound {len(missing)} missing parameter combinations:")
        for selection, model, n, rho in missing:
            print(f"  - {selection} Model {model}, N={n}, ρ={rho}")
        print()
    
    results = read_simulation_results()
    
    # Generate LaTeX table
    print("Generating LaTeX table...")
    output_file = generate_latex_table(results, output_path)
    
    print(f"✓ LaTeX table successfully generated: {output_file}")
    print()
    
    # Summary of data availability
    print("Data availability summary:")
    print("-" * 40)
    
    models = ["A", "B"]
    sample_sizes = [500, 5000]
    rho_values = [0.25, 0.50, 0.75]
    selection_types = ["nonendogenous", "endogenous"]
    
    total_expected = len(models) * len(sample_sizes) * len(rho_values) * len(selection_types)
    available = 0
    
    for selection in selection_types:
        for model in models:
            for N in sample_sizes:
                for rho in rho_values:
                    key = (selection, model, N, rho)
                    if key in results and not pd.isna(results[key].get('ab_bias', pd.NA)):
                        available += 1
                    else:
                        print(f"Missing: {selection} model {model}, N={N}, rho={rho}")
    
    print(f"Available: {available}/{total_expected} parameter combinations")
    
    if available < total_expected:
        print("\nNote: Missing data will appear as em-dashes (—) in the table")
    
    print()
    print("=" * 70)
    print("Table generation complete!")
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
    print("Note: Numbers are formatted with phantom spacing for alignment and first significant")
    print("      decimal digits are underlined for better readability.")

if __name__ == "__main__":
    main()
