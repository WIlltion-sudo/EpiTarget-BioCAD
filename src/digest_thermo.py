import math
import re
import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# 1. THERMODYNAMIC PARAMETERS (SantaLucia 1998)
# Enthalpy (dH, kcal/mol) and Entropy (dS, cal/K*mol) for Nearest-Neighbor pairs
# ==============================================================================
NN_PARAMS = {
    'AA': {'dH': -7.9, 'dS': -22.2}, 'TT': {'dH': -7.9, 'dS': -22.2},
    'AT': {'dH': -7.2, 'dS': -20.4}, 'TA': {'dH': -7.2, 'dS': -21.3},
    'CA': {'dH': -8.5, 'dS': -22.7}, 'TG': {'dH': -8.5, 'dS': -22.7},
    'GT': {'dH': -8.4, 'dS': -22.4}, 'AC': {'dH': -8.4, 'dS': -22.4},
    'CT': {'dH': -7.8, 'dS': -21.0}, 'AG': {'dH': -7.8, 'dS': -21.0},
    'GA': {'dH': -8.2, 'dS': -22.2}, 'TC': {'dH': -8.2, 'dS': -22.2},
    'CG': {'dH': -10.6, 'dS': -27.2}, 'GC': {'dH': -9.8, 'dS': -24.4},
    'GG': {'dH': -8.0, 'dS': -19.9}, 'CC': {'dH': -8.0, 'dS': -19.9}
}

# Terminal initiation penalties
INIT_GC = {'dH': 0.1, 'dS': -2.8}
INIT_AT = {'dH': 2.3, 'dS': 4.1}

# Standard Restriction Enzyme Database (Name: (Recognition Site, Cut Offset from 5'))
RESTRICTION_ENZYMES = {
    'EcoRI':   ('GAATTC', 1),   # G^AATTC
    'BamHI':   ('GGATCC', 1),   # G^GATCC
    'HindIII': ('AAGCTT', 1),   # A^AGCTT
    'NotI':    ('GCGGCCGC', 2), # GC^GGCCGC
    'TaqI':    ('TCGA', 1),     # T^CGA
    'HaeIII':  ('GGCC', 2)      # GG^CC
}

# Standard Molecular Weight DNA Ladder (in Base Pairs)
STANDARD_LADDER = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 750, 500, 250, 100, 50, 10]


# ==============================================================================
# 2. THERMODYNAMIC NEAREST-NEIGHBOR ENGINE
# ==============================================================================
def calculate_tm_nearest_neighbor(sequence, dnac_nM=500.0, na_mM=50.0):
    """
    Calculates Tm (°C) using SantaLucia (1998) nearest-neighbor thermodynamics.
    """
    seq = sequence.upper().strip()
    if not set(seq).issubset(set("ATCG")):
        raise ValueError("Sequence contains non-canonical nucleotides.")
    if len(seq) < 2:
        raise ValueError("Sequence length must be at least 2 base pairs for NN analysis.")

    total_dH = 0.0
    total_dS = 0.0

    # 1. Sum nearest-neighbor dinucleotide step values
    for i in range(len(seq) - 1):
        dinucleotide = seq[i:i+2]
        total_dH += NN_PARAMS[dinucleotide]['dH']
        total_dS += NN_PARAMS[dinucleotide]['dS']

    # 2. Add terminal base pair initiation parameters
    for terminal_base in [seq[0], seq[-1]]:
        if terminal_base in ['G', 'C']:
            total_dH += INIT_GC['dH']
            total_dS += INIT_GC['dS']
        else:
            total_dH += INIT_AT['dH']
            total_dS += INIT_AT['dS']

    # Convert units: dH to cal/mol, concentration to Molar
    dH_cal = total_dH * 1000.0
    dS_cal = total_dS
    dnac_molar = dnac_nM * 1e-9
    na_molar = na_mM * 1e-3
    R = 1.9872  # Gas constant cal/(K*mol)

    # Check for self-complementarity factor
    comp_map = str.maketrans("ATCG", "TAGC")
    rev_comp = seq.translate(comp_map)[::-1]
    x_factor = 1.0 if seq == rev_comp else 4.0

    # Calculate Tm in Kelvin before salt correction
    tm_kelvin = dH_cal / (dS_cal + R * math.log(dnac_molar / x_factor))

    # Apply SantaLucia Salt Correction for [Na+]
    tm_celsius = tm_kelvin - 273.15 + (16.6 * math.log10(na_molar))

    return {
        "Enthalpy_dH_kcal": round(total_dH, 1),
        "Entropy_dS_cal": round(total_dS, 1),
        "Tm_Celsius": round(tm_celsius, 2)
    }


# ==============================================================================
# 3. IN-SILICO RESTRICTION DIGEST ENGINE
# ==============================================================================
def perform_restriction_digest(sequence, enzyme_name):
    """
    Scans DNA sequence for restriction enzyme cut sites and computes fragment sizes.
    """
    seq = sequence.upper().strip()
    if enzyme_name not in RESTRICTION_ENZYMES:
        raise ValueError(f"Enzyme '{enzyme_name}' is not in database. Available: {list(RESTRICTION_ENZYMES.keys())}")

    site, cut_offset = RESTRICTION_ENZYMES[enzyme_name]
    
    # Locate site start positions using lookahead search for overlapping motifs
    site_starts = [m.start() for m in re.finditer(f"(?={site})", seq)]
    cut_positions = [start + cut_offset for start in site_starts]

    # Calculate fragment lengths based on cut boundaries
    boundaries = [0] + cut_positions + [len(seq)]
    fragments = [boundaries[i+1] - boundaries[i] for i in range(len(boundaries)-1)]

    return {
        "Enzyme": enzyme_name,
        "Recognition_Site": site,
        "Cut_Positions": cut_positions,
        "Fragment_Count": len(fragments),
        "Fragment_Lengths_bp": fragments
    }


# ==============================================================================
# 4. SIMULATED AGAROSE GEL ELECTROPHORESIS VISUALIZATION
# ==============================================================================
def plot_simulated_gel(fragment_dict, ladder=STANDARD_LADDER):
    """
    Renders an in-silico agarose gel plot comparing digest fragments against a MW ladder.
    """
    plt.figure(figsize=(5, 8), facecolor='#111111')
    ax = plt.gca()
    ax.set_facecolor('#050505')

    ladder_x = 1.0
    sample_x = 2.0

    def bp_to_y(bp):
        return -np.log10(bp)

    # Plot Molecular Weight Ladder Lane
    for band in ladder:
        y_pos = bp_to_y(band)
        ax.hlines(y=y_pos, xmin=ladder_x - 0.25, xmax=ladder_x + 0.25, colors='#88CCFF', linewidth=2.5, alpha=0.85)
        ax.text(ladder_x - 0.45, y_pos, f"{band} bp", color='#88CCFF', fontsize=7, va='center', ha='right')

    # Plot Digested Sample Lane
    sample_fragments = fragment_dict["Fragment_Lengths_bp"]
    for frag in sample_fragments:
        y_pos = bp_to_y(frag)
        ax.hlines(y=y_pos, xmin=sample_x - 0.25, xmax=sample_x + 0.25, colors='#00FFCC', linewidth=3.0, alpha=0.95)
        ax.text(sample_x + 0.35, y_pos, f"{frag} bp", color='#00FFCC', fontsize=8, va='center', ha='left')

    ax.set_xlim(0, 3.2)
    ax.set_xticks([ladder_x, sample_x])
    ax.set_xticklabels(['MW Ladder', f'Digest ({fragment_dict["Enzyme"]})'], color='white', fontweight='bold')
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('#444444')
    ax.spines['left'].set_color('#444444')

    plt.title(f"In-Silico Agarose Gel (0.8% Matrix)\nTarget: {fragment_dict['Enzyme']} Digest", color='white', pad=15)
    plt.tight_layout()
    plt.show()


# ==============================================================================
# 5. CLI EXECUTION INTERFACE
# ==============================================================================
if __name__ == "__main__":
    print("==========================================================")
    print(" EpiTarget-BioCAD: Digest & Thermo Engine ")
    print("==========================================================")

    # Standard default sequence matching the demo run
    default_seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGA"

    input_seq = input("Enter DNA sequence (Press Enter to use default demo): ").strip()
    if not input_seq:
        input_seq = default_seq

    print(f"\nAnalyzing Sequence (Length: {len(input_seq)} bp)...")

    # 1. Thermodynamic Melting Temp Calculation
    try:
        thermo_res = calculate_tm_nearest_neighbor(input_seq)
        print("\n--- THERMODYNAMIC METRICS ---")
        print(f"Nearest-Neighbor dH° : {thermo_res['Enthalpy_dH_kcal']} kcal/mol")
        print(f"Nearest-Neighbor dS° : {thermo_res['Entropy_dS_cal']} cal/(K*mol)")
        print(f"Calculated Tm (50mM Na+) : {thermo_res['Tm_Celsius']} °C")
    except ValueError as e:
        print(f"Thermodynamic Error: {e}")

    # 2. Restriction Digest Calculation
    print("\nAvailable Enzymes:", ", ".join(RESTRICTION_ENZYMES.keys()))
    chosen_enzyme = input("Choose enzyme for digest (default: EcoRI): ").strip()
    if not chosen_enzyme:
        chosen_enzyme = "EcoRI"

    try:
        digest_res = perform_restriction_digest(input_seq, chosen_enzyme)
        print(f"\n--- RESTRICTION DIGEST RESULTS ({chosen_enzyme}) ---")
        print(f"Recognition Motif : {digest_res['Recognition_Site']}")
        print(f"Cut Indices       : {digest_res['Cut_Positions']}")
        print(f"Total Fragments   : {digest_res['Fragment_Count']}")
        print(f"Fragment Sizes    : {digest_res['Fragment_Lengths_bp']} bp")

        # 3. Gel Plotting
        plot_simulated_gel(digest_res)

    except ValueError as e:
        print(f"Digest Error: {e}")