# Definition of functions for LAMP sequences
"""
Created on Fri Jan 24 16:58:05 2025
Version : January 24th 2025
Version : April 16th 2025
@author: Arnaud Buhot
"""
import numpy as np
import matplotlib.pyplot as plt

# %% Calculation of DS from DH and DG and DG from DH and DS at 37°C


def DStheo(DH, DG):
    # Entropy calculation from Enthalpy and Free energy at 37°C
    DS = (DH-DG)/0.31015
    return DS


def DGtheo(DH, DS):
    # Free energy at 37°C from Enthalpy and Entropy
    DG = DH-DS*0.31015
    return DG


# %% Complementary sequence of dna
def DNA_comp(dna):
    """
    DNA_comp: Complementary sequence of dna

    Parameters
    ----------
    dna : str
        DNA sequence with A, T, C, G letters from 5' to 3' ends.
        Transform in capital letters if necessary
        Check for consistancy (only A, T, C or G letters)

    Returns
    -------
    str
        Complementary sequence of dna from 5' to 3' end.

    """
    # Check errors in DNA sequence
    dna_letters = {'A', 'C', 'G', 'T'}
    dna = dna.upper()
    for let in dna:
        if let not in dna_letters:
            return print(f"Invalid DNA letter found: {let}")

    # Define the base pairing rules
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    # Generate the complementary sequence
    dna_comp = "".join(comp[base] for base in dna)
    dna_comp = dna_comp[::-1]
    # Return the complementary sequence in reverse (antiparallel 5'--3')
    return dna_comp


# %% Thermodynamic parameters for melting analysis
# Nearest-neighbor parameters for DNA pairs in kcal/mol and cal/mol/K
# Reference: Hicks and SantaLucia Annual. Review (2004) Table 1. Errors?
# First NN 5'-3' and second NN 3'-5'
NN_params = {
    # No mismatch SantaLucia Proc. Natl. Acad. Sci. USA, 1998, 95, 1460–1465
    ('AA', 'TT'): {'H': -7.9, 'S': -22.2, 'G': -1.00},
    ('TT', 'AA'): {'H': -7.9, 'S': -22.2, 'G': -1.00},
    ('AT', 'TA'): {'H': -7.2, 'S': -20.4, 'G': -0.88},
    ('TA', 'AT'): {'H': -7.2, 'S': -21.3, 'G': -0.58},
    ('CA', 'GT'): {'H': -8.5, 'S': -22.7, 'G': -1.45},
    ('TG', 'AC'): {'H': -8.5, 'S': -22.7, 'G': -1.45},
    ('GT', 'CA'): {'H': -8.4, 'S': -22.4, 'G': -1.44},
    ('AC', 'TG'): {'H': -8.4, 'S': -22.4, 'G': -1.44},
    ('CT', 'GA'): {'H': -7.8, 'S': -21.0, 'G': -1.28},
    ('AG', 'TC'): {'H': -7.8, 'S': -21.0, 'G': -1.28},
    ('GA', 'CT'): {'H': -8.2, 'S': -22.2, 'G': -1.30},
    ('TC', 'AG'): {'H': -8.2, 'S': -22.2, 'G': -1.30},
    ('CG', 'GC'): {'H': -10.6, 'S': -27.2, 'G': -2.17},
    ('GC', 'CG'): {'H': -9.8, 'S': -24.4, 'G': -2.24},
    ('GG', 'CC'): {'H': -8.0, 'S': -19.9, 'G': -1.84},
    ('CC', 'GG'): {'H': -8.0, 'S': -19.9, 'G': -1.84},
    # Mismatch G-T Biochemistry 1997, 36, 10581-10594
    ('AG', 'TT'): {'H': 1.0, 'S': 0.9, 'G': 0.71},
    ('TT', 'GA'): {'H': 1.0, 'S': 0.9, 'G': 0.71},
    ('AT', 'TG'): {'H': -2.5, 'S': -8.3, 'G': 0.07},
    ('GT', 'TA'): {'H': -2.5, 'S': -8.3, 'G': 0.07},
    ('CG', 'GT'): {'H': -4.1, 'S': -11.7, 'G': -0.47},
    ('TG', 'GC'): {'H': -4.1, 'S': -11.7, 'G': -0.47},
    ('CT', 'GG'): {'H': -2.8, 'S': -8.0, 'G': -0.32},
    ('GG', 'TC'): {'H': -2.8, 'S': -8.0, 'G': -0.32},
    ('GG', 'CT'): {'H': 3.3, 'S': 10.4, 'G': 0.08},
    ('TC', 'GG'): {'H': 3.3, 'S': 10.4, 'G': 0.08},
    ('GG', 'TT'): {'H': 5.8, 'S': 16.3, 'G': 0.74},
    ('TT', 'GG'): {'H': 5.8, 'S': 16.3, 'G': 0.74},
    ('GT', 'CG'): {'H': -4.4, 'S': -12.3, 'G': -0.59},
    ('GC', 'TG'): {'H': -4.4, 'S': -12.3, 'G': -0.59},
    ('GT', 'TG'): {'H': 4.1, 'S': 9.5, 'G': 1.15},
    ('TG', 'AT'): {'H': -0.1, 'S': -1.7, 'G': 0.43},
    ('TA', 'GT'): {'H': -0.1, 'S': -1.7, 'G': 0.43},
    ('TG', 'GT'): {'H': -1.4, 'S': -6.2, 'G': 0.52},
    ('TT', 'AG'): {'H': -1.3, 'S': -5.3, 'G': 0.34},
    ('GA', 'TT'): {'H': -1.3, 'S': -5.3, 'G': 0.34},
    # Mismatch A-C Biochemistry 1998, 37, 9435-9444
    ('AA', 'TC'): {'H': 2.3, 'S': 4.6, 'G': 0.88},
    ('CT', 'AA'): {'H': 2.3, 'S': 4.6, 'G': 0.88},
    ('AC', 'TA'): {'H': 5.3, 'S': 14.6, 'G': 0.77},
    ('AT', 'CA'): {'H': 5.3, 'S': 14.6, 'G': 0.77},
    ('CA', 'GC'): {'H': 1.9, 'S': 3.7, 'G': 0.75},
    ('CG', 'AC'): {'H': 1.9, 'S': 3.7, 'G': 0.75},
    ('CC', 'GA'): {'H': 0.6, 'S': -0.6, 'G': 0.79},
    ('AG', 'CC'): {'H': 0.6, 'S': -0.6, 'G': 0.79},
    ('GA', 'CC'): {'H': 5.2, 'S': 14.2, 'G': 0.81},
    ('CC', 'AG'): {'H': 5.2, 'S': 14.2, 'G': 0.81},
    ('GC', 'CA'): {'H': -0.7, 'S': -3.8, 'G': 0.47},
    ('AC', 'CG'): {'H': -0.7, 'S': -3.8, 'G': 0.47},
    ('TA', 'AC'): {'H': 3.4, 'S': 8.0, 'G': 0.92},
    ('CA', 'AT'): {'H': 3.4, 'S': 8.0, 'G': 0.92},
    ('TC', 'AA'): {'H': 7.6, 'S': 20.2, 'G': 1.33},
    ('AA', 'CT'): {'H': 7.6, 'S': 20.2, 'G': 1.33},
    # Mismatch C-T Nucleic Acids Research, 1998, 26, 2694–2701
    ('AC', 'TT'): {'H': 0.7, 'S': 0.2, 'G': 0.64},
    ('TT', 'CA'): {'H': 0.7, 'S': 0.2, 'G': 0.64},
    ('AT', 'TC'): {'H': -1.2, 'S': -6.2, 'G': 0.73},
    ('CT', 'TA'): {'H': -1.2, 'S': -6.2, 'G': 0.73},
    ('CC', 'GT'): {'H': -0.8, 'S': -4.5, 'G': 0.62},
    ('TG', 'CC'): {'H': -0.8, 'S': -4.5, 'G': 0.62},
    ('CT', 'GC'): {'H': -1.5, 'S': -6.1, 'G': 0.40},
    ('CG', 'TC'): {'H': -1.5, 'S': -6.1, 'G': 0.40},
    ('GC', 'CT'): {'H': 2.3, 'S': 5.4, 'G': 0.62},
    ('TC', 'CG'): {'H': 2.3, 'S': 5.4, 'G': 0.62},
    ('GT', 'CC'): {'H': 5.2, 'S': 13.5, 'G': 0.98},
    ('CC', 'TG'): {'H': 5.2, 'S': 13.5, 'G': 0.98},
    ('TC', 'AT'): {'H': 1.2, 'S': 0.7, 'G': 0.97},
    ('TA', 'CT'): {'H': 1.2, 'S': 0.7, 'G': 0.97},
    ('TT', 'AC'): {'H': 1.0, 'S': 0.7, 'G': 0.75},
    ('CA', 'TT'): {'H': 1.0, 'S': 0.7, 'G': 0.75},
    # Mismatch G-A Biochemistry 1998, 37, 2170-2179
    ('AA', 'TG'): {'H': -0.6, 'S': -2.3, 'G': 0.14},
    ('GT', 'AA'): {'H': -0.6, 'S': -2.3, 'G': 0.14},
    ('AG', 'TA'): {'H': -0.7, 'S': -2.3, 'G': 0.02},
    ('AT', 'GA'): {'H': -0.7, 'S': -2.3, 'G': 0.02},
    ('CA', 'GG'): {'H': -0.7, 'S': -2.3, 'G': 0.03},
    ('GG', 'AC'): {'H': -0.7, 'S': -2.3, 'G': 0.03},
    ('CG', 'GA'): {'H': -4.0, 'S': -13.2, 'G': 0.11},
    ('AG', 'GC'): {'H': -4.0, 'S': -13.2, 'G': 0.11},
    ('GA', 'CG'): {'H': -0.6, 'S': -1.0, 'G': -0.25},
    ('GC', 'AG'): {'H': -0.6, 'S': -1.0, 'G': -0.25},
    ('GG', 'CA'): {'H': 0.5, 'S': 3.2, 'G': -0.52},
    ('AC', 'GG'): {'H': 0.5, 'S': 3.2, 'G': -0.52},
    ('TA', 'AG'): {'H': 0.7, 'S': 0.7, 'G': 0.42},
    ('GA', 'AT'): {'H': 0.7, 'S': 0.7, 'G': 0.42},
    ('TG', 'AA'): {'H': 3.0, 'S': 7.4, 'G': 0.74},
    ('AA', 'GT'): {'H': 3.0, 'S': 7.4, 'G': 0.74},
    # Mismatch A-A Biochemistry 1999, 38, 3468-3477
    ('AA', 'TA'): {'H': 1.2, 'S': 1.7, 'G': 0.61},
    ('AT', 'AA'): {'H': 1.2, 'S': 1.7, 'G': 0.61},
    ('CA', 'GA'): {'H': -0.9, 'S': -4.2, 'G': 0.43},
    ('AG', 'AC'): {'H': -0.9, 'S': -4.2, 'G': 0.43},
    ('GA', 'CA'): {'H': -2.9, 'S': -9.8, 'G': 0.17},
    ('AC', 'AG'): {'H': -2.9, 'S': -9.8, 'G': 0.17},
    ('TA', 'AA'): {'H': 4.7, 'S': 12.9, 'G': 0.69},
    ('AA', 'AT'): {'H': 4.7, 'S': 12.9, 'G': 0.69},
    # Mismatch C-C Biochemistry 1999, 38, 3468-3477
    ('AC', 'TC'): {'H': 0.0, 'S': -4.4, 'G': 1.33},
    ('CT', 'CA'): {'H': 0.0, 'S': -4.4, 'G': 1.33},
    ('CC', 'GC'): {'H': -1.5, 'S': -7.2, 'G': 0.70},
    ('CG', 'CC'): {'H': -1.5, 'S': -7.2, 'G': 0.70},
    ('GC', 'CC'): {'H': 3.6, 'S': 8.9, 'G': 0.79},
    ('CC', 'CG'): {'H': 3.6, 'S': 8.9, 'G': 0.79},
    ('TC', 'AC'): {'H': 6.1, 'S': 16.4, 'G': 1.05},
    ('CA', 'CT'): {'H': 6.1, 'S': 16.4, 'G': 1.05},
    # Mismatch G-G Biochemistry 1999, 38, 3468-3477
    ('AG', 'TG'): {'H': -3.1, 'S': -9.5, 'G': -0.13},
    ('GT', 'GA'): {'H': -3.1, 'S': -9.5, 'G': -0.13},
    ('CG', 'GG'): {'H': -4.9, 'S': -15.3, 'G': -0.11},
    ('GG', 'GC'): {'H': -4.9, 'S': -15.3, 'G': -0.11},
    ('GG', 'CG'): {'H': -6.0, 'S': -15.8, 'G': -1.11},
    ('GC', 'GG'): {'H': -6.0, 'S': -15.8, 'G': -1.11},
    ('TG', 'AG'): {'H': 1.6, 'S': 3.6, 'G': 0.44},
    ('GA', 'GT'): {'H': 1.6, 'S': 3.6, 'G': 0.44},
    # Mismatch T-T Biochemistry 1999, 38, 3468-3477
    ('AT', 'TT'): {'H': -2.7, 'S': -10.8, 'G': 0.69},
    ('TT', 'TA'): {'H': -2.7, 'S': -10.8, 'G': 0.69},
    ('CT', 'GT'): {'H': -5.0, 'S': -15.8, 'G': -0.12},
    ('TG', 'TC'): {'H': -5.0, 'S': -15.8, 'G': -0.12},
    ('GT', 'CT'): {'H': -2.2, 'S': -8.4, 'G': 0.45},
    ('TC', 'TG'): {'H': -2.2, 'S': -8.4, 'G': 0.45},
    ('TT', 'AT'): {'H': 0.2, 'S': -1.5, 'G': 0.68},
    ('TA', 'TT'): {'H': 0.2, 'S': -1.5, 'G': 0.68},
}

NN_init = {'Init': {'H': 0.2, 'S': -5.6, 'G': 1.96},
           # 'Init': {'H': 0.2, 'S': -5.7, 'G': 1.96},
           'Sym': {'H': 0.0, 'S': -1.4, 'G': 0.43},
           'A': {'H': 2.2, 'S': 6.9, 'G': 0.05},
           'T': {'H': 2.2, 'S': 6.9, 'G': 0.05},
           'C': {'H': 0.0, 'S': 0.0, 'G': 0.0},
           'G': {'H': 0.0, 'S': 0.0, 'G': 0.0}
           }

# NN_dend Bommarito Nucleic Acids Research, 2000, 28, 1929-1934
# 5' dangling end: First pair NN 5'-3' and second nucleotide NN 3'-5'
NN_dend5 = {('AA', 'T'): {'H': 0.2, 'S': 2.3, 'G': -0.51},
            ('AC', 'G'): {'H': -6.3, 'S': -17.1, 'G': -0.96},
            ('AG', 'C'): {'H': -3.7, 'S': -10.0, 'G': -0.58},
            ('AT', 'A'): {'H': -2.9, 'S': -7.6, 'G': -0.50},
            ('CA', 'T'): {'H': 0.6, 'S': 3.3, 'G': -0.42},
            ('CC', 'G'): {'H': -4.4, 'S': -12.6, 'G': -0.52},
            ('CG', 'C'): {'H': -4.0, 'S': -11.9, 'G': -0.34},
            ('CT', 'A'): {'H': -4.1, 'S': -13.0, 'G': -0.02},
            ('GA', 'T'): {'H': -1.1, 'S': -1.6, 'G': -0.62},
            ('GC', 'G'): {'H': -5.1, 'S': -14.0, 'G': -0.72},
            ('GG', 'C'): {'H': -3.9, 'S': -10.9, 'G': -0.56},
            ('GT', 'A'): {'H': -4.2, 'S': -15.0, 'G': 0.48},
            ('TA', 'T'): {'H': -6.9, 'S': -20.0, 'G': -0.71},
            ('TC', 'G'): {'H': -4.0, 'S': -10.9, 'G': -0.58},
            ('TG', 'C'): {'H': -4.9, 'S': -13.8, 'G': -0.61},
            ('TT', 'A'): {'H': -0.2, 'S': --0.5, 'G': -0.10},
            }

# 3' dangling end: First nucleotide NN 5'-3' and second pair NN 3'-5'
NN_dend3 = {('A', 'AT'): {'H': -0.7, 'S': -0.8, 'G': -0.48},
            ('C', 'AG'): {'H': -2.1, 'S': -3.9, 'G': -0.92},
            ('G', 'AC'): {'H': -5.9, 'S': -16.5, 'G': -0.82},
            ('T', 'AA'): {'H': -0.5, 'S': -1.1, 'G': -0.12},
            ('A', 'CT'): {'H': 4.4, 'S': 14.9, 'G': -0.19},
            ('C', 'CG'): {'H': -0.2, 'S': -0.1, 'G': -0.23},
            ('G', 'CC'): {'H': -2.6, 'S': -7.4, 'G': -0.31},
            ('T', 'CA'): {'H': 4.7, 'S': 14.2, 'G': 0.28},
            ('A', 'GT'): {'H': -1.6, 'S': -3.6, 'G': -0.50},
            ('C', 'GG'): {'H': -3.9, 'S': -11.2, 'G': -0.44},
            ('G', 'GC'): {'H': -3.2, 'S': -10.4, 'G': -0.01},
            ('T', 'GA'): {'H': -4.1, 'S': -13.1, 'G': -0.01},
            ('A', 'TT'): {'H': 2.9, 'S': 10.4, 'G': -0.29},
            ('C', 'TG'): {'H': -4.4, 'S': -13.1, 'G': -0.35},
            ('G', 'TC'): {'H': -5.2, 'S': -15.0, 'G': -0.52},
            ('T', 'TA'): {'H': -3.8, 'S': -12.6, 'G': 0.13},
            }


# %% Create dictionary for Watson Crick Franklin pairs
# Data from NNDB https://rna.urmc.rochester.edu/NNDB/
bases = ['A', 'C', 'G', 'T']

# Read dg from file 'NNBD/dna_terminal_mismatch_dg.txt'
dna_wcf_dg = list()
with open('NNBD/dna_watson_crick_stack_dg.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#') or not line:
            continue  # Skip comment or empty lines
        parts = line.split('\t')
        dna_wcf_dg.append(parts)

# Read dg from file 'NNBD/dna_terminal_mismatch_dh.txt'
dna_wcf_dh = list()
with open('NNBD/dna_watson_crick_stack_dh.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#') or not line:
            continue  # Skip comment or empty lines
        parts = line.split('\t')
        dna_wcf_dh.append(parts)

# Create dictionary NNDB_termis for terminal mismatches
nline = 0
nl = len(dna_wcf_dh)
nd = int(nl/7)
NNDB_wcf = dict()

for ndd in range(nd):
    n11 = dna_wcf_dh[nline][0][0]
    nline += 1
    n21 = dna_wcf_dh[nline][0][0]
    nline += 2
    for n12 in bases:
        for i in range(4):
            n22 = bases[i]
            pair = (n11+n12, n21+n22)
            pairc = (n22+n21, n12+n11)
            if dna_wcf_dg[nline][i+1] != '.':
                dg = float(dna_wcf_dg[nline][i+1])
                dh = float(dna_wcf_dh[nline][i+1])
                NNDB_wcf[pair] = dict()
                NNDB_wcf[pair]['G'] = dg
                NNDB_wcf[pair]['H'] = dh
                NNDB_wcf[pair]['S'] = DStheo(dh, dg)
                NNDB_wcf[pairc] = dict()
                NNDB_wcf[pairc]['G'] = dg
                NNDB_wcf[pairc]['H'] = dh
                NNDB_wcf[pairc]['S'] = DStheo(dh, dg)
        nline += 1


# %% Create dictionary for terminal mismatches (double dangling ends)
# Data from NNDB https://rna.urmc.rochester.edu/NNDB/
bases = ['A', 'C', 'G', 'T']

# Read dg from file 'NNBD/dna_terminal_mismatch_dg.txt'
dna_tm_dg = list()
with open('NNBD/dna_terminal_mismatch_dg.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#') or not line:
            continue  # Skip comment or empty lines
        parts = line.split('\t')
        dna_tm_dg.append(parts)

# Read dg from file 'NNBD/dna_terminal_mismatch_dh.txt'
dna_tm_dh = list()
with open('NNBD/dna_terminal_mismatch_dh.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#') or not line:
            continue  # Skip comment or empty lines
        parts = line.split('\t')
        dna_tm_dh.append(parts)

# Create dictionary NNDB_termis for terminal mismatches
nline = 0
nl = len(dna_tm_dh)
nd = int(nl/7)
NNDB_termis = dict()

for ndd in range(nd):
    n11 = dna_tm_dh[nline][0][0]
    nline += 1
    n21 = dna_tm_dh[nline][0][0]
    nline += 2
    for n12 in bases:
        n1 = n11 + n12
        for i in range(4):
            n22 = bases[i]
            n2 = n21 + n22
            pair = (n1, n2)
            dg = float(dna_tm_dg[nline][i+1])
            dh = float(dna_tm_dh[nline][i+1])
            NNDB_termis[pair] = dict()
            NNDB_termis[pair]['G'] = dg
            NNDB_termis[pair]['H'] = dh
            NNDB_termis[pair]['S'] = DStheo(dh, dg)
        nline += 1


# %% Salt corrections Peyret and Owczarzy
# Tris buffer: Assume half of the Tris ions are free (true for all buffers?)
# dNTP presence: If dNTP < Mg2+, cMg = cMg - cdNTP
# Peyret salt correction from Owczarzy 2008.
# Peyret value 3.3 Eq.6 or Ashen value 3.79 Eq.7
def SaltP(nd, cNa, cMg):
    if cMg > 0.:
        DSs = 0.368*nd*np.log(cNa + 3.3*np.sqrt(cMg))
    elif cMg == 0.:
        if cNa > 0.:
            DSs = 0.368*nd*np.log(cNa)
        else:
            DSs = 0.
            print('Problem with salt concentration.')
            print('Na+ = {cNa} and Mg2+ = {cMg}')
    else:
        DSs = 0.
        print(f'Problem with negative Mg2+ concentration: {cMg}')
    return DSs


# Owczarzy salt correction Owczarzy 2008 Eq.4, 16, 18-20 and 22
def SaltO(DH, nd, fGC, cNa, cMg):
    a = 3.92e-5
    b = -9.11e-6
    c = 6.26e-5
    d = 1.42e-5
    e = -4.82e-4
    f = 5.25e-4
    g = 8.31e-5
    DS = 0.
    # Calculate ln[cMg]
    if cMg > 0.:
        lncMg = np.log(cMg)

    # Calculate Salt correction
    if cNa == 0.:
        if cMg == 0.:
            print(r'Problem with no salt: Na+ = {cNa} and Mg2+ = {cMg}')
        else:
            DS += DH*(a+b*lncMg)
            DS += DH*fGC*(c+d*lncMg)
            DS += DH*(e + f*lncMg + g*lncMg**2)/(2.*(nd-1))
    else:
        lncNa = np.log(cNa)
        Rval = np.sqrt(cMg)/cNa
        if Rval < 0.22:
            # Only Monovalent salt correction Eq.4 Owczarzy 2008
            DS += DH*((4.29*fGC-3.95)*1e-5*lncNa + 9.4e-6*lncNa**2)
        elif Rval < 6.0:
            a2 = 3.92e-5*(0.843 - 0.352*np.sqrt(cNa)*lncNa)
            d2 = 1.42e-5*(1.279 - 4.03e-3*lncNa - 8.03e-3*lncNa**2)
            g2 = 8.31e-5*(0.486 - 0.258*lncNa + 5.25e-3*lncNa**3)
            DS += DH*(a2+b*lncMg)
            DS += DH*fGC*(c+d2*lncMg)
            DS += DH*(e + f*lncMg + g2*lncMg**2)/(2.*(nd-1))
        else:
            DS += DH*(a+b*lncMg)
            DS += DH*fGC*(c+d*lncMg)
            DS += DH*(e + f*lncMg + g*lncMg**2)/(2.*(nd-1))
    return DS*1000.


# %% Thermodynamic parameters for DH, DS and DG
def ThermoParam(seq1, seq2, cNa=1., cMg=0., salt='Owczarzy'):
    """
    Thermodynamic parameters from hybridization of seq1 with seq2.
    Nearest Neighbor parameters from SantaLucia (see Dictionaries)
    Constraint: Developed for seq1 = F-FIP and seq2 = Fd-Q
    Thus, 5' end of seq1 should be complementary to 3' end of seq2

    Parameters
    ----------
    seq1 : str
        Sequence of DNA (generally FIP).
    seq2 : str
        Sequence of DNA (generally Fd).
    cNa : float, optional
        Na+ (or monovalent) salt concentration in M.
        The default value is 1. M.
    cMg : float, optional
        Mg2+ (or divalent) salt concentration in M.
        Should be corrected by cdNTP if present: cMg - cdNTP
        The default value is 0. M.
    salt: str, optional
        Salt correction 'Peyret' or 'Owczarzy'
        From Owczarzy et al., Biochemistry 2008, 47, 5336–5353
        'Peyret': Eq.(6)
        'Owczarzy': Eq.(22)
        The deflault value is 'Owczarzy'

    Returns
    -------
    DH : float
        Enthalpy in kcal/mol.
    DS : flaot
        Entropy in cal/mol/K.
    DG : flaot
        Free energy in kcal/mol.

    """
    # Default values: Na+ = 1M and Mg2+ = 0M
    # Default salt correction: Owczarzy
    # seq1 = F-FIP and seq2 = Fd-Q to have perfect quenching we assume that
    # seq1 and seq2 form a duplex with fully hybridized at 5'-seq1 and 3'-seq2
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    print(f"Sequence 1 : 5'-{seq1}-3'")
    # Sequence  2: 3'-5'
    seq2 = seq2[::-1]
    print(f"Sequence 2 : 3'-{seq2}-5'")
    # Check that duplex is hybridized at 5'-seq1 and 3'-seq2:
    if seq1[0] != comp[seq2[0]]:
        print(f'Problem of hybridization: {seq1[0]} {comp[seq2[0]]}')
    # Smallest seq length
    nb1 = len(seq1)
    nb2 = len(seq2)
    print(f'Lengths (seq1 and seq2): {nb1} and {nb2} bases')
    nb = min(nb1, nb2)
    # Determine existence of dangling ends
    nd = nb
    while seq1[nd-1] != comp[seq2[nd-1]]:
        # print(seq1[nd-1], comp[seq2[nd-1]])
        nd -= 1
        if nd < 0:
            print('Problem no hybridization duplex')
            break
    print(f'Duplex length: {nd} base pairs')

    # Determine mismatches poitions
    nm = 0
    nFC = 0
    mm = list()
    for i in range(nd):
        if seq1[i] != comp[seq2[i]]:
            nm += 1
            mm.append(i+1)
            print(f'Mismatch {seq1[i]}-{seq2[i]} in position {i+1}')
        elif seq1[i] in ['C', 'G'] and seq2[i] in ['C', 'G']:
            # Count for GC pairs
            nFC += 1
    if nm > 0:
        print(f'Total number of mismatches: {nm}')
        print(f'List of mismatch positions: {mm}')
    elif nm == 0:
        print('No mismatches')
    else:
        print(f'Error in the determination of mismatches {nm}')
    fGC = float(nFC)/float(nd)
    print(f'Fraction of GC pairs in duplex: {100.*fGC:.1f}%')
    # TODO: Check positions of mismatches : not to close to end of duplex
    # not two consecutive mismatches (only single point mutations)

    # Add Initiation DH and DS
    DH = NN_init['Init']['H']
    DS = NN_init['Init']['S']
    DG = NN_init['Init']['G']
    # Add Terminal A or T penalty at seq1 5' and seq2 3' ends
    DH += NN_init[seq1[0]]['H']
    DS += NN_init[seq1[0]]['S']
    DG += NN_init[seq1[0]]['G']
    DH += NN_init[seq1[nd-1]]['H']
    DS += NN_init[seq1[nd-1]]['S']
    DG += NN_init[seq1[nd-1]]['G']
    # Add symetry penalty
    if seq1 == seq2:
        DH += NN_init['Sym']['H']
        DS += NN_init['Sym']['S']
        DG += NN_init['Sym']['G']
        print('Symetry penalty added')

    # Add Nearest-Neighbor DH and DS
    for i in range(nd-1):
        pair1 = seq1[i:i+2]
        pair2 = seq2[i:i+2]
        nn = (pair1, pair2)
        # print(nn)
        if nn in NN_params:
            DH += NN_params[nn]['H']
            DS += NN_params[nn]['S']
            DG += NN_params[nn]['G']
        else:
            print('Problem of sequence')
    # Add dangling ends or A-T penalty in seq1 3' and seq2 5' ends
    if nd < nb:
        print('Presence of two dangling ends')
        # From NNDB_termis
        pair1 = seq1[nd-1:nd+1]
        pair2 = seq2[nd-1:nd+1]
        nn = (pair1, pair2)
        print(f'Dangling end pairs: {nn}')
        if nn in NN_params:
            DH += NNDB_termis[nn]['H']
            DS += NNDB_termis[nn]['S']
            DG += NNDB_termis[nn]['G']
        else:
            print('Problem of sequence')
    elif nb1 > nb2:
        print("Presence of a 3' dangling end")
        pair1 = seq1[nd]+seq1[nd-1]
        pair2 = seq2[nd-1]
        nn = (pair2, pair1)
        print(f'Dangling end pairs: {nn}')
        if nn in NN_dend3:
            DH += NN_dend3[nn]['H']
            DS += NN_dend3[nn]['S']
            DG += NN_dend3[nn]['G']
            print(NN_dend3[nn])
        else:
            print('Problem of sequence')
    elif nb1 < nb2:
        print("Presence of a 5' dangling end")
        pair1 = seq1[nd-1]
        pair2 = seq2[nd]+seq2[nd-1]
        nn = (pair2, pair1)
        print(f'Dangling end pairs: {nn}')
        if nn in NN_dend5:
            DH += NN_dend5[nn]['H']
            DS += NN_dend5[nn]['S']
            DG += NN_dend5[nn]['G']
            print(NN_dend5[nn])
        else:
            print('Problem of sequence')
    elif nb1 == nb2:
        DH += NN_init[seq1[nd-1]]['H']
        DS += NN_init[seq1[nd-1]]['S']
        DG += NN_init[seq1[nd-1]]['G']
        print('No dangling end')

    # Salt correction from Owczarzy 2008
    # Peyret Eq.6 (Na and Mg)
    # Owczarzy Eq.4, 16, 18-20 and 22
    if salt == 'Peyret':
        DScorr = SaltP(nd, cNa, cMg)
    elif salt == 'Owczarzy':
        DScorr = SaltO(DH, nd, fGC, cNa, cMg)
    DS += DScorr                # Salt correction in cal/mol/K
    DG -= 310.15*DScorr/1000.   # T*DS with T = 310.15 K and DG in kcal/mol
    return DH, DS, DG, nd, fGC


def Tmelt(DH, DS, cA=1.0e-6, cB=1.0e-6):
    """
    Tmelt: Melting temperature from Enthalpy and Entropy

    Parameters
    ----------
    DH : float
        Enthalpy in kcal/mol.
    DS : float
        Entropy in cal/mol/K.
    cA : float, optional
        Concentration of sequence 1 in M.
        The default value is 1.0e-6 M or 1 uM.
    cB : float, optional
        Concentration of sequence 2 in M.
        The default value is 1.0e-6 M or 1 uM.

    Returns
    -------
    tm : float
        Melting temperature in °C.

    """
    # Tm calculation (using 1uM strand concentration as default)
    if cA > cB:
        cDNA = cA-cB/2.
    else:
        cDNA = cB-cA/2.
    R = 1.9872  # gas constant in cal/mol·K
    tm = 1000.*DH/(DS + R*np.log(cDNA))-273.15
    return tm


def Theta(DH, DS, T, cA=1.0e-6, cB=1.0e-6):
    R = 1.9872          # gas constant in cal/mol·K
    TK = T + 273.15     # Temperature in Kelvin
    KD = np.exp(-(1000.*DH-TK*DS)/(R*TK))
    theta = 1. + KD*(cA+cB)
    theta -= np.sqrt(1 + 2*KD*(cA+cB) + (KD*(cA-cB))**2)
    if cA > cB:
        theta = theta/(2*KD*cB)
    else:
        theta = theta/(2*KD*cA)
    return theta


def dTheta(DH, DS, T, cA=1.0e-6, cB=1.0e-6):
    R = 1.9872          # gas constant in cal/mol·K
    TK = T + 273.15     # Temperature in Kelvin
    KD = np.exp(-(1000.*DH-TK*DS)/(R*TK))
    theta = Theta(DH, DS, T, cA, cB)
    sq = np.sqrt(1 + 2*KD*(cA+cB) + (KD*(cA-cB))**2)
    dtheta = theta*1000.*DH/(sq*R*TK*TK)
    return dtheta


def PlotTmelt(DH, DS, cA=1.0e-6, cB=1.0e-6):
    Temp = [45. + 0.1*i for i in range(451)]
    theta = [Theta(DH, DS, T, cA, cB) for T in Temp]
    dtheta = [dTheta(DH, DS, T, cA, cB) for T in Temp]
    # Create the plot
    plt.plot(Temp, theta, color='Blue')
    # Labels and title
    plt.xlabel('Temperature (°C)')
    plt.xlim(45, 90)
    plt.ylabel(r'$\theta$')
    plt.ylim(0, 1)
    plt.title('Hybridized fraction and derivative')
    plt.show()
    # Create the plot
    plt.plot(Temp, dtheta, color='Green')
    # Labels and title
    plt.xlabel('Temperature (°C)')
    plt.xlim(45, 90)
    plt.ylabel(r'$d\theta / dT$')
    plt.title('Hybridized fraction and derivative')
    plt.show()


def Hybridization(seq1, seq2, cA=1.0e-6, cB=1.0e-6, cNa=1., cMg=0.):
    # Calculation of thermodynamic parameters
    DH, DS, DG, nd, fGC = ThermoParam(seq1, seq2, cNa, cMg)
    print(f"ΔH: {DH:.2f} kcal/mol")
    print(f"ΔS: {DS:.2f} cal/mol/K")
    print(f"ΔS(bis): {DStheo(DH, DG):.2f} cal/mol/K")
    print(f"ΔG(37°C): {DG:.2f} kcal/mol")
    print(f"ΔG(bis): {DGtheo(DH, DS):.2f} kcal/mol")

    Tm = Tmelt(DH, DS, cA, cB)
    print(f"Melting Temperature (Tm): {Tm:.1f}°C")

    Tmbis = Tmelt(DH, DStheo(DH, DG), cA, cB)
    print(f"Melting Temperature (Tmbis): {Tmbis:.1f}°C")

    # Plot of hybridized fraction vs temperature (range [25, 95] in °C)
    PlotTmelt(DH, DS, cA, cB)
    return DH, DS, DG, Tm


# %% LAMP primers and dumbbell
def FIP(F1, F2):
    # Give the FIP sequence FIP = F1c + F2
    F1c = DNA_comp(F1)
    FIP = F1c + F2
    return FIP


# Give the dumbbel sequence for a set of primers B1, B2, F1, F2
def Dumbbell(B1, B2, F1, F2):
    # Dumbbell = F1c-F2-F1-B1c-B2c-B1
    B1c = DNA_comp(B1)
    B2c = DNA_comp(B2)
    F1c = DNA_comp(F1)
    dumb = F1c+F2+F1+B1c+B2c+B1
    return dumb
