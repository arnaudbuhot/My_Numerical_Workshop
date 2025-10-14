# Nearest Neighbor thermodynamics and LAMP primers and dumbbell
# Example of usage of functions from Definitions
"""
LAMP sequences, NN thermodynamical parameters
Author: A. Buhot
Version January 24th 2025
Version April 16th 2025
"""
import matplotlib.pyplot as plt
import numpy as np
from Definitions import DNA_comp, Dumbbell, Hybridization, SaltO, SaltP
from Definitions import ThermoParam, Tmelt, PlotTmelt, DStheo, DGtheo

# %% Complementary strand from a DNA sequence
seq1 = "ATGCGTAC"
help(DNA_comp)      # Help for DNA_comp function
print(f"DNA Sequence : 5'-{seq1}-3'")
seq2 = DNA_comp(seq1)
print(f"Complementary: 5'-{seq2}-3'")

# Transform in upper letters inside function but do not correct seq1
seq1 = "atcg"
print(f"DNA Sequence : 5'-{seq1}-3'")
seq2 = DNA_comp(seq1)
print(f"Complementary: 5'-{seq2}-3'")

# Alert if wrong letters included
seq1 = 'atcgb'
seq2 = DNA_comp(seq1)

# %% DS vs DG(37°C) and salt correction coefficient
DH = -4.2
DG = 0.48
DS = DStheo(-4.2, 0.48)
print(f'DS = {DS:.2f}')
DG = DGtheo(DH, DS)
print(f'DG = {DG:.2f}')

SaltDS = 0.368      # in cal/mol/K
SaltDG = -0.114     # in kcal/mol
Temp = 310.15           # in Kelvin at 37°C
SaltcorrDS = 1000.*SaltDG/Temp
print(f'SaltcorrDS = {SaltcorrDS:.5f}')
SaltcorrDG = Temp*SaltDS/1000.
print(f'SaltcorrDG = {SaltcorrDG:.5f}')

# %% Calculation of thermodynamic parameters
seq1 = "ATCGATCG"
seq1 = 'GATGAATATGGAGGCGGCCAA'
seq2 = DNA_comp(seq1)
DH, DS, DG, nd, fGC = ThermoParam(seq1, seq2)
print(f"ΔH: {DH:.2f} kcal/mol")
print(f"ΔS: {DS:.2f} cal/mol/K")
print(f"ΔS(bis): {DStheo(DH, DG):.2f} cal/mol/K")
print(f"ΔG(37°C): {DG:.2f} kcal/mol")
print(f"ΔG(bis): {DGtheo(DH, DS):.2f} kcal/mol")

Tm = Tmelt(DH, DS)
print(f"Melting Temperature (Tm): {Tm:.1f}°C")

Tmbis = Tmelt(DH, DStheo(DH, DG))
print(f"Melting Temperature (Tmbis): {Tmbis:.1f}°C")

# Plot of hybridized fraction vs temperature (range [25, 95] in °C)
PlotTmelt(DH, DS)

# %% Idem with simple function Hybridization
seq1 = "ATCGATCG"
seq1 = 'GATGAATATGGAGGCGGCCAA'
DH, DS, DG, Tm = Hybridization(seq1, seq2)

# %% Example usage
seq1 = "ATGCATGCATGCATGCATGCATGC"
# seq1 = "ATCGATCG"
seq2 = DNA_comp(seq1)
seq2 = 'GCATGCATGCATGCATGCATGCAT'   # Comp seq1
seq2 = 'AAATGCATGCATGCATGCATGCAT'   # 5' end GC remplaced by AA
seq2 = 'GTGCATGTATGCATGCATGCATGCAT'   # Pos 6 (5' end): Mutation C -> T
print(seq2)
DH, DS, DG, nd, fGC = ThermoParam(seq1, seq2)
Tm = Tmelt(DH, DS)
print(f"ΔH: {DH:.2f} kcal/mol")
print(f"ΔS: {DS:.2f} cal/mol/K")
print(f"ΔG(37°C): {DG:.2f} kcal/mol")
print(f"ΔS(bis): {DStheo(DH, DG):.2f} cal/mol/K")
print(f"ΔG(bis): {DGtheo(DH, DS):.2f} kcal/mol")
print(f"Melting Temperature (Tm): {Tm:.1f} °C")
Tm = Tmelt(DH, DStheo(DH, DG))
print(f"Melting Temperature (Tm) bis: {Tm:.1f} °C")
PlotTmelt(DH, DS)

# %% Zika Virus From Pierre Garneret PhD Table 6.5
ZIKV_F3 = 'CAGCGTTCACATTCACCA'
ZIKV_FIP = 'TTGCATGTCCACCGCCATAGTCACAGTGGAGGTACAG'
ZIKV_B3 = 'CCGACTCCTATGACAATGTAAG'
ZIKV_BIP = 'GGGAGGTTGATAACCGCTAACCGGATCAAGTTCCAGCATCAT'
ZIKV_LF = 'CTGAGCTGGAACCTTGCA'
ZIKV_LB = 'ATCACTGAAAGCACTGAGAACT'

# %% West Nile Virus sequences
# From Ball et al. Anal. Chem. 88, 3562 2016)
WNV_F3 = 'TGGATTTGGTTCTCGAAGG'
WNV_F2 = 'CAGCTGCGTGACTATCATGT'
WNV_F1c = 'TTGGCCGCCTCCATATTCATCA'
WNV_B3 = 'GGTCAGCACGTTTGTCATT'
WNV_B2 = 'TGAGCTTCTCCCATGGTCG'
WNV_B1c = 'TGCTATTTGGCTACCGTCAGCG'
WNV_LF = 'CATCGATGGTAGGCTTGTC'
WNV_LB = 'TCTCCACCAAAGCTGCGT'

print("WNV_B1c")
WNV_B1 = DNA_comp(WNV_B1c)

print("WNV_F1c")
WNV_F1 = DNA_comp(WNV_F1c)

print("WNV_F2c")
WNV_F2c = DNA_comp(WNV_F2)

# Give the FIP and BIP sequences FIP = F1c + F2 and BIP = B1c + B2
WNV_FIP = WNV_F1c + WNV_F2
WNV_BIP = WNV_B1c + WNV_B2
WNV_Dumbbell = Dumbbell(WNV_B1, WNV_B2, WNV_F1, WNV_F2)

print(f"FIP sequence: 5'-{WNV_FIP}-3'")
print(f"FIP length: {len(WNV_FIP)}")
print(f"BIP sequence: 5'-{WNV_BIP}-3'")
print(f"BIP length: {len(WNV_BIP)}")
print(f"Dumbbell sequence: 5'-{WNV_Dumbbell}-3'")
print(f"Dumbbell length: {len(WNV_Dumbbell)}")

# %% Fd sequences (Complementary to FIP)
# DE3 = non-complementary sequence
# TmCp for 157mM NaCl and 5 mM Mg2+ at 1uM for both seq.
WNV_Fd21 = 'GATGAATATGGAGGCGGCCAA'      # 73.8
WNV_Fd21m16 = 'GATGAATATGGAGGCAGCCAA'   # 63.5 Mutation G en A idem Ball et al.
WNV_Fd21m6 = 'GATGAGTATGGAGGCGGCCAA'    # 71.3 Mutation A en G en position 6
WNV_Fd10 = 'AGGCGGCCAA'                 # 60.5
WNV_Fd10DE3 = 'ACCAGGCGGCCAA'           # 62.9
WNV_Fd9DE1 = 'TGGCGGCCAA'               # 56.8
WNV_Fd9DE3 = 'CCTGGCGGCCAA'             # 57.2
WNV_Fd8DE2 = 'TAGCGGCCAA'
WNV_Fd8DE3 = 'CTCGCGGCCAA'
WNV_Fd7DE3 = 'TCCCGGCCAA'
# WNV_Fd8DE2, WNV_Fd9DE1 and WNV_Fd10 already present in Ball et al.
Fd_list = [WNV_Fd21, WNV_Fd21m16, WNV_Fd21m6, WNV_Fd10, WNV_Fd10DE3,
           WNV_Fd9DE1, WNV_Fd9DE3, WNV_Fd8DE2, WNV_Fd8DE3, WNV_Fd7DE3]
ThermoParam = dict()
for Fd in Fd_list:
    ThermoParam[Fd] = Hybridization(WNV_FIP, Fd, cNa=1., cMg=0.0)

for Fd in Fd_list:
    print(f'Tm({str(Fd)}) = {ThermoParam[Fd][3]}')

# %% Salt concentraiton effects
seq1 = WNV_FIP
seq2 = WNV_Fd21
DH, DS, DG, nd, fGC = ThermoParam(seq1, seq2)
Tm = Tmelt(DH, DS)
print(f'Tm = {Tm:.1f}°C')
cNa = 0.1
cMg = 0.005
DScor = SaltO(DH, nd, fGC, cNa, cMg)
print(DScor)
Tm = Tmelt(DH, DS+DScor)
print(f'Tm = {Tm:.1f}°C')
DScor = SaltP(nd, cNa, cMg)
print(DScor)
Tm = Tmelt(DH, DS+DScor)
print(f'Tm = {Tm:.1f}°C')

# Na+ concentration below which R > 6.0 for cMg = 10mM
cNa0 = np.sqrt(0.01)/6.0
print(f'Concentration Na+ = {1000.*cNa0:.1f} mM')

cNa = 1.
cMg1 = (0.22*cNa)**2
cMg2 = (6.0*cNa)**2
print(cMg1, cMg2)
cMg_list = np.linspace(0., 0.01, 81)
TmO = list()
TmP = list()
for cMg in cMg_list:
    DScor = SaltO(DH, nd, fGC, cNa, cMg)
    TmO.append(Tmelt(DH, DS+DScor))
    DScor = SaltP(nd, cNa, cMg)
    TmP.append(Tmelt(DH, DS+DScor))

cMg_list = cMg_list*1000.
xmin = np.min(cMg_list)
xmax = np.max(cMg_list)
xlim = [xmin, xmax]
# Create the plot
plt.plot(cMg_list, TmO, color='Blue', label='Owczarzy')
plt.plot(cMg_list, TmP, color='Green', label='Peyret')
plt.xlim(xlim)
plt.legend()
# Labels and title
plt.xlabel(r'$Mg^{2+}$ concentration (in mM)')
plt.ylabel('Melting temperature (in °C)')
plt.title(r'Owczarzy vs Peyret for $[Na^+] = 1 M$')
plt.show()

cNa = 0.1
cMg1 = (0.22*cNa)**2
cMg2 = (6.0*cNa)**2
print(cMg1, cMg2)
cMg_list = np.linspace(0., 0.01, 81)
TmO = list()
TmP = list()
for cMg in cMg_list:
    DScor = SaltO(DH, nd, fGC, cNa, cMg)
    TmO.append(Tmelt(DH, DS+DScor))
    DScor = SaltP(nd, cNa, cMg)
    TmP.append(Tmelt(DH, DS+DScor))

cMg_list = cMg_list*1000.
xmin = np.min(cMg_list)
xmax = np.max(cMg_list)
xlim = [xmin, xmax]
# Create the plot
plt.plot(cMg_list, TmO, color='Blue', label='Owczarzy')
plt.plot(cMg_list, TmP, color='Green', label='Peyret')
plt.xlim(xlim)
plt.legend()
# Labels and title
plt.xlabel(r'$Mg^{2+}$ concentration (in mM)')
plt.ylabel('Melting temperature (in °C)')
plt.title(r'Owczarzy vs Peyret for $[Na^+] = 0.1 M$')
plt.show()

cNa = 0.01
cMg1 = (0.22*cNa)**2
cMg2 = (6.0*cNa)**2
print(cMg1, cMg2)
cMg_list = np.linspace(0., 0.01, 81)
TmO = list()
TmP = list()
for cMg in cMg_list:
    DScor = SaltO(DH, nd, fGC, cNa, cMg)
    TmO.append(Tmelt(DH, DS+DScor))
    DScor = SaltP(nd, cNa, cMg)
    TmP.append(Tmelt(DH, DS+DScor))

cMg_list = cMg_list*1000.
xmin = np.min(cMg_list)
xmax = np.max(cMg_list)
xlim = [xmin, xmax]
# Create the plot
plt.plot(cMg_list, TmO, color='Blue', label='Owczarzy')
plt.plot(cMg_list, TmP, color='Green', label='Peyret')
plt.xlim(xlim)
plt.legend()
# Labels and title
plt.xlabel(r'$Mg^{2+}$ concentration (in mM)')
plt.ylabel('Melting temperature (in °C)')
plt.title(r'Owczarzy vs Peyret for $[Na^+] = 0.01 M$')
plt.show()

# %% Thermodynamic parameters
FIP = 'TTGGCCGCCTCCATATTCATCACAGCTGCGTGACTATCATGT'
# DynaMelt Two-state hybridization melting
# Paramters Salt concentration 1M Na+ and no Mg2+, Strands 1uM
# TP_TwoState = (DG, DH, DS Tm)
print(Fd_list)
print('ThermoParam :', ThermoParam)
TP_TwoState = dict()
TP_TwoState[WNV_Fd21] = (-27.6, -168.2, -453.4, 74.6)
TP_TwoState[WNV_Fd21m16] = (-22.0, -147.1, -403.3, 66.2)
TP_TwoState[WNV_Fd21m6] = (-24.9, -154.6, -418.2, 71.7)
TP_TwoState[WNV_Fd10] = (-14.0, -69.4, -178.5, 59.3)
TP_TwoState[WNV_Fd10DE3] = (-14.4, -68.8, -175.2, 61.7)
TP_TwoState[WNV_Fd9DE1] = (-13.7, -78.3, -208.2, 55.3)
TP_TwoState[WNV_Fd9DE3] = (-13.7, -78.3, -208.2, 55.3)
TP_TwoState[WNV_Fd8DE2] = (-11.7, -66.5, -176.8, 48.1)
TP_TwoState[WNV_Fd8DE3] = (-11.4, -66.8, -178.6, 46.8)
TP_TwoState[WNV_Fd7DE3] = (-9.3, -55.0, -147.4, 36.5)
print('TP_TwoState:', TP_TwoState)

# Partition function model Dynamelt with strands 10uM concentration
# TP_Partition = (DG, DH, DH, Tm(Cp))
TP_Partition = dict()
TP_Partition[WNV_Fd21] = (-42.0, -187.5, -536.4, 82.1)
TP_Partition[WNV_Fd21m16] = (-34.6, -166.6, -486.8, 73.4)
TP_Partition[WNV_Fd21m6] = (-38.0, -174.1, -501.6, 80.3)
TP_Partition[WNV_Fd10] = (-18.7, -101.6, -306.1, 71.6)
TP_Partition[WNV_Fd10DE3] = (-19.2, -99.1, -296.0, 74.0)
TP_Partition[WNV_Fd9DE1] = (-18.8, -108.4, -331.4, 68.2)
TP_Partition[WNV_Fd9DE3] = (-19.1, -113.5, -348.7, 68.1)
TP_Partition[WNV_Fd8DE2] = (-15.8, -101.4, -316.1, 64.8)
TP_Partition[WNV_Fd8DE3] = (-16.1, -107.3, -336.5, 64.3)
TP_Partition[WNV_Fd7DE3] = (-15.5, -113.8, -362.8, 67.7)
print(TP_Partition)

DG1 = [ThermoParam[Fd][2] for Fd in Fd_list]
DG2 = [TP_TwoState[Fd][0] for Fd in Fd_list]
DG3 = [TP_Partition[Fd][0] for Fd in Fd_list]
xmin = np.min(DG1)
xmax = np.max(DG1)
xlim = [1.05*xmin, 0.95*xmax]
ymin = np.min(DG3)
ylim = [1.05*ymin, 0.95*xmax]
# Create the plot
plt.scatter(DG1, DG2, color='Blue')
plt.scatter(DG1, DG3, color='Green')
plt.plot(xlim, xlim, color='Black')
plt.xlim(xlim)
plt.ylim(ylim)
# Labels and title
plt.xlabel('Model two-state')
plt.ylabel('DynaMelt Two-state')
plt.title('Free energy (kcal/mol)')
plt.show()

DH1 = [ThermoParam[Fd][0] for Fd in Fd_list]
DH2 = [TP_TwoState[Fd][1] for Fd in Fd_list]
DH3 = [TP_Partition[Fd][1] for Fd in Fd_list]
xmin = np.min(DH1)
xmax = np.max(DH1)
xlim = [1.05*xmin, 0.95*xmax]
ymin = np.min(DH3)
ylim = [1.05*ymin, 0.95*xmax]
# Create the plot
plt.scatter(DH1, DH2, color='Blue', label='Python vs TwoState')
plt.scatter(DH1, DH3, color='Green', label='Python vs Partition')
plt.plot(xlim, xlim, color='Black')
plt.xlim(xlim)
plt.ylim(ylim)
plt.legend()
# Labels and title
plt.xlabel('Model two-state')
plt.ylabel('DynaMelt Two-state')
plt.title('Enthalpy (kcal/mol)')
plt.show()

DS1 = [ThermoParam[Fd][1] for Fd in Fd_list]
DS2 = [TP_TwoState[Fd][2] for Fd in Fd_list]
DS3 = [TP_Partition[Fd][2] for Fd in Fd_list]
xmin = np.min(DS1)
xmax = np.max(DS1)
xlim = [1.05*xmin, 0.95*xmax]
ymin = np.min(DS3)
ylim = [1.05*ymin, 0.95*xmax]
# Create the plot
plt.scatter(DS1, DS2, color='Blue', label='Python vs TwoState')
plt.scatter(DS1, DS3, color='Green', label='Python vs Partition')
plt.plot(xlim, xlim, color='Black')
plt.xlim(xlim)
plt.ylim(ylim)
plt.legend()
# Labels and title
plt.xlabel('Model two-state')
plt.ylabel('DynaMelt Two-state')
plt.title('Entropy (cal/mol/K)')
plt.show()

Tm1 = [ThermoParam[Fd][3] for Fd in Fd_list]
Tm2 = [TP_TwoState[Fd][3] for Fd in Fd_list]
Tm3 = [TP_Partition[Fd][3] for Fd in Fd_list]
xmin = np.min(Tm1)
xmax = np.max(Tm1)
xlim = [0.95*xmin, 1.05*xmax]
ymax = np.max(Tm3)
ylim = [0.95*xmin, 1.05*ymax]
# Create the plot
plt.scatter(Tm1, Tm2, color='Blue', label='Python vs TwoState')
plt.scatter(Tm1, Tm3, color='Green', label='Python vs Partition')
plt.plot(xlim, xlim, color='Black')
plt.xlim(xlim)
plt.ylim(ylim)
plt.legend()
# Labels and title
plt.xlabel('Model two-state')
plt.ylabel('DynaMelt Two-state')
plt.title('Melting temperature (°C)')
plt.show()

print([Tm1[i]-Tm2[i] for i in range(len(Tm1))])
print([DG1[i]-DG2[i] for i in range(len(DG1))])
