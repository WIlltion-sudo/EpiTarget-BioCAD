
<div align="center">

# 🧬 EpiTarget-BioCAD
### *In-Silico Restriction Digest & Thermodynamic Melting Temperature (T_m) Engine*

[![Python Version](https://img.shields.io/badge/Python-3.8%2B-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](LICENSE)
[![Dependencies](https://img.shields.io/badge/Dependencies-Zero--Biopython-FF6F61?style=for-the-badge)](https://matplotlib.org/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey?style=for-the-badge)](#)
[![Field](https://img.shields.io/badge/Field-Bioinformatics-purple?style=for-the-badge&logo=dna)](https://en.wikipedia.org/wiki/Bioinformatics)

**An advanced biophysical computation engine for precise thermodynamic duplex stability profiling, non-overlapping/overlapping restriction site mapping, and virtual agarose gel electrophoresis simulation.**

<p align="center">
  <a href="#-key-features">Key Features</a> •
  <a href="#-biophysical-mechanics">Biophysical Mechanics</a> •
  <a href="#-installation--setup">Installation</a> •
  <a href="#-usage-and-execution">Quick Start</a> •
  <a href="#-license">License</a>
</p>


</div>

## 🎯 Overview

Moving beyond static base-frequency analytics, **EpiTarget-BioCAD** bridges raw sequence data with biophysical thermodynamics[cite: 1, 2]. By executing unified nearest-neighbor thermodynamic matrices alongside lookahead endonuclease scanning, this engine enables researchers and developers to simulate enzymatic cleavages and evaluate oligonucleotide annealing behavior *in-silico*[cite: 1].

```text
       5' -- [G A A T T C] -- 3'           --->  Enzymatic Cleavage Engine
       3' -- [C T T A A G] -- 5'                 |-- Cut Coordinates Map
                                                 |-- Fragment Migration Vectors
                                                 `-- Logarithmic Gel Plotting

```

---

## ⚡ Key Features

| Feature | Description | Technical Implementation |
| --- | --- | --- |
| **Nearest-Neighbor $T_m$ Engine** | Computes precise oligonucleotide melting temperatures.

 | Implements **SantaLucia (1998)** unified $\Delta H^\circ, \Delta S^\circ$ parameters with salt $[Na^+]$ correction.

 |
| **Endonuclease Cut Scanner** | Scans DNA for standard restriction recognition sites (*EcoRI*, *BamHI*, *HindIII*, *NotI*, *TaqI*, *HaeIII*).

 | Zero-width lookahead regex assertion `(?=MOTIF)` to capture overlapping motifs.

 |
| **In-Silico Gel Plotter** | Simulates agarose gel electrophoresis.

 | Dark-mode matplotlib rendering with logarithmic migration distance relative to a molecular ladder.

 |
| **Self-Complementary Logic** | Detects self-palindromic oligonucleotides.

 | Automatically applies symmetry factor adjustments ($x=1$ vs $x=4$) to thermodynamic concentration terms.

 |

---

## 🔬 Biophysical Mechanics

### 1. Thermodynamic Nearest-Neighbor Model (SantaLucia 1998)

Standard $2(A+T) + 4(G+C)$ estimations fail for non-canonical sequence lengths and salt variations. This module evaluates adjacent dinucleotide steps ($\Delta H^\circ_{\text{NN}}, \Delta S^\circ_{\text{NN}}$) alongside terminal initiation penalties:

$$T_m = \frac{\Delta H^\circ_{\text{total}} \times 1000}{\Delta S^\circ_{\text{total}} + R \ln(C_T / x)} - 273.15 + 16.6 \log_{10}([Na^+])$$

* **$\Delta H^\circ, \Delta S^\circ$**: Enthalpy ($\text{kcal/mol}$) and Entropy ($\text{cal/K}\cdot\text{mol}$).


* **$R$**: Universal Gas Constant ($1.9872 \text{ cal/K}\cdot\text{mol}$).


* **$C_T$**: Oligonucleotide concentration ($500 \text{ nM}$ standard).


* **$x$**: Symmetry factor ($x=4$ for non-self-complementary, $x=1$ for self-complementary).


* **$[Na^+]$**: Monovalent sodium ion concentration ($50 \text{ mM}$ default assay buffer).



### 2. Enzymatic Cut & Logarithmic Gel Separation

* **Cut Mapping:** Splitting a sequence of length $L$ at cut positions $[p_1, p_2, \dots, p_n]$ yields $n+1$ fragments.


* **Migration Velocity:** Gel migration distance correlates inversely with the logarithmic base-pair lengt.



---

## ⚙️ Installation & Setup

### 1. Clone the Repository

```bash
git clone https://github.com/Willtion-sudo/EpiTarget-BioCAD.git
cd EpiTarget-BioCAD

```

### 2. Set Up a Virtual Environment

* **On Windows (Command Prompt):**
```cmd
python -m venv test_env
test_env\Scripts\activate

```


* **On Linux / macOS:**
```bash
python3 -m venv test_env
source test_env/bin/activate

```



### 3. Install Requirements

```bash
pip install -r requirements.txt

```

---

## 🚀 Usage and Execution

Run the core engine script directly from your terminal:

```bash
python src/digest_thermo.py

```

### Terminal Walkthrough

```text
==========================================================
 EpiTarget-BioCAD: Digest & Thermo Engine 
==========================================================
Enter DNA sequence (Press Enter to use default demo): ATGCGATCGATCGATCGATCGATCGATCGATCGATCGA
Analyzing Sequence (Length: 38 bp)...

--- THERMODYNAMIC METRICS ---
Nearest-Neighbor dH° : -289.4 kcal/mol
Nearest-Neighbor dS° : -784.2 cal/(K*mol)
Calculated Tm (50mM Na+) : 71.45 °C

Available Enzymes: EcoRI, BamHI, HindIII, NotI, TaqI, HaeIII
Choose enzyme for digest (default: EcoRI): TaqI

--- RESTRICTION DIGEST RESULTS (TaqI) ---
Recognition Motif : TCGA
Cut Indices       : [5, 9, 13, 17, 21, 25, 29, 33]
Total Fragments   : 9
Fragment Sizes    : [5, 4, 4, 4, 4, 4, 4, 4, 5] bp
```
## 📄 License

Distributed under the **MIT License**. See [`LICENSE`](https://www.google.com/search?q=MIT+LICENSE) for complete details.
