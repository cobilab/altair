<div align="center">

![Conda](https://img.shields.io/conda/dn/bioconda/altair-mf)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![Speed](https://img.shields.io/static/v1.svg?label=Testing&message=High-speed%20&color=green)](#)
[![HF](https://img.shields.io/static/v1.svg?label=Testing&message=High-flexibility&color=blue)](#)
[![AFM](https://img.shields.io/static/v1.svg?label=Method&message=alignment-free&color=yellow)](#)

</div>

<p align="center"><img src="imgs/altair.png" alt="AltaiR" width="250" border="0" /></p>

<p align="center">
<b>AltaiR</b><br/>
Alignment-free and temporal analysis of multi-FASTA data (C toolkit)
</p>

* * *

## ✨ What is AltaiR?

AltaiR is a fast, alignment-free toolkit for **temporal analysis and characterization of multi-FASTA datasets**, targeting large-scale collections such as **genomes** and **proteomes**.

It is particularly useful for scenarios with many sequences collected over time (e.g., epidemic/pandemic datasets), where alignment-based workflows can be slow, brittle, or unnecessary for the desired analyses. AltaiR is implemented in **multi-threaded C**, is **highly flexible**, and is designed to run **without external dependencies** (core toolkit). It accepts any sequence(s) in **(multi-)FASTA** format.

### ✅ Highlights

- ⚡ High speed (multi-threaded C implementation)
- 🧩 High flexibility (multiple independent analysis modules)
- 🧬 Alignment-free methods (compression-based and k-mer/word-based analyses)
- 📦 No external dependencies for the core toolkit
- 🗂️ Works with any (multi-)FASTA input

* * *

## Contents

- Commands
- ⚙️ Installation
- Quickstart
- Help and parameters
- Reproducing experiments (pipeline)
- Citation
- Issues
- License

* * *

## Commands

AltaiR provides a single entry point (`AltaiR`) with **six subcommands**:

- `average` — moving average filter for a float column in a CSV file (column index is a parameter)
- `filter` — filter FASTA records by alphabet, completeness, length, CG content, presence/absence of string patterns
- `frequency` — compute alphabet frequencies per FASTA record (optionally with alphabet filtering)
- `nc` — compute **Normalized Compression (NC)** per FASTA record (configurable compression level)
- `ncd` — compute **Normalized Compression Distance (NCD)** for each record relative to a reference
- `raw` — compute **Relative Absent Words (RAWs)** with CG% estimation per RAW

* * *

## ⚙️ Installation

### Option A — Conda (recommended)

Install Miniconda (or Mambaforge), then create an environment and install from Bioconda:

```bash
mamba create -n altair -c conda-forge -c bioconda altair-mf
conda activate altair
AltaiR -h
````

To install into an existing environment:

```bash
conda install -y -c bioconda altair-mf
```

### Option B — Build from source (CMake)

Requirements: `cmake`, `git`, and a C compiler toolchain.

```bash
sudo apt-get install -y cmake git build-essential
git clone https://github.com/cobilab/altair.git
cd altair
cmake -S src -B build
cmake --build build -j
./build/AltaiR -h
```

> If you prefer the in-tree build used by some minimal setups:
>
> ```bash
> cd altair/src
> cmake .
> make
> ```

### Optional — Additional tools for pipeline scripts (Gto)

Some scripts in `pipeline/` require the **Gto** toolkit.

Conda:

```bash
conda install -c cobilab gto --yes
```

Manual:

```bash
git clone https://github.com/cobilab/gto.git
cd gto/src/
make
export PATH="$HOME/gto/bin:$PATH"
```

---

## Quickstart

1. Show top-level help:

```bash
AltaiR -h
```

2. Run a subcommand (see each module’s `-h` for required parameters):

```bash
AltaiR filter -h
AltaiR frequency -h
AltaiR nc -h
AltaiR ncd -h
AltaiR raw -h
AltaiR average -h
```

---

## Help and parameters

Top-level help:

```bash
AltaiR
# or
AltaiR -h
```

Per-subcommand help:

```bash
AltaiR average -h
AltaiR filter -h
AltaiR frequency -h
AltaiR nc -h
AltaiR ncd -h
AltaiR raw -h
```

---

## Reproducing experiments (pipeline)

Assuming AltaiR was compiled and you are working under `pipeline/`:

```bash
cp ../src/AltaiR .
```

> Some steps require `python3`, `bash`, and (optionally) `gto` (see “Additional tools”).

### Filtering sequences

```bash
python3 Histogram.py
bash Filter.sh 29885 29921
```

### Similarity profiles (NCD)

```bash
bash Simulation.sh
bash Similarity.sh ORIGINAL.fa
bash SimProfile.sh sim-data.csv 2 0 1.2
mv NCDProfilesim-data.csv.pdf NCD_P1.pdf
```

### Phylogenetic tree construction

```bash
python3 tree.py sim-data.csv -N 50
```

### Complexity profiles (NC)

```bash
bash ComplexitySars.sh
python3 CompProfileSars.py comp-data.csv sorted_output.fa 0.961 0.9617
mv NCProfilecomp-data.csv.pdf NC.pdf
```

### Frequency profiles

```bash
bash FrequencySars.sh
python3 combine_freq_and_date.py
mv base_frequencies_plot.pdf Freq.pdf
```

### Relative singularity (RAWs) profiles

```bash
bash RawSars.sh
python3 RawSarsProfile.py sorted_output.fa
mv relativeSingularityProfile.pdf RAWProfiles.pdf
```

---

## Citation

If you use AltaiR in your research, please cite:

Silva, Jorge M., Armando J. Pinho, and Diogo Pratas. **“AltaiR: a C toolkit for alignment-free and temporal analysis of multi-FASTA data.”**
*GigaScience* 13 (2024): giae086.

[https://doi.org/10.1093/gigascience/giae086](https://doi.org/10.1093/gigascience/giae086)

---

## Issues

Please report bugs and feature requests in the repository issue tracker:

* `https://github.com/cobilab/altair/issues`

---

## License

AltaiR is licensed under **GNU GPL v3**. See [LICENSE](LICENSE).
More information: `http://www.gnu.org/licenses/gpl-3.0.html`
