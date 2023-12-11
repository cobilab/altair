<div align="center">

![Conda](https://img.shields.io/conda/dn/bioconda/altair-mf)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![Speed](https://img.shields.io/static/v1.svg?label=Testing&message=High-speed%20&color=green)](#)
[![HF](https://img.shields.io/static/v1.svg?label=Testing&message=High-flexibility&color=blue)](#)
[![AFM](https://img.shields.io/static/v1.svg?label=Method&message=alignment-free&color=yellow)](#)

</div>

<p align="center"><img src="imgs/altair.png" alt="AltaiR" width="250" border="0" /></p>
<p align="center">
<b>AltaiR: a C toolkit for alignment-free and spatial-temporal analysis of multi-FASTA data</b>. 
</p>

<p align="justify">
This method provides alignment-free and spatial-temporal analysis of multi-FASTA data through the implementation of a C toolkit highly flexible and with characteristics covering large-scale data, namely extensive collections of genomes/proteomes. This toolkit is ideal for scenarios entangling the presence of multiple sequences from epidemic and pandemic events. AlcoR is implemented in C language using multi-threading to increase the computational speed, is flexible for multiple applications, and does not contain external dependencies. The tool accepts any sequence(s) in (multi-) FASTA format.

The AltaiR toolkit contains one main menu (command: <b>AltaiR</b>) with the six sub menus for computing the features that it provides, namely
<ul>
<li><b>average</b>: moving average filter of a column float CSV file (the column to use is a parameter);</li>
<li><b>filter</b>: filters FASTA reads by characteristics: alphabet, completeness, length, CG quantity, multiple string patterns and pattern absence; </li>
<li><b>frequency</b>: computes the alphabet frequencies for each FASTA read (it enables alphabet filtering);</li>
<li><b>nc</b>: computes the Normalized Compression (NC) for all FASTA reads according to a compression level;</li>
<li><b>ncd</b>: computes the Normalized Compression Distance (NCD) for all FASTA reads according to a reference;</li>
<li><b>raw</b>: computes Relative Absent Words (RAWs) with CG quantity estimation for all RAWs.</li>
</ul>
</p>

## INSTALLATION ##

### Conda
Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html), then run the following:
```bash
conda install -y -c bioconda altair-mf
```

Otherwise, CMake is needed for installation (http://www.cmake.org/). You can download it directly from http://www.cmake.org/cmake/resources/software.html or use an appropriate packet manager. In the following instructions we show the procedure to install, compile and run AltaiR:

<pre>
sudo apt-get install cmake git
git clone https://github.com/cobilab/altair.git
cd altair/src/
cmake .
make
</pre>

### Additional Tools
For certain scripts, the Gto toolkit is required, installable via Conda:
```bash
conda install -c cobilab gto --yes
```
Or manually:
```bash
git clone https://github.com/cobilab/gto.git
cd gto/src/
make
export PATH="$HOME/gto/bin:$PATH"
```

## PARAMETERS

To see the possible options type
<pre>
AltaiR
</pre>
or
<pre>
AltaiR -h
</pre>

If you are not interested in viewing each sub-program option, type 
<pre>
AltaiR average -h
AltaiR filter -h
AltaiR frequency -h
AltaiR nc -h
AltaiR ncd -h
AltaiR raw -h
</pre>

## Reproducing Experiments
Assuming AltaiR is compiled under the `src/` folder, and you are in the `pipeline/` folder.
```bash
cp ../src/AltaiR .
```

### Filtering Sequences
To filter sequences use the following command:
```bash
python3 Histogram.py
bash Filter.sh 29885 29921
```

### Similarity Profiles (NCD)
To simulate and measure similarity profiles:

```bash
bash Simulation.sh
bash Similarity.sh ORIGINAL.fa
bash SimProfile.sh sim-data.csv 2 0 1.2
mv NCDProfilesim-data.csv.pdf NCD_P1.pdf
```

### Phylogenetic Tree Construction
Use the `tree.py` script to construct a phylogenetic tree from NCD values:
```bash
python3 tree.py sim-data.csv -N 50
```

### Complexity Profiles (NC)
Run the following script to generate complexity profiles:
```bash
bash ComplexitySars.sh
python3 CompProfileSars.py comp-data.csv sorted_output.fa 0.961 0.9617
mv NCProfilecomp-data.csv.pdf NC.pdf
```

### Frequency Profiles
Generate frequency profiles using the following commands:
```bash
bash FrequencySars.sh
python3 combine_freq_and_date.py
mv base_frequencies_plot.pdf Freq.pdf
```

### Relative Singularity (RAWs) Profiles
To calculate RAWs profiles:
```bash
bash RawSars.sh
python3 RawSarsProfile.py sorted_output.fa
mv relativeSingularityProfile.pdf RAWProfiles.pdf
```

## Citation

If you use AltaiR in your research, please cite:
- *pending*

## Issues

For any issues, please report at [AltaiR Issues](https://github.com/cobilab/altair/issues).

## License

AltaiR is licensed under GPL v3. For more information, visit [GPL v3 License](http://www.gnu.org/licenses/gpl-3.0.html).
