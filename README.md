<div align="center">

![Conda](https://img.shields.io/conda/dn/bioconda/altair)
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

## CITATION ##

On using this software/method please cite:

* pending

## ISSUES ##

For any issue let us know at [issues link](https://github.com/cobilab/altair/issues).

## LICENSE ##

GPL v3.

For more information:
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>

