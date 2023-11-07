<p align="center"><img src="images/MITNANEX.png" alt="MITNANEX"></p>

<h1 align="center">MITochondrial NANopore reads EXtractor</h3>

<div align="center">

  [![Status](https://img.shields.io/badge/status-active-success.svg)]() 
  [![GitHub Issues](https://img.shields.io/github/issues/kylelobo/The-Documentation-Compendium.svg)](https://github.com/juanjo255/MITNANEX/issues)
  [![GitHub Pull Requests](https://img.shields.io/github/issues-pr/kylelobo/The-Documentation-Compendium.svg)](https://github.com/juanjo255/MITNANEX/pulls)
  [![License](https://img.shields.io/badge/license-MIT-blue.svg)](/LICENSE)

</div>

---

## Table of Contents
+ [About](#about)
+ [Getting Started](#getting_started)
+ [Usage](#usage)
+ [Contributing](../CONTRIBUTING.md)

## 🧐 About <a name = "about"></a>
This tool's main purpose is to extract mitocondrial nanopore reads from the WGS _De novo_. However, it will also assembly a draft mitogenome using [Flye](https://github.com/fenderglass/Flye.git).

## 🏁 Getting Started <a name = "getting_started"></a>

### Prerequisites
The tool has a module written in Rust. I have to find out if it is necessary to download rust and/or run masturin

### Installing

#### Conda/mamba

The best way to use the program is throught a beautiful conda enviroment.

**Note** This has only been tested on MacOS M1

```
git clone https://github.com/juanjo255/MITNANEX.git
cd MITNANEX
mamba env create --name mitnanex --file environment_mac.yml
mamba activate mitnanex
pip install -r requirements_macOS.txt
```

## 🎈 Usage <a name="usage"></a>

* Quick start
  ```
  ./mitnanex_cli.sh -i path/to/fastQ  -p 15000 -m 1000 -t 8 -s 0.7 -g 22k -w path/to/output
  ```
* For help message
  ```
  ./mitnanex_cli.sh -h
  ```
  ```
    Options:
        -i        Input file. [required]
        -t        Threads. [4].
        -p        Proportion. For sampling. It can be a proportion or a number of reads (0.3|10000). [0.3].
        -m        Min-len. Filter reads by minimun length. Read seqkit seq documentation. [-1].
        -M        Max-len. Filter reads by maximun length. Read seqkit seq documentation. [-1].
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [./mitnanex_results].
        -r        Prefix name add to every produced file. [input file name].
        -c        Coverage. Minimum coverage per cluster accepted. [-1].
        -d        Different output directory. Create a different output directory every run (it uses the date and time).
        -s        Mapping identity. Minimun identity between two reads to be store in the same cluster.[0.6]
        -q        Min mapping quality (>=). This is for samtools. [-1].
        -f        Flye mode. [--nano-hq]
        -g        GenomeSize. This is your best estimation of the mitogenome for read correction with Canu. [required]
        *         Help.
  
  ```
