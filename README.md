<p align="center"><img src="images/MITNANEX.png" alt="MITNANEX"></p>

<h1 align="center">MITochondrial NANopore reads EXtratore</h3>

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
This is intended to be a extractor for mitocondrial nanopore reads from the whole genome sequencing, in order to be able to assembly mitogenome or do any other proccess where researchers think mitocrondrial reads could be useful.

## 🏁 Getting Started <a name = "getting_started"></a>

### Prerequisites
The tool has a module written in Rust. I have to find out if it is necessary to download rust and/or run masturin
```
Give examples
```

### Installing

#### Conda/mamba

The best way to use the program is throught a beautiful conda enviroment.

**Note** this has only been tested on MacOS M1

```
git clone https://github.com/juanjo255/MITNANEX.git && cd MITNANEX && mamba env create --name mitnanex --file environment.yml
```

## 🎈 Usage <a name="usage"></a>

* Quick start
  ```
  ./mitnanex_cli.sh -i path/to/readsFile -p 0.4 -m 300 -t 4 -s 0.7 -w path/to/output
  ```
* For help message
  ```
  ./mitnanex_cli.sh -h
  ```
  ```
  MITNANEX
      Version: 1.0
      https://github.com/juanjo255/MITNANEX_PROJECT.git
  
      Usage: mitnanex.sh [options] FASTQ
  
      Options:
          -i        Input file.
          -t        Threads. [4].
          -p        Proportion. For sampling with seqkit. Read seqkit sample documentation. [0.4].
          -m        Min-len. Filter reads by minimun length. Read seqkit seq documentation. [-1].
          -M        Max-len. Filter reads by maximun length. Read seqkit seq documentation. [-1].
          -w        Working directory. Path to create the folder which will contain all mitnanex information. [./mitnanex_results].
          -r        Prefix name add to every produced file. [input file name].
          -c        Coverage. Minimum coverage per cluster accepted. [-1].
          -d        Different output directory. Create a different output directory every run (it uses the date and time).
          -s        Mapping identity. Minimun identity between two reads to be store in the same cluster.[0.6]
          -q        Min mapping quality (>=). This is for samtools. [-1].
          *         Help.
  
  ```
