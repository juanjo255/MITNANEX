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
+ [Algorithm overview](#algorithm_overview)
+ [Contributing](../CONTRIBUTING.md)

## 🧐 About <a name = "about"></a>
MITNANEX's main purpose is to extract mitocondrial Nanopore reads **_De novo_** from the **WGS**, with no need for seeds or reference sequences. However, it will also returned a draft assembly of the mitogenome using [Flye](https://github.com/fenderglass/Flye.git).

## 🏁 Getting Started <a name = "getting_started"></a>

### Installing
* First, you need to clone this repository and add to PATH:
```
  git clone https://github.com/juanjo255/MITNANEX.git; cd MITNANEX; export PATH=$(pwd):$PATH
```

#### Conda/mamba

The best way to install MITNANEX's dependencies is throught a beautiful conda/mamba enviroment, additionally you need to have Rust installed (https://www.rust-lang.org/tools/install).

* For Mac M1 using mamba (you can change it for conda):
  ```
  CONDA_SUBDIR=osx-64; mamba create -n mitnanex -c conda-forge -c bioconda seqkit seqtk fpa minimap2 miniasm flye gfastats samtools
  mamba activate mitnanex
  pip install pandas maturin biopython scikit-learn utils-mitnanex
  ```
* For Linux
  ```
  conda create -n mitnanex -c conda-forge -c bioconda Seqkit Seqtk fpa Minimap2 Miniasm Flye Gfastats Samtools
  conda activate mitnanex
  pip install pandas maturin biopython scikit-learn utils-mitnanex
  ```
   
### Dependencies
MITNANEX needs the following tools:
1. Seqkit
2. Seqtk
3. fpa
4. Minimap2
5. Miniasm
6. Flye
7. Pandas
8. Gfastats
9. Samtools
10. Filtlong
11. Maturin
12. Biopython
13. scikit-learn
14. utils-mitnanex 

**Notes:** 

+ This has only been tested on MacOS M1 using a x86 env architecture.
  
+ ```setup.sh``` will create a **mamba** enviroment with all the dependencies in the  ```.yml``` file.

## 🎈 Usage <a name="usage"></a>

* Quick start
  ```
  ./mitnanex_cli.sh -i path/to/fastQ  -p 15000 -m 1000 -t 8 -s 0.6 -g GenomeSize(g|m|k) -w path/to/output
  ```
**Notes:** 
+ It only receives fastQ files.

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
        -d        Different output directory. Create a different output directory every run (it uses the date and time). [False]
        -s        Mapping identity. Minimun identity between two reads to be store in the same cluster.[0.6]
        -q        Min mapping quality (>=). This is for samtools. [-1].
        -f        Flye mode. [--nano-hq]
        -g        GenomeSize. This is your best estimation of the mitogenome for read correction with Canu. [required]
        *         Help.
  
  ```

  ## Algorithm overview <a name = "algorithm_overview"></a>
  ### **How does MITNANEX work?**
  + MITNANEX is a pipeline that depends on other open source tools (see [dependencies](#getting_started)).
  + Through this I will show the results that belong to the assemble of Talaromyces santanderensis mitogenome using MITNANEX from a Nanopore run performed at EAFIT university.
  + First, it will use seqkit and seqkt to subsample the reads, after that  MITNANEX starts with minimap2 finding overlaps between reads. MITNANEX will group reads that have at least certain level of identity (tweakable parameter), each read will be counted for the "coverage" of the group and each cluster will be represented only by its largest read.
  + Once all reads are grouped, MITNANEX will only keep at least 3 groups with the highest coverage (tweakable parameter). Given the short length of the mitchondrial genome and its high coverage during WGS, we expect to have most of it in these clusters.
  + <p align="center"><img src="images/Cluster_filter_by_cov.png" alt="Cluster_filter_by_cov"></p>
  + Now with the selected clusters MITNANEX will use the representative read of each cluster and get its trinucleotidic composition (codon) which will be reduce is normalized by the read length, and reduce its dimensionality to 2 with a PCA such as the classic strategy during metagenomic binning. Here, given the difference between mitochondrial and the nuclear genome, we expect the mitochondrial reads to have an oligocomposition different enough to be separated from the nuclear. The known weakness of Kmean for outliers made the selection of this clustering algortihm attractive. Thus, using the clustering algorithm Kmeans, with a k set to 2, is selected the cluster with the highest coverage. Below the cluster in yellow was selected.
  + <p align="center"><img src="images/Kmeans_on_pca.png" alt="Kmeans_on_pca"></p>
  + With the reads collected from the selected clusters, miniasm will assemble unitigs, where we expect to assemble most of the mitogenome (or even longer given the problems that miniasm has). Miniasm is useful in this steps for 2 main reasons:
    1. It can work with low coverage.
    2. It's extremely fast and the unitigs produced are enough for the next step.
  + The unitigs are used to collect more reads from the total of reads to perform a final assembly with Flye.
  + Flye is almost the assembler par excellence for Nanopore reads and it's among the best at circularizing genomes. An important characteristic for mitochondrial genomes. 
     
 






