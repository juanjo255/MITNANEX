### MEANGS

1. Continuing my journey through mitochondria assembly. Today I want to talk about a tool released in 2021 that promises to be a different seed-free *de novo* assembler of mitogenome from animal whole genome NGS data; lets talk about **MEANGS** #Bioinformatics #Science 
	
	Thread 🧵
---
2.  Warning ⚠️: All this tweet is based on what the authors wrote in the paper and my understanding about it. In case you want to read the original publication and deepen the tool, here is the link: https://doi.org/10.1093/bib/bbab538
 ----
3. Why is MEANGS different from other seed-free assembler like Norgal? Well, the main disadvantage that have seed-free is that since these tools exploit of the high depth of mitochondrial reads, they despise reads of normal/low depth, resulting in incomplete assemblies. (GIF)
----
4. Then, to solve this problem and achieve an improvement in integrity and accuracy, MEANGS authors though in a solution as simple as simply not filtering low sequencing reads 🥴. But then,  How do they filter mitochondrial reads from other other reads? 
	
	Let's see the step-by-step. (GIF)
----
5.  **Seed generation** 🫘 

	MEANGS uses profile HMMs to select mitochondrial. MEANGS takes from MitoZ a profile Hidden Markov Models (HMM) database that contains mitogenome sequence conservation modelling information from different taxonomic clades in Metazoans adopted.
----
6. Thus, MEANGS uses HMMER to filter candidate mitochondrial reads using a confidence score S. Reads with a high score were selected and, using SSAKE, assembled to generate candidate mitochondrial Protein-Coding Genes (PCG) fragments which will be used as sequence seed.
----
7. If you want all can end here. MEANGS offers 2 modes: default mode and deepin mode. The default mode is limited to assemble only PCG fragments, whereas deepin mode return the complete mitogenome. (GIF)
----
8.  ***De novo* assembly** 🧬
	If deepin mode is selected, MEANGS assemble mitogenome using SSAKE. NGS data generally contain a high depth of mitochondrial sequences, thus the robust seed-and-extend algorithm can extend seed sequences into a complete mitogenome.
----
9.  Finally, MEANGS is not perfect, but it is a good start, as it has demostrated been as good or better than other tools. I take my leave with the histogram below that shows a benchmarking between MEANGS and other open-source software tools using 16 different mitogenomes.
![[Screenshot 2022-12-08 at 8.19.01 PM.png]]
### Resultado
