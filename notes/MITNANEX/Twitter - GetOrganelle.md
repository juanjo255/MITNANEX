
### GetOrganelle

In my journey learning about mitochondrial genome assembly, I have discovered a few tools. Today let's talk about the awesome mechanism of GetOrganelle. #Bioinformatics #Science


Thread 🧵

---

Let's get into matter. GetOrganelle is a pipeline that recruits organelles reads from WGS. It has two main Python scripts: "get_organelle_from_reads.py" and "get_organelle_from_assembly.py". This scripts exploit tools such as 𝗕𝗼𝘄𝘁𝗶𝗲2, 𝗕𝗟𝗔𝗦𝗧 and 𝗦𝗣𝗔𝗱𝗲𝘀.

---
Now I am going to summarize the step-by-step of how GetOrganelle works. The caveat would be that this is based on what the authors wrote in the paper and my understanding about it.

---
Step 1.  

In this first step GetOrganelle uses Bowtie2 to map WGS reads to a 𝙎𝙚𝙚𝙙𝘿𝙖𝙩𝙖𝙗𝙖𝙨𝙚 which is a database with complete reference organelle genomes or organelle fragments. In this way, it recruits the first target-reads which will be used as new seed.

---
Using the initial mapped-reads, GetOrganelle recruits more target-reads and replace the seed with new overlapped reads and so on until there aren't more new reads.

---
Step 2.

The recruited reads are used to assemble the genome using SPAdes. GetOrganelle takes advantage of the gradient of k-mers performed by SPAdes to assess the assembly from the smallest to the largest k-mer. Since, for example, long k-mers are preferable for repetitive genome

----
Step 3.

At this point, in a perfect world this could be almost the end of the story, we would have an assembly graph with organelle-associated contigs and we could spent the rest of the day watching Netflix or doing something else

---

However, the reality is that organelles and nuclear genome can share sequences between them, therefore non-target reads can leaked what can gives us some troubles in the assembly, since the assembly graph can have non-target-associated contigs

---
In order to filter correct contigs, GetOrganelle uses BLAST to map contigs with a 𝐋𝐚𝐛𝐞𝐥𝐃𝐚𝐭𝐚𝐛𝐚𝐬𝐞, which is a database with coding regions of the organelle genome. Thus, by recording gene identities and type of organelle, we can leave only true target contigs

---

But that is not enough. Identifying true contigs based only on a BLAST hit is risky, so GetOrganelle groups the resulting contigs with the BLAST information and coverage values of them.

---
 I won't go in depth with the concepts of this last step, since are a little more complex and would do this tweet longer than I want, but if you are interested in knowing more about it, I recommend you to read the paper.

----

Step 4.

Finally, using contig multiplicity, GetOrganelle searches for any possible path between contigs. If it fails to export a circular graph assembly (it assumes circular genome) for some reason (e. g. Insufficient target assembly graph), then it will export the target-contigs.

---

The End 🙇‍♂️.

The simplicity and completeness of this tool makes it, for me, one of the best tools to assemble the mitochondrial genome from DNA, although, as always, there is room for improvement...


### **Resultado**

<iframe border=0 frameborder=0 height=500 width=500 src="https://twitframe.com/show?url=https://twitter.com/Juanpicon255/status/1576240896244260865"> </iframe>