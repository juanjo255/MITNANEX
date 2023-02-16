# [GetOrganelle](https://doi.org/10.1186/s13059-020-02154-5)

This script exploits **Bowtie2**, **BLAST**, and **SPAdes**, as well as the Python libraries **Numpy**, **Scipy**, and **Sympy** as dependencies.

Bases de datos:

**LabelDatabase:** Es usada para identificar los contigs, por ende, esta conformada de regiones conservadas, por ejemplo, genes codificantes para proteina. 

**SeedDatabase**: Es la primera base de datos que se usa para reclutar los primeros reads asociados a organelas como mitocondria, pero despues de dicho reclutamiento no se usa más.
 

#### Step 1: Mapping reads to seed and assembling seed-mapped reads for parameter estimation

+ En este primer paso, se usa bowtie2 para mapear los reads contra una **seedDatabase** la cual contiene una mitocondria de referencia u fragmentos de la organella. Los reads mapeados ahora serán usados para mapear más reads. Estos reads mapeados a su vez son ensamblados en contigs para la estimación de parametros.

#### Step 2. Recruiting more target-associated reads through extending iterations

Aqui primero GetOrganelle usa un algoritmo de pre-agrupamiento para acelerar el proceso de mapeo, este pre-agrupamiento se realiza al comparar reads y juntar aquellos duplicados, pues la logica es que los reads duplicados, como la mitocondria posee mayor cobertura, seran mas probables de ser reads de mitocondria. Al final, tenemos varios grupos de reads duplicados y no duplicados junto con substrings de todos los reads de estos grupos. si dos grupos comparten alguna de estas substring son combinados los grupos, formando de esta forma pseudo-contigs.

Se realizan iteraciones de mapeos con reads semilla que se van actualizando con cada nuevo mapeo y extendiendolos a traves de overlaps usando substrings. Para estas iteraciones los reads semillas son cortados en substring con los que se compara cada read de cada grupo formado en el pre-agrupamiento, si un read de un determinado grupo coincide todo el grupo de reads es tratado como acetable para la mitocondria. el tamaño de la palabra es importante y es determinado por las caracteristicas de los datos, este se ve afectado por la longitud de la lectura, la calidad de la lectura, el número total de lecturas, el porcentaje de lecturas del genoma de los orgánulos, la heterogeneidad de la cobertura de la base de los orgánulos y otros factores.

#### Step 3. Conducting de novo assembly

Reclutados los reads se usa SPAdes para ensamblar el genoma. GetOrganelle se aprovecha del gradiente de Kmers realizado por SPAdes para checkear desde el kmer más pequeño hasta el más largo. Por ejemplo, con trazos repetitivos y una cobertura suficiente es preferible un ensamble de Kmers largo.

#### Step 4. Roughly filtering for target-like contigs

Dado que el genoma de las organelas puede compartir secuencias con otras organelas o con el genoma nuclear, es muy probable que reads que no pertenezcan al genoma se cuelen y se formen contigs ajenos. Para resolver esta problematica GetOrganelle utiliza la **LabelDatabase** con lo que se determina a través de un **BLAST** los genes conservados y su origen (organela o nuclear). Así, por medio de secuencias conservadas en las organelas se elimina aquellos contigs que no presenten hits.

#### Step 5. Identifying target contigs and exporting all configurations

1. Some contigs, including mitochondrial contigs that have short sequence of plastome origin or target-like shallow-depth contaminant contigs, would be labeled incorrectly as target-hit-contigs **(false positive)**. On the other hand, some sequences might be true target contigs but are too short or divergent from sequences in the label database to be labeled as target-hit-contigs **(false negative)**. 

	Para abarcar la problematica anterior, getOrganelle utiliza un sistema de puntajes, por ejemplo, dado unos hits a traves del BLAST a la **LabelDatabase** ya que solo un contig representará un gen, solo el *record* con el mejor **Hit Weight (HW)** value es conservado. El HW se calcula como el producto del tamaño del hit del query y la cobertura del query del contig. Luego, getOrganelle calcula el **Contig Weight (CW)** value sumando los HW de todos sus genes, los contig con mayor puntaje serán más probables de ser contigs relacionados con las organelas; con esto podemos deshacernos del los falsos positivos, pues estos de acuerdo con los autores, suelen venir en queries cortos y de poca profundidad.

2. Si tenemos un contig desconocido solo será clasificado como verdadero, si conecta con la cabeza y cola de contigs verdaderos.




