+ Se realizó la preparacion de libreria usando el kit de ligacion SQK-LSK110, se siguio el siguiente protocolo:
![[genomic-dna-by-ligation-sqk-lsk110_10Nov2020-flongle.pdf]]



+ Los metadatos, que contienen los calculos de Qubit y concentracion de muestra, están en la siguiente tabla de excel: ![[Metadata_nanopore_12_12_2022 1.xlsx]]
+ La celda que se escogio contenia 85 poros.
+ Se optuvo una concentracion final de 4.89 ng/uL para 6 uL de muestra final. Estos se le sumo 4 uL para un total de 10 uL con el fin de obtener suficiente libreria para cargar 2 celdas (5 y 5). Se estimo el tamaño de la libreria en 2kb lo que nos daba un total de 21 fmol. 

## UPDATE
+ La secuenciación se fue a la porra. A pesar de que existian 85 poros solo hubo un promedio de 5 a 10 poros trabajando.
+ Posibles fallas:
	1. ***ADN muy diluido*:** Aunque podria ser refutable, pues al final se encontró una concentración de 4,89 ng/uL lo que daría 29.34 ng. Dado que lo recomendado es entre 3-20 fmol. Siguiendo está tabla brindada por Nanopore:
	![[Screenshot 2022-12-13 at 8.15.16 AM.png|400]]
		El tapeStation nos indico que el tamaño no estaba por debajo de los 1500, por lo que cualquier asumsion de tamaño promedio de fragmento por menor a 8kb, que es lo más seguro, daria una medida correcta de fmoles para la Flongle. Por tanto, tomamos la decision de diluirla hasta los 10 uL con el objetivo de poder tener carga para 2 celdas (5 uL cada una) con lo que quedaría aprox. 2.9 ng/uL, lo que reduciria los fmoles casi a la mitad y seguiría cumpliendo con el rango.
	
	2. ***Realización de la libreria*:** Durante el protocolo se usaron AMPure beads vencidas. Tambien tal vez el paso de los adaptadores falló, pues los poros nos presentaron sintomas de contaminacion por Etanol.
	3. ***Cargada de la celda*:** Esta es la que más me suena, pues durante la cargada de la celda un poco del liquido quedó por fuera del pozo, como una piscinita, por lo que se pudo no llenar completamente con el buffer. Esto no pasó en un videos que he observado.

