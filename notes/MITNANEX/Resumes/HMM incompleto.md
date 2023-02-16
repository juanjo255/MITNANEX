Es un modelo estadistico que describe la evolucion de eventos observables que dependen de factores internos.

Se le llama a los eventos observables **simbolos**, mientras que el factor invisible se le conoce como **Estado** 

Los simbolos son representados por **Xn**, mientras que los estados con **Yn**, donde Yn es el evento invisible que le subyace a Xn.

**X** toma un finito numero de posibles valores del conjunto de observaciones; **O** = {O1...ON}.

**Y** toma un valor del conjunto de estados; **S** = {1 ... M}.

**N** y **M** representan el diferente numero de observaciones y estados del modelo, respectivamente.

El HMM esta compuesto por la medicion de **3** probabilidades:
+ Probabilidad de transicion: se mide la probabilidad de transicion de un estado a otro, denotado **t(i,j)** donde **i** es el estado actual.
+ Probabilidad de estado inicial: Es la probabilidad del estado inicial Y1, denotado por  **π(_i_)**
+ Probabilida de emision: La probabilida de que una observacion **n** sea Xn = X estara determinado por el estado Yn, denotado por **_e_(_x_ | _i_)**

Todos las anteriores parametos se resumen en el simbolo  **Θ**. 

## PERFIL DE MODELO OCULTO DE MARKOV

Esta es una variante del HMM usado para representar el perfil de un alinemiento multiple en el que se quieren determinar propiedades estadisticas como patrones y motifs.

Posee una arquitectura de izquierda a derecha sin ciclos.

Usa 3 tipos de estados ocultos: Match states (Mk), Insert states (Ik) y delete states (Dk)