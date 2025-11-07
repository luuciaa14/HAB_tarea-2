
# Tarea 2: Network Propagation

## ¿Qué es el *Network Propagation*?

El **Network Propagation**, o propagación en redes, es una técnica computacional ampliamente utilizada en la biología de sistemas. Se utiliza para analizar la estructura y la dinámica de redes complejas, como las interacciones proteína-proteína o las redes génicas. La idea fundamental es la de difundir información a través de la red, de manera que los nodos conectados influyan entre sí. Este proceso permite extender el conocimiento obtenido de un conjunto inicial de genes hacia otros genes conectados, que podrían estar relacionados funcional o biológicamente. 

En el contexto biológico, esta técnica ayuda a identificar genes potencialmente implicados en una enfermedad, proceso celular o ruta metabólica, incluso si no se conocían previamente.
Entre los algoritmos más usados para implementar propagación en red destacan:

* **Random Walk with Restart (RWR)**: simula un "paseante aleatorio" que se mueve por la red y periódicamente regresa a los nodos semilla.
* **Heat Diffusion (HD)**: modela la propagación de información como si fuera calor que se difunde a través de la red.
* **GUILD y DIAMOnD**: frameworks que aplican variantes de estas ideas para priorizar genes candidatos o detectar módulos asociados a enfermedades.

## Contexto de la práctica

En esta práctica se implementará un ejemplo de **propagación en redes utilizando los algoritmos GUILD y DIAMOnD**, dos enfoques complementarios para analizar redes biológicas.

A partir de una red de interacción entre proteínas (obtenidas de **STRING**), se usarán como *genes semilla* los genes **ENO1**, **PGK1** y **HK2**, todos ellos relacionados con la **glucólisis**, una vía metabólica esencial en la producción de energía celular.

Con ello, buscaremos:

* **Priorizar genes candidatos** que podrían estar funcionalmente asociados a los iniciales.
* **Visualizar la difusión de información** en una red biológica real.
* **Explorar diferencias** entre los métodos GUILD (difusión iterativa tipo RWR) y DIAMOnD (expansión modular progresiva).

Los resultados esperados son archivos con los **genes ordenados por relevancia** o **su posición dentro del módulo**, lo que nos permitirá interpretar qué regiones de la red están más relacionadas con las semillas seleccionadas.

## Procesamiento de la red de STRING

Antes de aplicar los algoritmos de propagación, es necesario **preparar la red de interacciones proteína–proteína** para que tenga un formato adecuado.

El archivo original descargado de **STRING** contiene millones de interacciones codificadas mediante identificadores de proteínas, junto con un valor de combined score que indica la confianza de la interacción.

El script `process_STRING.py` realiza el **preprocesamiento y filtrado de esta red** para dejarla lista para el análisis.

### Objetivos del script

* **Filtrar las interacciones** con un nivel de confianza alto (`combined_score ≥ 800`), de forma que solo se mantengan aquellas con evidencia robusta.
* **Convertir los identificadores de proteínas** (ENSP) en **identificadores de genes (ENTREZ)** utilizando la librería `mygene`.
* **Guardar una red simplificada y limpia**, lista para los algoritmos de propagación.

### Pasos del procesamiento

1. **Lectura del archivo de STRING** en fragmentos grandes (*chunks*) para manejar mejor el uso de memoria.
2. **Filtrado de interacciones** por puntuación (`combined_score >= 800`).
3. **Recopilación de los identificadores ENSP únicos** presentes en las interacciones.
4. **Mapeo de ENSP a ENTREZ** mediante consultas a la API de MyGene.
5. **Eliminación de filas sin correspondencia** (aquellas proteínas que no pudieron mapearse).
6. **Exportación final** de un archivo `.txt` o `.tsv` con tres columnas:
```
protein1_entrez   protein2_entrez   combined_score
```

### Ejemplo de uso

El script puede ejecutarse desde la línea de comandos:
```
python scripts/process_STRING.py
```
### Resultado esperado

Un archivo en data/ con nombre similar a:

```
string12_9606_800_entrez.txt
```

que representa la red filtrada y mapeada, lista para ser utilizada por los algoritmos **GUILD** o **DIAMOnD** en el siguiente paso del análisis.

## Carga y representación de la red

En primer lugar, se implementa la **funcionalidad básica para leer y representar la red de interacciones biológicas**.
El objetivo es comprobar que los archivos de red están correctamente formateados y que el script puede construir un **grafo** que posteriormente servirá como base para los algoritmos de propagación.

El código utiliza las librerías **pandas** y **networkx**, dos herramientas muy comunes en bioinformática para manejar datos tabulares y redes biológicas respectivamente.
A partir de un archivo que contiene pares de genes o proteínas conectados (por ejemplo, una interacción en la red STRING), el script crea un **grafo no dirigido**, donde:

* Los **nodos** representan genes o proteínas.
* Las **aristas** representan las interacciones entre ellos.

Una vez cargada la red, el script imprime por pantalla el número total de nodos y aristas, lo que permite validar que los datos se han leído correctamente antes de continuar con los pasos siguientes (lectura de semillas y propagación).

De esta forma, asentamos las bases del código, ya que **toda la propagación de información en versiones posteriores dependerá de este grafo inicial**.
En bioinformática, este tipo de representación es fundamental para estudiar cómo los genes o proteínas se relacionan entre sí y para poder difundir señales o puntuaciones a través de la red.

## Lectura de genes semilla y compatibilidad con formatos de red

En la segunda versión del script se añadió la funcionalidad para **leer los genes semilla** y comprobar su presencia dentro de la red cargada.
Los genes semilla son el **punto de partida** del proceso de propagación: a partir de ellos se difundirá la información para identificar otros genes relacionados funcionalmente.

El archivo `data/genes_seed.txt` contiene los nombres de los genes semilla (uno por línea):

```
ENO1
PGK1
HK2
```

Estos genes pertenecen a la vía glucolítica y se utilizan como referencia para iniciar la propagación dentro de la red biológica.

El script incorpora ahora la función `load_seeds(path)`, que lee el archivo y devuelve una lista de nombres.
Posteriormente, compara esos genes con los nodos de la red para verificar si están presentes, mostrando un resumen con:

* El número total de semillas leídas,
* Cuántas están realmente en la red,
* Cuáles no se han encontrado.

Esto permite **comprobar la coherencia entre los identificadores** de la red y los de las semillas antes de aplicar cualquier algoritmo.

## Compatibilidad con distintos tipos de red

En este punto también se adaptó la función `load_network(path)` para poder reconocer automáticamente **distintos formatos de red**:

* Redes con identificadores numéricos (ENTREZ), como `network_guild.txt` y `network_diamond.txt`.
* Redes con identificadores de genes en formato HUGO, como `string_network_filtered_hugo-400.tsv`.

El script detecta si el archivo contiene las columnas `protein1_hugo` y `protein2_hugo` (propias del formato STRING) y ajusta la lectura en consecuencia, evitando que la primera línea de cabecera se interprete como un nodo más.
De esta forma, el código es más flexible y puede usarse indistintamente con diferentes redes biológicas.

## Propagación tipo GUILD (Random Walk with Restart)

En esta versión se implementa el primer algoritmo de propagación: **Random Walk with Restart (RWR)**, una de las bases conceptuales de los métodos de priorización usados en GUILD.

El enfoque se basa en un modelo de **paseo aleatorio sobre la red**, donde un “walker” (o señal) se mueve de un nodo a sus vecinos con cierta probabilidad. En cada paso, existe una probabilidad `r` (llamada restart) de volver a los **nodos semilla**.
De este modo, los nodos cercanos a las semillas en la topología de la red acumulan más señal, obteniendo una puntuación más alta.

### Funcionamiento general del algoritmo

1. Se asigna una puntuación inicial de 1/N a cada semilla (donde N es el número de semillas).
2. En cada iteración:
    * La señal de cada nodo se reparte equitativamente entre sus vecinos.
    * Una fracción r del total se “reinicia” volviendo a las semillas.
3. El proceso se repite hasta que las puntuaciones se estabilizan (convergencia).

### Parámetros principales

* `restart`: controla cuánto peso se devuelve a las semillas en cada paso (por defecto 0.5).
* `max_iter`: número máximo de iteraciones permitidas.
* `tol`: tolerancia usada para detener el proceso si las puntuaciones no cambian significativamente.

### Resultados esperados

El script devuelve una puntuación (score) para cada gen o proteína de la red, que refleja su **proximidad funcional a las semillas**.
Tras la ejecución, se muestran los **5 genes con mayor puntuación**, representando aquellos más conectados o relevantes según el algoritmo.

## Expansión de módulo estilo DIAMOnD

En esta versión se incorpora una segunda forma de analizar la red, inspirada en el algoritmo **DIAMOnD (Disease Module Detection)**. A diferencia del enfoque basado en difusión (RWR), este método trabaja de forma **constructiva**: parte del conjunto de genes semilla y va añadiendo, uno a uno, los genes que están más conectados con ese conjunto.

La intuición es la siguiente: si varios genes de una misma enfermedad o proceso biológico están cerca en la red, entonces los genes que se conectan fuertemente con ellos también son buenos candidatos a pertenecer al mismo módulo funcional.

### Cómo funciona la versión simplificada

1. Se inicializa el **módulo** con los genes semilla presentes en la red.
2. Para cada nodo que aún no está en el módulo, se cuenta cuántos vecinos tiene dentro del módulo actual.
3. Se añade al módulo el nodo con **mayor número de vecinos dentro del módulo**.
4. Se repite el proceso hasta alcanzar un tamaño máximo (`k`) o hasta que ya no haya nodos conectados.

El resultado es un **orden de entrada** de los genes al módulo, que indica qué genes parecen más cercanos (topológicamente) al conjunto inicial.

Este comportamiento es similar al de DIAMOnD original, pero sin la parte estadística completa (p-valores y modelo nulo). Para la práctica es suficiente porque muestra la idea de “crecer” un módulo de enfermedad a partir de unas semillas.

### Parámetros

- `k`: número máximo de genes que queremos en el módulo. En la práctica se pueden usar valores entre 30 y 200 según el tamaño de la red.
- semillas: se usan las mismas que en la propagación RWR (`ENO1`, `PGK1`, `HK2`), por lo que ambas partes del script son comparables.

### Salida

Al ejecutar el script, además de la salida del RWR, se mostrará una lista ordenada de los primeros genes añadidos al módulo, por ejemplo:

```text
[OK] Expansión estilo DIAMOnD completada. Primeros genes del módulo:
  01. ENO1
  02. PGK1
  03. HK2
  04. GAPDH
  05. ALDOA
  06. TPI1
  ...
