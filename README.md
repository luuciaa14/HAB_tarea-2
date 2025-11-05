
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