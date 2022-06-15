# Nucleos

## Separar en canales una imagen .czi
En esta sección se describen los pasos a seguir para separa una imagen .czi en imágenes de los canales de color correspondientes en formato .tif.
El programa hace la aisgnación de canales de forma atomática a partir de la lectura de los metadatos.
Las imágenes resultantes de la separación tendrán el mismo nombre que su imagen original, con la excepción de que al inal se les agregará el nombre del marcador biológica que aparece en los metadatos.
1. Ir a la carpeta __Color_channels__
2. Colocar dentro de la carpeta __Images__ la imagen en .czi de interés.
3. Ir al archivo _SepararCanalesColor.m_ dentro de la carpeta __Color_channels__ y escribir el nombre de archivo de la imagen (sin la extensión .czi).
4. Este código mostrará un mosaico de los canales ya separados.
Las imágenes de canales de color se guardarán la carpeta __cell_images__ con la extensión .tif

## Procesar imágenes de los canales
Esta sección es la más compleja y poderosa de todas, ya que uno le indica al programa la nomenclatura base de las imágenes de los canales de color y este se encarga de analizar todas automáticamente.
* Realiza detección y segmentación de núcleos a partir de la imagen de marcaje nuclear.
* Conteo de nćleos con combinatrias de dupletes de marcadores.
* Realiza segmentación de las regiones ventricular, ventriculo lateral y estriado.
* Enumera y obtiene características físicas de c/u de los núcleos a partir de la imagen de marcaje nuclear.
* Se obtienen diversas métricas, conteo y proporciones de interés.
A continuación se enlistan los pasos para utilizar este programa.
1. Ir a la carpeta __scripts__
2. Accedar al archivo _Celulas_Universal.m_
3. Indicar el nombre de las imágenes de los canales que se tiene como objetivo analizar (si la extensión .tif).
4. Correr el programa.
5. En caso de ser necesario, uno puede ajustar el valor de "area_celula" para mejorar resultados.

## Programa de transparencia
Si uno lo desea, se pueden superponer 2 iágenes y ajustar la transparencia lentamente para verificar diferencias entre sí.

Esto es útil cuando se quiere comparar por ejemplo:
* Una imagen original de marcaje nuclear con una imagen de núcleos segmentados, para poder observar y comparar que tan óptima ha sido la segmetación.
* Que las regiones de interés hayan sido segmentadas correctamente.
Para utilzar esta herramienta los pasos son:
1. Ir a la carpeta __tables_graphics__
2. Colocar ahí las imágenes de interés.
3. Abrir el notebook _Notebook_Transparencia.ipynb_
4. Especificar las dos imágenes que se desean comparar y correr los bloques de código.
5. Para ajustar el nivel de transparencia basta con deslizar el btón que aparece al correr los bloques.


Autores:
* [Edgar Lara Arellano](https://github.com/Edgar-La)
* [David Ángel Alba Bonilla](https://github.com/DavidAlba2627)
* [Uriel Cárdenas Aguilar](https://github.com/Uriel148)
