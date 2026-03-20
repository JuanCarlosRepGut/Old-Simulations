import numpy as np
import matplotlib.pyplot as plt

#Mismos comentarios que en programas similares-
#Leemos los datos del archivo solapamientoT.dat

nombre_archivo= 'solapamientoT.dat'
temperaturas=[]
solapamientos=[]
with open(nombre_archivo, 'r') as archivo:
    datos=archivo.readlines()

    for linea in datos:
            if linea.strip()== "":
                if temperaturas and solapamientos:
                    plt.plot(temperaturas, solapamientos, marker='.', linestyle='-')

                temperaturas=[]
                solapamientos=[]
            else:
                datos=linea.strip().split(",")
                temperatura=float(datos[0])
                solapamiento=float(datos[1])
                temperaturas.append(temperatura)
                solapamientos.append(solapamiento)

    if temperaturas and solapamientos:
         plt.plot(temperaturas, solapamientos, marker='.', linestyle='-')


plt.xlabel('Temperatura')
plt.ylabel('Solapamiento')
plt.title('Solapamiento en funcion de la Temperatura')
plt.grid(True)

etiquetas=["Patrón 1","Patrón 2","Patrón 3","Patrón 4"]
plt.legend(etiquetas,loc='best')

plt.show()