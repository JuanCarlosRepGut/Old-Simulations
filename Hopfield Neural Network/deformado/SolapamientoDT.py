import numpy as np
import matplotlib.pyplot as plt
#Igual que en programas de graficación anteriores
nombre_archivo='solapamientoDvT.dat'

temperaturas=[]
solapamientos=[]

with open(nombre_archivo, 'r') as archivo:
    datos=archivo.readlines()

    for linea in datos:
            if linea.strip()== "":
                if temperaturas and solapamientos:
                    plt.plot(temperaturas, solapamientos, marker='.',linestyle='-')

                temperaturas=[]
                solapamientos=[]
            else:
                datos=linea.strip().split(",")
                Temperatura=float(datos[0])
                Solapamiento=float(datos[1])
                temperaturas.append(Temperatura)
                solapamientos.append(Solapamiento)

    if temperaturas and solapamientos:
         plt.plot(temperaturas, solapamientos, marker='.',linestyle='-')
    
    plt.xlabel('Temperatura')
    plt.ylabel('Solapamiento')
    plt.title('Solapamiento a diferentes temperaturas. Patron Deformado.')
    plt.grid(True)

    etiquetas=["D=0.10", "D=0.20", "D=0.30", "D=0.40"]
    plt.legend(etiquetas)

    plt.show()