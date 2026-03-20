import matplotlib.pyplot as plt
import numpy as np

nombre_archivo='solapamiento.dat'

pasos=[]
solapamientos=[]

with open(nombre_archivo, 'r') as archivo:
    datos=archivo.readlines()

    for linea in datos:
            if linea.strip()== "":
                if pasos and solapamientos:
                    plt.plot(pasos, solapamientos, marker='.', linestyle='-')

                pasos=[]
                solapamientos=[]
            else:
                datos=linea.strip().split(",")
                paso=int(datos[0])
                solapamiento=float(datos[1])
                pasos.append(paso)
                solapamientos.append(solapamiento)

    if pasos and solapamientos:
         plt.plot(pasos, solapamientos, marker='.', linestyle='-')
    
    plt.xlabel('Pasos Monte-Carlo')
    plt.ylabel('Solapamiento')
    plt.title('Solapamiento frente a pasos Monte-Carlo para varios patrones. T=0.0001')
    plt.grid(True)
    plt.xlim(0,10)

    etiquetas=["Patrón 1","Patrón 2","Patrón 3","Patrón 4"]
    plt.legend(etiquetas,loc='best')
    
    

    plt.show()