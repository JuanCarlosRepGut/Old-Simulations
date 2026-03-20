import numpy as np
import matplotlib.pyplot as plt


nombre_archivo='solapamientoD.dat'
#Los comentarios son similares a los anteriores (Apartado aleatorio)
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

    if pasos and solapamiento:
         plt.plot(pasos, solapamientos,  marker='.', linestyle='-')
         
    
    plt.xlabel('Pasos Monte-Carlo')
    plt.ylabel('Solapamiento')
    plt.title('Solapamiento frente a pasos Monte-Carlo. Patrón Deformado.')
    plt.grid(True)
    plt.xlim(0,10)

    etiquetas=["D=0.1", "D=0.2", "D=0.3", "D=0.4"]
    plt.legend(etiquetas)

    plt.show()