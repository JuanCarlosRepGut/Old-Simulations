import numpy as np
import matplotlib.pyplot as plt

#Igual que en programas de graficación anteriores
nombre_archivo='solapamientoD2.dat'

pasos=[]
solapamientos=[]

with open(nombre_archivo, 'r') as archivo:
    datos=archivo.readlines()

    for linea in datos:
            if linea.strip()== "":
                if pasos and solapamientos:
                    plt.plot(pasos, solapamientos, marker='.',linestyle='-')

                pasos=[]
                solapamientos=[]
            else:
                datos=linea.strip().split(",")
                paso=float(datos[0])
                Solapamiento=float(datos[1])
                pasos.append(paso)
                solapamientos.append(Solapamiento)

    if pasos and solapamientos:
         plt.plot(pasos, solapamientos, marker='.',linestyle='-')
    
    plt.xlabel('Temperatura')
    plt.ylabel('Solapamiento')
    plt.title('Solapamiento a diferentes temperaturas. Patron Deformado.')
    plt.grid(True)
   

    etiquetas=["T=0.015; D=0.1", "T=0.03; D=0.1", "T=0.045; D=0.1", "T=0.06; D=0.1", "T=0.015; D=0.2", "T=0.03; D=0.2", "T=0.045; D=0.2", "T=0.06; D=0.2", "T=0.015; D=0.3", "T=0.03; D=0.3", "T=0.045; D=0.3", "T=0.06; D=0.3"]
    #La leyenda ocupa mucho espacio, con esta linea la reubicamos para que no estorbe
    plt.legend(etiquetas,loc='upper right', bbox_to_anchor=(1.4,1))
    
    

    plt.show()