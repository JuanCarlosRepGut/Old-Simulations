import numpy as np
import matplotlib.pyplot as plt

#Nombre del archivo que vamos a abrir y leer los datos
nombre_archivo='solapamientovT.dat'

#Arrays donde guardamos los datos
pasos=[]
solapamientos=[]
#Abrimos el archivo
with open(nombre_archivo, 'r') as archivo:
    datos=archivo.readlines()

    #Se tiene en cuenta que el archivo tiene varios sets de datos diferentes que graficar en la misma figura
    for linea in datos:
            if linea.strip()== "":
                if pasos and solapamientos:
                    #Datos que se grafican y estilo
                    plt.plot(pasos, solapamientos, marker='.', linestyle='-')

                pasos=[]
                solapamientos=[]
            else:
                datos=linea.strip().split(",")
                paso=int(datos[0])
                solapamiento=float(datos[1])
                pasos.append(paso)
                solapamientos.append(solapamiento)
#Datos que se grafican y estilo
    if pasos and solapamientos:
         plt.plot(pasos, solapamientos, marker='.', linestyle='-')
    
    #Nombramos los ejes, el grafico y activamos la cuadricula
    plt.xlabel('Pasos Monte-Carlo')
    plt.ylabel('Solapamiento')
    plt.title('Solapamiento frente a pasos Monte-Carlos a diferentes temperaturas')
    plt.grid(True)

    #Escribimos la leyenda y la ubicamos en el grafico
    etiquetas=["T=0.0001", "T=0.01", "T=0.02", "T=0.03", "T=0.04", "T=0.05", "T=0.06", "T=0.07"]
    plt.legend(etiquetas,loc='lower left')
    
    
    #Mostramos la grafica
    plt.show()


