import numpy as np
import matplotlib.pyplot as plt

#Nombre del archivo que vamos a abrir y leer los datos
nombre_archivo= 'solapamiento.dat'
#Abrimos el archivo
with open(nombre_archivo, 'r') as archivo:
    datos=archivo.readlines()

#Arrays donde guardamos los datos
paso=[]
solapamiento=[]
for linea in datos:
    paso, solapamiento=linea.strip().split(',')
    paso.append(int(paso))
    solapamiento.append(float(solapamiento))

paso=np.array(paso)
solapamiento=np.array(solapamiento)
#Datos que se grafican y estilo
plt.plot(paso,solapamiento, marker='s', linestyle='-',color='g')
#Nombramos los ejes, el grafico y activamos la cuadricula
plt.xlabel('Tiempo(paso Monte-Carlo)')
plt.ylabel('Solapamineto')
plt.title('Solapamiento en funcion de cada paso Monte-Carlo. T=0.0001')
plt.grid(True)
plt.xlim(0,10)

#Mostramos la grafica
plt.show()