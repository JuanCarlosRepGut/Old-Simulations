import numpy as np
import matplotlib.pyplot as plt
#Leemos los datos del archivo solapamiento.dat
nombre_archivo= 'solapamientoT.dat'
#Abrimos el archivo
with open(nombre_archivo, 'r') as archivo:
    datos=archivo.readlines()

#Arrays donde guardamos los datos
temp=[]
solapamiento=[]
for linea in datos:
    temp_i, solapamiento_i=linea.strip().split(',')
    temp.append(float(temp_i))
    solapamiento.append(float(solapamiento_i))

temp=np.array(temp)
solapamiento=np.array(solapamiento)
#Datos que se grafican y estilo
plt.plot(temp,solapamiento, marker='.', linestyle='-',color='r')
#Nombramos los ejes, el grafico y activamos la cuadricula
plt.xlabel('Temperatura')
plt.ylabel('Solapamiento')
plt.title('Solapamiento en funcion de la temperatura')
plt.grid(True)

#Mostramos la grafica
plt.show()