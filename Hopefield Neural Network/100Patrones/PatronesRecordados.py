import numpy as np
import matplotlib.pyplot as plt

#nombre del archivo
nombre_archivo= 'solapamiento.dat'
#abrimos el archivo con los datos
with open(nombre_archivo, 'r') as archivo:
    datos=archivo.readlines()

#Definimos Arrays donde guardaremos los datos del fichero
PA=[] #Patrones almacenados
PC=[] #Patrones recordados
for linea in datos:
    PA1, PC1=linea.strip().split(',') #Los datos se separan mediante comas
    PA.append(int(PA1)) #Añadimos los datos del fichero en los arrays
    PC.append(int(PC1))

PA=np.array(PA) #Obtenemos los arrays con los datos del fichero
PC=np.array(PC)

#Los mostramos en la grafica y definimos el estilo.
plt.plot(PA,PC, marker='.', linestyle='',color='b')
plt.xlabel('Patrones almacenados')
plt.ylabel('Patrones recordados')
plt.title('Patrones recordados frente patrones almacenados.')
plt.grid(True)

plt.show()