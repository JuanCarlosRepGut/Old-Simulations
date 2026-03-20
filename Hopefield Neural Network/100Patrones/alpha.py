import numpy as np
import matplotlib.pyplot as plt

#nombre del archivo
nombre_archivo= 'alpha.dat'
#abrimos el archivo con los datos
with open(nombre_archivo, 'r') as archivo:
    datos=archivo.readlines()

#Definimos Arrays donde guardaremos los datos del fichero
PA=[] #Patron Almacenado
alpha=[] #Fracción de patrones recordados alpha
for linea in datos:
    PA1, PC1=linea.strip().split(',') #Los datos vienen separados por comas
    PA.append(float(PA1)) #Añadimos los datos al array.
    alpha.append(float(PC1))

PA=np.array(PA) #Los arrays con todos los datos
alpha=np.array(alpha)

#Mostramos en la gráfica los datos.
plt.plot(PA,alpha, marker='.', linestyle='',color='r')
plt.xlabel('Patrones almacenados')
plt.ylabel('Patrones recordados')
plt.title('Patrones recordados frente patrones almacenados.')
plt.grid(True)

plt.show()