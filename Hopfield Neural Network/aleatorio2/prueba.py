import re
import matplotlib.pyplot as plt

# Nombre del archivo de datos
nombre_archivo = 'solapamientoT.dat'

# Leer los datos del archivo
with open(nombre_archivo, 'r') as archivo:
    lineas = archivo.readlines()

# Almacenar los datos en listas separadas
x = []
curva1 = []
curva2 = []
curva3 = []
curva4 = []

for linea in lineas:
    datos = re.split('\s*,\s*', linea.strip())
    x.append(float(datos[0]))
    curva1.append(float(datos[1]))
    curva2.append(float(datos[2]))
    curva3.append(float(datos[3]))
    curva4.append(float(datos[4]))

# Crear la gráfica
plt.plot(x, curva1, label='Curva 1')
plt.plot(x, curva2, label='Curva 2')
plt.plot(x, curva3, label='Curva 3')
plt.plot(x, curva4, label='Curva 4')

# Configurar la leyenda y los ejes
plt.legend()
plt.xlabel('Eje X')
plt.ylabel('Valores')

# Mostrar la gráfica
plt.show()