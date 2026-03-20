import matplotlib.pyplot as plt
from scipy import optimize

nombre_archivo='solapamientovT.dat'

datos_enteros=[]
datos_reales=[]

with open(nombre_archivo, 'r') as archivo:
    lineas=archivo.readlines()

    for linea in lineas:
            if linea.strip()== "":
                if datos_enteros and datos_reales:
                    plt.plot(datos_enteros, datos_reales, marker='.', linestyle='-')

                datos_enteros=[]
                datos_reales=[]
            else:
                datos=linea.strip().split(",")
                entero=int(datos[0])
                real=float(datos[1])
                datos_enteros.append(entero)
                datos_reales.append(real)

    if datos_enteros and datos_reales:
         plt.plot(datos_enteros, datos_reales, marker='.', linestyle='-')
    
    plt.xlabel('Pasos Monte-Carlo')
    plt.ylabel('Solapamiento')
    plt.title('Solapamiento frente a pasos Monte-Carlos a diferentes temperaturas')
    plt.grid(True)

    etiquetas=["T=0.0001", "T=0.01", "T=0.02", "T=0.03", "T=0.04", "T=0.05", "T=0.06", "T=0.07"]
    plt.legend(etiquetas,loc='lower left')
    
    

    plt.show()


