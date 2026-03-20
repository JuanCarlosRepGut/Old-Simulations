#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"


gsl_rng *tau;

#define MAX 30
#define MAX_ITER 10000

//Funciones auxiliares:
void randspininicial(int s[][MAX],int N); 
double funcion_solapamiento(int s[][MAX], int N, int xi[][MAX], double a);

int main(void)
{
    //Declaración de variables.
    int i, j, k, l, n, m, iter, paso; //Contadores que recorren las matrices s, w, theta y el bucle del algoritmo.
    int N,  pasos; //Dimension de la red y numero de pasos Monte Carlo.
    char c; // Esta variable caracter sirve para que el programa reconozca el archivo de entrada de los patrones.
    double a, p, difH, x, suma, solapamiento; //Variables reales para hacer calculos. a=fracción de neuronas en estado activo o inactivo; p=probabiliddad de Boltzman;
                                              //difH=diferencia del hamiltoniano entre iteraciones; x=numero aleatorio para comprobar si se hace el cambio;
                                              //suma=variable auxiliar para el calculo e difH;solapamiento=solapamiento e la red, cuanto se parece al patrón almacenado.
    int s[MAX][MAX], xi[MAX][MAX]; //Matrices de enteros. s representa el estado de la red y xi el patrón almacenado
    double w[MAX][MAX][MAX][MAX], theta[MAX][MAX], T; //Arrays de reales. w representa los pesos sinápticos y theta la activación de las neuronas. T es la temperatura del sistema.
    FILE *fspin, *fpatron, *fsolapamiento; //Punteros a los ficheros
    extern gsl_rng *tau; //Puntero para generar los números aleatorios
    int semilla=1395734759; //Se elige una semilla para los numeros aleatorios
    
    tau=gsl_rng_alloc(gsl_rng_taus); 
    gsl_rng_set(tau,semilla); //Se establece dicha semilla

    N=30; //Dimension de red 30*30
    T=0.0001; //Temperatura=0.0001
    pasos=50; //pasos Monte Carlo del algoritmo

    //Llamamos a la función que establece una configuración inicial.
    randspininicial(s,N);

//Abrimos el fichero donde se guarda el patrón, leemos los datos y lo guardamos en un array. Debemos pasar de caracteres a numeros enteros.
fpatron=fopen("patron.txt", "r");
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            fscanf(fpatron, "%c", &c);
            while(c== ' ' || c== '\n')
                fscanf(fpatron, "%c", &c);
            xi[i][j]=c-'0';
        }
    }
//Cerramos el fichero para no tener problemas de memoria.    
fclose(fpatron);  

 //Calculamos a según las expresiones teóricas.

a=0.0;
for(i=0;i<N;i++)
{
    for(j=0;j<N;j++)
        a+=xi[i][j];
}
a/=(N*N);

//Calculamos matriz w según las expresiones teóricas.

for(i=0;i<N;i++)
{
    for(j=0;j<N;j++)
    {
        for(k=0;k<N;k++)
        {
            for(l=0;l<N;l++)
            {
                if((i==k)&&(j==l)) w[i][j][k][l]=0.0;
                else w[i][j][k][l]=(xi[i][j]-a)*(xi[k][l]-a)/(N*N);
            }
        }
    }
}

//Calculamos matriz theta según las expresiones teóricas.
//Iniciamos a 0 porque theta se obtiene a partir de una suma.
for(i=0;i<N;i++)
{
    for(j=0;j<N;j++)
    {
        theta[i][j]=0.0;
    }
}

for(i=0;i<N;i++)
{
    for(j=0;j<N;j++)
    {
        for(k=0;k<N;k++)
        {
            for(l=0;l<N;l++)
                theta[i][j]+=w[i][j][k][l]/2;
        }
    }
}

//Abrimos los ficheros donde guardamos los datos del estado de la red y el solapamiento.
fspin=fopen("ising_data.dat","w");
fsolapamiento=fopen("solapamiento.dat","w");

//Comenzamos el bucle iterativo desde la primera iteración hasta el número de pasos Monte Carlo que se requiera.
for(iter=0;iter<=pasos*N*N;iter++)
{
    n=gsl_rng_uniform_int(tau,N); //Elegimos un punto aleatorio de la red.
    m=gsl_rng_uniform_int(tau,N);

    //Calculamos el valor de difH en ese punto.
    suma=0.0;
    for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                suma+=(w[i][j][n][m]+w[n][m][i][j])*(1-2*s[n][m])*s[i][j];
            }
        }
    difH=-suma/2 + (theta[n][m]*(1-2*s[n][m]));

    //Con difH podemos obtener la probabilidad de boltzmann.

    p=exp(-difH/T);
    if(p>1.0) p=1.0; //La probabilidad es como máximo 1.

    //Generamos un numero aleatorio para comparar.
    x=gsl_rng_uniform(tau);

    //Si el numero aleatorio x es menor que la probabilidad p cambiamos el estado de la neurona (n,m). 
    if(x<p)
        s[n][m]=1-s[n][m];

    //Cuando pasa un numero de iteraciones equivalente a un paso Monte Carlo, guardamos los datos del estado de la red, calculamos el solapamiento y guardamos los datos
    //del solapamientoen función del paso. Cada uno en su respectivo fichero.
    if((iter%(N*N)==0))
    {
    for(i=0;i<N;i++)
    {
    for(j=0;j<N;j++)
        {
            fprintf(fspin, " %i", s[i][j]);
            if(j<N-1) fprintf(fspin, ",");
    }
    fprintf(fspin, "\n");
    if((i!=0)&&(i%(N-1)==0))
            fprintf(fspin, "\n");
    }
    solapamiento=funcion_solapamiento(s,N,xi,a);
     paso=iter/(N*N);
    fprintf(fsolapamiento, " %i,%f\n", paso, solapamiento);
    
    }    
}
//Cerramos ambos ficheros
fclose(fspin);
fclose(fsolapamiento);
    return 0;
}

//Función que genera una configuración inicial aleatoria.
void randspininicial(int s[][MAX],int N)

{
    int i, j;
    int n;
    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=1395734759;
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);
    
    //Para cada punto de la red se genera un numero aleatorio entre 0 y 1.
    for(i=0;i<N;i++)
    for(j=0;j<N;j++)
    {
        n=gsl_rng_uniform_int(tau,2);
        s[i][j]=n;
    }

    return;
}

//Función que calcula el solapamiento a partir de ña expresión teórica.
double funcion_solapamiento(int s[][MAX], int N, int xi[][MAX], double a)
{
    int i, j;
    double m;
    //Iniciamos el solapamiento a 0.
    m=0.0;
    //Calculamos la suma de los terminos.
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            m+=(xi[i][j]-a)*(s[i][j]-a);
        }
    }
    //Finalmente devolvemos el valor debidamente normalizado.
    return m/(N*N*a*(1.0-a));
}


