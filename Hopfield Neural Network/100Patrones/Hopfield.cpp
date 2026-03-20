#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"


gsl_rng *tau;

#define MAX 30
#define MAXP 100

//Ahora hay una nueva funcion que genera patrones aleatorios para almacenarnos en memoria.
void randspininicial(int s[][MAX],int N, int P);
void funcion_solapamiento(int s[][MAX], int N, int xi[][MAX][MAXP], double a[MAXP],double solapamiento[MAXP],int P);
void Patronaleatorio(int xi[][MAX][MAXP], int N, int P);

int main(void)
{
    int i, j, k, l, n, m, N,PA, iter, paso ,pasos, mu, P,u,patronrecordado; //PA es un contador que significa patrón almacenado. patronrecordado es un patron que cuantifica los patrones recordados
                                                                            // Y u es otro contador auxiliar que recorre el codigo para 1 patron almacenado, 2 patrones almacenados,... hasta P patrones almacenados.
    char c;
    double a[MAXP], y, p, difH, x, suma, solapamiento[MAXP], sumamu, alpha; //alpha es la fraccion de patrones recordados definida en la memoria.
    int s[MAX][MAX], xi[MAX][MAX][MAXP];
    double w[MAX][MAX][MAX][MAX], theta[MAX][MAX], T;
    FILE *fspin, *fsolapamiento, *falpha;
    extern gsl_rng *tau;
    int semilla=1395734759;
    
    tau=gsl_rng_alloc(gsl_rng_taus); 
    gsl_rng_set(tau,semilla);

    N=20; //Ahora la red se toma de 20*20 (400 neuronas).
    P=60; //Tomaremos un total de 60 patrones.
    T=0.0001; //Temperatura baja para que la red funcione correctamente
    pasos=50; //50 pasos para que la red tenga tiempo de estabilizarse.
//Generamos los 60 patrones aleatorios.

Patronaleatorio(xi,N,P);

randspininicial(s,N,P);

fspin=fopen("ising_data.dat","w");
fsolapamiento=fopen("solapamiento.dat","w");
falpha=fopen("alpha.dat","w");

//Bucle que recorre cada numero de patrones almacenados
for(PA=1;PA<=P;PA++)
{
    alpha=0.0; //Se inician los contadores.
    patronrecordado=0;
    randspininicial(s,N,P);
//Calculamos a para cada valor de patrones almacenados cada vez
for(mu=0;mu<PA;mu++) 
{
    a[mu]=0.0;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            a[mu]+=xi[i][j][mu];
    }
    a[mu]=a[mu]/(N*N);
}

//Calculamos matriz w para cada valor de patrones almacenados cada vez
for(i=0;i<N;i++)
{
    for(j=0;j<N;j++)
    {
        for(k=0;k<N;k++)
        {
            for(l=0;l<N;l++)
            {
                if((i==k)&&(j==l)) w[i][j][k][l]=0.0;
                else 
                {
                    sumamu=0.0;
                    for(mu=0;mu<PA;mu++)
                    {
                        sumamu+=((xi[i][j][mu]-a[mu])*(xi[k][l][mu]-a[mu]));
                    }
                    w[i][j][k][l]=sumamu/(N*N);
    
                }
            }
        }
    }
}

//Calculamos matriz theta

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

//Bucle de la red

for(iter=0;iter<=pasos*N*N;iter++)
{
    //Se elige un punto aleatorio de la red
    n=gsl_rng_uniform_int(tau,N);
    m=gsl_rng_uniform_int(tau,N);

    //Se calcula la diferencia de Hamiltoniano
    suma=0.0;
    for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                suma+=(w[i][j][n][m]+w[n][m][i][j])*(1-2*s[n][m])*s[i][j];
            }
        }
    difH=-suma/2 + (theta[n][m]*(1-2*s[n][m]));

    //Se calcula la probabilidad según la distribución de Boltzmann
    p=exp(-difH/T);
    if(p>1.0) p=1.0;

    //Numero aleatorio
    x=gsl_rng_uniform(tau);
    //Si el numero aleatorio es menor que la probabilidad se cambia el estado de la neurona n,m
    if(x<p)
        s[n][m]=1-s[n][m];

    //Cuando pasa un paso Monte Carlo guardamos en el fichero de spin el estado de la red
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

     if(iter==pasos*N*N)
     { 
     //Llamamos a la funcion solapamiento en funcion del numero de patrones almacenados.
     funcion_solapamiento(s,N,xi,a,solapamiento,PA);
     //Este bucle comprueba si almenos uno de los patrones se considera recordado y aumenta los contadores
     for(u=0;u<PA;u++)
     {
        if(abs(solapamiento[u])>=0.75)
            {
                patronrecordado++;
                alpha++;
                alpha/=(double)(N*N);
            }
     }
     
     }
        //Cuando se guardan los datos para un patrón se da un salto de linea.
        
}

} //Guardamos los datos en los ficheros.
fprintf(fsolapamiento, " %i,%i\n", PA, patronrecordado);
fprintf(falpha, "%i,%f\n", PA, alpha);
}


//Cerramos los ficheros
fclose(fspin);
fclose(fsolapamiento);
fclose(falpha);
    return 0;
}

//Función que cambia aleatoriamente el estado de la red. Toda la red debe estar en estado activo.
void randspininicial(int s[][MAX],int N, int P)

{
    int i, j, t,semillanueva;
    int n;
    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=1395734759;
    tau=gsl_rng_alloc(gsl_rng_taus);
    //Este bucle hace que la semilla que genera las configuraciones iniciales aleatorias cambie, de forma que en cada iteración la red inicial es diferente.
    for(t=1;t<=P;t++)
        semilla=int(semilla*1.01);
    gsl_rng_set(tau,semilla);
    
    for(i=0;i<N;i++)
    for(j=0;j<N;j++)
    {
        n=gsl_rng_uniform_int(tau,2);
        s[i][j]=n;
    }

    return;
}

//Esta función se encarga de calcular el solapamiento
void funcion_solapamiento(int s[][MAX], int N, int xi[][MAX][MAXP], double a[MAXP],double solapamiento[MAXP], int P)
{
    int i, j,k;
    for(k=0;k<P;k++)
    {
        solapamiento[k]=0.0;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            solapamiento[k]+=((xi[i][j][k]-a[k])*(s[i][j]-a[k]));
        }
    }
        solapamiento[k]=solapamiento[k]/(N*N*a[k]*(1-a[k]));
    }

    return;
}
//Función que genera un patrón aleatorio. Es muy similar a la funcion randspininicial. Pero sin la necesidad de cambiar la semilla.
void Patronaleatorio(int xi[][MAX][MAXP], int N, int P)
{
    int i, j,mu;
    int n;
    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=1395734759;
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);
    for(mu=0;mu<P;mu++)
    for(i=0;i<N;i++)
    for(j=0;j<N;j++)
    {
        n=gsl_rng_uniform_int(tau,2);
        xi[i][j][mu]=n;
    }

    return;
}