#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"


gsl_rng *tau;

#define MAX 30
#define MAXP 10
#define MAX_ITER 100000

void randspininicial(int s[][MAX],int N, double deformacion);
void funcion_solapamiento(int s[][MAX], int N, int xi[][MAX][MAXP], double a[MAXP],double solapamiento[MAXP],int P);
//Comentarios similares a los codigos de Un Patron/ Deformado.
int main(void)
{
    int i, j, k, l, n, m, N, iter, paso, pasos, mu, P,t; //P=dimension array de patrones (Numero de patrones que hay en la red) mu=contador de nºpatrones
    char c;
    double a[MAXP], p, difH, x, suma, solapamiento[MAXP], sumamu, deformacion; //sumamu variable auxiliar para el calculo de w
    int s[MAX][MAX], xi[MAX][MAX][MAXP];
    double w[MAX][MAX][MAX][MAX], theta[MAX][MAX], T;
    FILE *fspin, *fpatron, *fpatron1 , *fpatron2 , *fpatron3, *fsolapamiento;
    extern gsl_rng *tau;
    int semilla=1395734759;
    
    tau=gsl_rng_alloc(gsl_rng_taus); 
    gsl_rng_set(tau,semilla);

    N=30;
    P=4;
    T=0.0001;
    pasos=50;


//Abrimos los archivos que guardan los patrones.

fpatron=fopen("patron.txt", "r");
fpatron1=fopen("patron1.txt", "r");
fpatron2=fopen("patron2.txt", "r");
fpatron3=fopen("patron3.txt", "r");

    //Guardamos en el array xi los patrones.xi[i][j][mu] donde i representa la fila, j la columna, y mu el número del patron que se guarda.
   for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            fscanf(fpatron, "%c", &c);
            while(c== ' ' || c== '\n')
                fscanf(fpatron, "%c", &c);
            xi[i][j][0]=c-'0';
        }
    }

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            fscanf(fpatron1, "%c", &c);
            while(c== ' ' || c== '\n')
                fscanf(fpatron1, "%c", &c);
            xi[i][j][1]=c-'0';
        }
    }

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            fscanf(fpatron2, "%c", &c);
            while(c== ' ' || c== '\n')
                fscanf(fpatron2, "%c", &c);
            xi[i][j][2]=c-'0';
        }
    }

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            fscanf(fpatron3, "%c", &c);
            while(c== ' ' || c== '\n')
                fscanf(fpatron3, "%c", &c);
            xi[i][j][3]=c-'0';
        }
    }
    
//Cerramos los ficheros
fclose(fpatron);
fclose(fpatron1);
fclose(fpatron2);
fclose(fpatron3);  
      
      //Abrimos los ficheros donde se guarda el estado de la red y el solapamiento
fspin=fopen("ising_data.dat","w");
fsolapamiento=fopen("solapamientoD.dat","w");
 //Calculamos a


for(mu=0;mu<P;mu++)
{
    a[mu]=0.0;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            a[mu]+=xi[i][j][mu];
    }
    a[mu]=a[mu]/(N*N);
    //printf("%f",a[mu]);
    //printf("%c", '\t');
}

//Calculamos matriz w
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
                    for(mu=0;mu<P;mu++)
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
        //printf("%f",theta[i][j]);
    }
    //printf("%c",'\n');
}

//Bucle que recorre los distintos parametros de deformacion, esta vez se han tomado solo dos valores distintos en vez de 4 para ahorra tiempo de compilación
    for(t=1;t<=2;t++)
    {
    deformacion=0.2*t;

for(mu=0;mu<P;mu++)
{
    //Iniciamos la red a una configuracion aleatoria:
    //Primero ponemos en 1 todas las neuronas
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            s[i][j]=xi[i][j][0]; //Cambiando el valor de la tercera componente por 0, 1, 2 o 3 da como resultado que la configuración inicial sea cualquiera de esos patrones deformados
    }
    //Con esta función se cambian aleatoriamente los estados de 1 a 0.
    randspininicial(s,N,deformacion);

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
        
     paso=iter/(N*N);
     funcion_solapamiento(s,N,xi,a,solapamiento,P);
    
     fprintf(fsolapamiento, " %i,%f\n", paso, solapamiento[mu]);
     if(paso==pasos) fprintf(fsolapamiento,"%c", '\n');

}
}
}
}

//Cerramos los ficheros
fclose(fspin);
fclose(fsolapamiento);
    return 0;
}

//Función que cambia aleatoriamente el estado de la red. Toda la red debe estar en estado activo.
void randspininicial(int s[][MAX],int N, double deformacion)
{
    int i, j;
    int n;
    double x;
    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=1395734759;
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);
    
    for(i=0;i<N;i++)
    for(j=0;j<N;j++)
    {
        x=gsl_rng_uniform(tau);
        if(x<deformacion)
        {
            s[i][j]=1-s[i][j];
        }
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
