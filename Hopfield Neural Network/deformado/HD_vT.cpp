#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"

gsl_rng *tau;

#define MAX 30
#define MAX_ITER 10000

void randspininicial(int s[][MAX],int N, double deformacion);
double funcion_solapamiento(int s[][MAX], int N, int xi[][MAX], double a);
//Comentarios similares a Hopfield_variasT.cpp /Un_Patron/Aleatorio.
int main(void)
{
    int i, j, k, l, n, m, N, iter, paso,t, t1,pasos;
    char c;
    double a, T, y, p, difH, x, suma, solapamiento,solapamientomedio, deformacion;
    int s[MAX][MAX], xi[MAX][MAX];
    double w[MAX][MAX][MAX][MAX], theta[MAX][MAX];
    FILE *fspin, *fpatron, *fsolapamiento, *fsolapamientoT;
    extern gsl_rng *tau;
    int semilla=1395734759;
    
    tau=gsl_rng_alloc(gsl_rng_taus); 
    gsl_rng_set(tau,semilla);

    N=30;
    pasos=50;
    

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
    
fclose(fpatron);  

//Calculamos a

a=0.0;
for(i=0;i<N;i++)
{
    for(j=0;j<N;j++)
        a+=xi[i][j];
}
a/=(N*N);

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
                else w[i][j][k][l]=(xi[i][j]-a)*(xi[k][l]-a)/(N*N);
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


fspin=fopen("ising_data.dat","w");
fsolapamiento=fopen("solapamientoD2.dat","w");
//Estudiamos solo 3 deformaciones para distinguirlas debidamente en la grafica de datos
for(t1=1;t1<=3;t1++)
    {

for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            s[i][j]=xi[i][j];
    }
    
    deformacion=0.1*t1;
    randspininicial(s,N,deformacion);
    //Tomamos 4 temperaturas diferentes. Esta vez obviamos T=0.0001
    for(t=1;t<5;t++)
    {
        T=0.015*t;

for(iter=1;iter<=pasos*N*N;iter++)
{

    paso=iter/(N*N);
    n=gsl_rng_uniform_int(tau,N);
    m=gsl_rng_uniform_int(tau,N);

    suma=0.0;
    for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                suma+=(w[i][j][n][m]+w[n][m][i][j])*(1-2*s[n][m])*s[i][j];
            }
        }
    difH=-suma/2 + (theta[n][m]*(1-2*s[n][m]));

    //Calculamos Solapamiento

    p=exp(-difH/T);
    if(p>1.0) p=1.0;

    x=gsl_rng_uniform(tau);

    if(x<p)
        s[n][m]=1-s[n][m];


    if((iter%(N*N)==0))
    {
        solapamiento=funcion_solapamiento(s,N,xi,a);
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
    
    fprintf(fsolapamiento, " %i,%f\n", paso, solapamiento);
    
    }    
}
    fprintf(fspin, "\n");
    fprintf(fsolapamiento, "\n");
    }
    
    }
fclose(fspin);
fclose(fsolapamiento);

    return 0;
}

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
    {   //Genereamos un número real aleatorio, si este es menor que la deformación
        //Cambiamos el estado de la neurona i,j.
        x=gsl_rng_uniform(tau);
        if(x<deformacion)
        {
            s[i][j]=1-s[i][j];
        }
    }

    return;
}

double funcion_solapamiento(int s[][MAX], int N, int xi[][MAX], double a)
{
    int i, j;
    double m;
    m=0.0;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            m+=(xi[i][j]-a)*(s[i][j]-a);
        }
    }

    return m/(N*N*a*(1.0-a));
}