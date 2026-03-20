#include <stdio.h>
#include <math.h>

        
#define G 6.67384E-11     // Constante de gravitación universal 
#define PI 3.14159265359    //Numero Pi

#define dST 1.496E11      // Distancia entre la Tierra y el Sol en m
#define dSM 2.279E11       // Distancia entre Marte y el Sol en m


#define MS 1.989E30         // Masa del Sol en kg
#define MT 5.97E24         // Masa de la Tierra en kg
#define MM 6.42E23         // Masa de Marte en kg


//Se obtiene a partir de w = 2pi/T 
#define wT 1.990595E-7       // Velocidad angular de la Tierra (rad/s)
#define wM 1.06224E-7        // Velocidad angular de Marte (rad/s)


#define desfase 0.79            // desfase inicial entre Marte y la Tierra. Valores utilizados (0.79 (Hohmann) ; 0.678 (orbita 2) ; 0.633 (orbita 3))
#define RT 6.37814E6          // Radio de la Tierra en m
#define m 10000                //Masa del cohete

//Funciones auxiliares para la resolución númerica del sistema de ecuaciones diferencias mediante Runge-Kutta orden 4.
double fr(double f[5]);     //Ecuación diferencial de r
double fphi(double f[5]);   //Ecuación diferencial de phi
double fpr(double f[5]);    //Ecuación diferencial de momento en dirección radial (pr)
double fpphi(double f[5]);  //Ecuación diferencial de momento en dirección tangencial (pphi)

int main(void)
{
    //El vector f guarda los valores de las ecuaciones de movimiento r, phi, pr y pphi y el tiempo en ese orden
    //El vector aux es un vector auxiliar que simplifica el calculo del algoritmo de Runge-Kutta
    //Los vectores x e y guardan las posiciones en cartesianas del Sol, Tierra, Marte y cohete respectivamente. (Animación trayectoria)
    //Los vectores xM e yM guardan  las posiciones en cartesianas de Marte y la nave desde el sistema de referencia centrado en Marte (Animación órbita)
    //La matriz k guarda los valores de cada sumando del método runge kutta para cada ecuación.
    double f[5], aux[5], x[4], y[4], xM[2], yM[2], k[4][4]; 
    //h indica el paso del método Runge-Kutta. Tiene unidades de segundos.
    double h;
    //Theta es el ángulo con el que sale disparado inicialmente el cohete, vr la velocidad radial, vtheta la velocidad tangencial y v el módulo de la velocidad de despuege.
    double theta, vr, vtheta, v;
    //Contadores: i cuenta las iteraciones, j recorre las posiciones de los vectores x,y,xM y yM, n es el número total de iteraciones y guardar sirve para guardar los datos en ficheros cada "guardar" iteraciones.
    int i, j, n, guardar;
    //Valores booleanos para detectar si el planeta se acerca a Marte y si orbita.
    bool orbita, entrada;
    //r_orbita es la distancia a partir de cual se activan los propulsores, disM la distancia de la nave a Marte, dis_entada es la distancia en la cual se considera que la nave esta cerca de Marte
    //vimpulso1 y vimpulso2 son el modulo de la velocidad de cada uno de los impulsos necesarios para orbitar Marte (Mirar sección transferencia de Hohmann de la memoria)
    //vx y vy son la velocidad de impulso esta vez expresadas en coordenadas cartesianas.
    double r_orbita, disM,dis_entrada, vimpulso1,vimpulso2,vx,vy;
    //T_0 y T_f son la energía cinética de la nave incial y final respectivamente. E y DifE lo mismo para la energía total. Sirven para el balance energético de la misión.
    double T_0, T_f, E, DifE;
    //Ficheros donde guardamos todos los datos. Posicion respecto al Sol y Marte
    FILE *fdatos, *fenergia, *forbitaM;
    //Abrimos los ficheros
    fdatos=fopen("planets_data.dat", "w");
    fenergia=fopen("Energia.dat", "w");
    forbitaM=fopen("orbita_marte.dat", "w");
    //Damos valores a los parametros iniciales:
    n=4000000;
    h=30;
    orbita=false;
    entrada=false;
    guardar=20000;

    //Vaalores iniciales de r, phi, pr y pphi

    f[0]=(dST+RT)/dST; //r
    f[1]=0.0; //phi
    
    vr=0.0; //velocidad radial inicial de cohete
    vtheta=32800.0; //velocidad tangencial inicial de cohete.Valores utilizados (32800 (Hohmann) ; 35000 (orbita 2) ; 37500 (orbita 3))
    theta=asin((vtheta)/(sqrt(vr*vr+vtheta*vtheta))); //Ángulo inicial de salida del cohete
    v=sqrt(vr*vr+vtheta*vtheta)/dST; //Módulo.
    f[2]=v*cos(theta-f[1]); //pr
    f[3]=v*f[0]*sin(theta-f[1]); //pphi
    f[4]=0.0; //t

    r_orbita=4.5E8/dST; //Valores utilizados (4.5E8/dST (Hohmann) ; 5.0E8/dST (orbita 2) ; 5.5E8/dST (orbita 3))
    
    DifE=0.0;
    T_0=m*(dST*wT*dST*wT)/2;
    T_f=m*(f[2]*f[2]+(f[3]*f[3])/(f[0]*f[0]))*dST*dST/2;
    E=abs(T_0-T_f);
    DifE+=E;
    printf("Energia despegue:" "%f\n", E);

    for(i=0;i<n;i++)
    {
        //Calcualmos k1
        k[0][0]=h*fr(f);
        k[0][1]=h*fphi(f);
        k[0][2]=h*fpr(f);
        k[0][3]=h*fpphi(f);

        //Calculamos k2
        for (j=0;j<4;j++)
        {
            aux[j]=f[j]+k[0][j]/2.0;
        }
        aux[4]=f[4]+h/2;

        k[1][0]=h*fr(f);
        k[1][1]=h*fphi(f);
        k[1][2]=h*fpr(f);
        k[1][3]=h*fpphi(f);

        //Calculamos k3
        for (j=0;j<4;j++)
        {
            aux[j]=f[j]+k[1][j]/2.0;
        }
        aux[4]=f[4]+h/2;;

        k[2][0]=h*fr(f);
        k[2][1]=h*fphi(f);
        k[2][2]=h*fpr(f);
        k[2][3]=h*fpphi(f);

        //Calculamos k4
        for (j=0;j<4;j++)
        {
            aux[j]=f[j]+k[2][j];
        }
       aux[4]=f[4]+h;

        k[3][0]=h*fr(f);
        k[3][1]=h*fphi(f);
        k[3][2]=h*fpr(f);
        k[3][3]=h*fpphi(f);

        //Avanzamos al siguiente paso temporal

        for(j=0;j<4;j++)
        {
            f[j]+=(k[0][j]+2.0*k[1][j]+2.0*k[2][j]+k[3][j])/6.0;
        }
        f[4]+=h;

            //Guardamos en en los vectores las posiciones de la siguiente iteración.
            x[0]=y[0]=0.0;
            x[1]=cos(wT*f[4]); y[1]= sin(wT*f[4]);
            x[2]=(dSM/dST)*cos(wM*f[4]+desfase); y[2]=(dSM/dST)*sin(wM*f[4]+desfase);
            x[3]=f[0]*cos(f[1]); y[3]=f[0]*sin(f[1]);

            xM[0]=yM[0]=0.0;
            xM[1]=f[0]*cos(f[1])-(1.52)*cos(wM*f[4]+desfase);
            yM[1]=f[0]*sin(f[1])-(1.52)*sin(wM*f[4]+desfase);

        //Guardamos los datos tras "guardar" iteraciones.
        if(i%guardar==0)
        {
            for(j=0;j<4;j++)
            {
                fprintf(fdatos, "%f,%f\n", x[j], y[j]);
                if(j!=0&&(j%3==0))
                fprintf(fdatos, "\n");
            }

            for(j=0;j<2;j++)
            {

                fprintf(forbitaM, "%f,%f\n", xM[j], yM[j]);
                if(j==1)
                fprintf(forbitaM, "\n");
            }

            fprintf(fenergia,"%f\t,%f\n", f[4], DifE);
        }
        //Cuando el cohete se acerca a Marte cambiamos el paso temporal.

        disM=sqrt((1.52*1.52)+f[0]*f[0]-2*1.52*f[0]*cos(f[1]-wM*f[4]-desfase));
        dis_entrada=1.0E9/dST;
        if((entrada==false)&&(disM<dis_entrada)) //En el primer momento que el cohete se encuentra cercana a Marte se cambia el paso temporal y se imprime el tiempo que se ha tardado en llegar en días.
        {
            h=10;
            guardar=15000;
            entrada=true;
            printf("Duracion del viaje en dias:\t"); printf("%f\n", (f[4]/86400));
        }

        if((orbita==false)&&(disM<r_orbita)) //En el primer momento que el cohete esta a una distancia valida para orbitar se da el impulso y se hace el balance energético.
        {
            //Energia cinetica inical antes del impulso
            T_0=m*(f[2]*f[2]+(f[3]*f[3])/(f[0]*f[0]))*dST*dST/2;

            //Modulo de velocidades del impulso
            vimpulso1=sqrt(G*MM/(dST*disM))/dST; 
            vimpulso2=dSM/dST*wM;

            //velocidad de impulso en cartesianas:

            vx=vimpulso1*(1.52*sin(wM*f[4]+desfase)-f[0]*sin(f[1]))/disM-vimpulso2*sin(wM*f[4]+desfase);
            vy=-vimpulso1*(1.52*cos(wM*f[4]+desfase)-f[0]*cos(f[1]))/disM+vimpulso2*cos(wM*f[4]+desfase);
        
            f[2]=vx*cos(f[1])+vy*sin(f[1]);
            f[3]=f[0]*(-vx*sin(f[1])+vy*cos(f[1]));

            T_f=m*(f[2]*f[2]+(f[3]*f[3])/(f[0]*f[0]))*dST*dST/2;
            E=abs(T_f-T_0);
            printf("Energia pulso:" "%f\n", E);
            DifE+=E;

            orbita=true;
        }
    }
    //Cerramos los ficheros.
    fclose(fdatos);
    fclose(fenergia);
    fclose(forbitaM);

    return 0;
}

double fr(double f[5])
{
    return f[2];
}

double fphi(double f[5])
{
    return (f[3]/(f[0]*f[0]));
}

double fpr(double f[5])
{
    double disT,disM,disT3,disM3,muM,muT,Delta; //Valores auxiliares para hacer el calculo más limpio.
    double term1,term2,term3,term4,resultado;

    disT=sqrt((1)*(1)+f[0]*f[0]-2*f[0]*(1)*cos(f[1]-wT*f[4])); //Distancia cohete-Tierra
    disM=sqrt((1.52)*(1.52)+f[0]*f[0]-2*(1.52)*f[0]*cos(f[1]-wM*f[4]-desfase)); //Distancai cohete-Marte
    disT3=pow(disT,3);
    disM3=pow(disM,3);
    muM=MM/(MS);
    muT=MT/(MS);
    Delta=(G*MS)/(dST*dST*dST);
    muM=3.21267E-7;
    muT=3.002514E-9;

    term1=(f[3]*f[3])/(f[0]*f[0]*f[0]);
    term2=1.0/(f[0]*f[0]);
    term3=((muT)*(f[0]-(1)*cos(f[1]-wT*f[4])))/(disT3);
    term4=((muM)*(f[0]-(1.52)*cos(f[1]-wM*f[4]-desfase)))/(disM3);
    
    resultado=term1-Delta*(term2+term3+term4);

    return resultado;

}

double fpphi(double f[5]){
    double disT,disM,disT3,disM3,muM,muT,Delta; //Valores auxiliares para hacer el calculo más limpio.
    double term1,term2,resultado;

    disT=sqrt((1)*(1)+f[0]*f[0]-2*f[0]*(1)*cos(f[1]-wT*f[4])); //Distancia cohete-Tierra
    disM=sqrt((1.52)*(1.52)+f[0]*f[0]-2*(1.52)*f[0]*cos(f[1]-wM*f[4]-desfase)); //Distancai cohete-Marte
    disT3=pow(disT,3);
    disM3=pow(disM,3);

    Delta=(G*MS)/(dST*dST*dST);
    muM=3.21267E-7;
    muT=3.002514E-9;

    term1=(muT*(1)*f[0]*sin(f[1]-wT*f[4]))/(disT3);
    term2=(muM*(1.52)*f[0]*sin(f[1]-wM*f[4]-desfase))/(disM3);
    resultado=-Delta*(term1+term2);

    return resultado;

}