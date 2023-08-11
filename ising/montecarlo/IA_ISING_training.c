#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


//#define termalizar
#define ising

const int L = 20;
const float beta_inicial = 0.0;   //muy por debajo de la beta crítica -> temperaturas altas (desorden)
const float paso_beta = 0.02;
const int Nbetas = 50; //Nbetas distintas
const int Nconfig = 3000; //configs distintas para cada beta
//const int Npasos = 105000;
const int Nter = 5000;

//Parisi
#define NormRANu (2.3283063671E-10F)
unsigned int irr[256];
unsigned char ind_ran = 0, ig1 = 0, ig2 = 0, ig3 = 0;
unsigned int ir1;

void ini_ran(int SEMILLA){
    int INI, FACTOR, SUM, i;
    srand(SEMILLA);
    INI = SEMILLA;
    FACTOR = 67397;
    SUM = 7364893;
    for (i=0; i<256; i++){
        INI = (INI*FACTOR+SUM);
        irr[i] = INI;
    }
    ind_ran = ig1 = ig2 = ig3 = 0;
}
double Random (){                    //genera numero Random
    double r;
    ig1 = ind_ran-24;
    ig2 = ind_ran-55;
    ig3 = ind_ran-61;
    irr[ind_ran] = irr[ig1]+irr[ig2];
    ir1 = (irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r = ir1*NormRANu;
    return r;
}
//fin Parisi

void inicializa_config(int *config){ // genera una configuración aleatoria(primera configuración)
    int i;
    double alea;

    for (i=0; i<(L*L); i++){
            alea = Random();

            if (alea<0.5)  config[i] = 1;
            else  config[i] = -1;
        }
}
//Definicion de contornos
    void Def_contornos(int L, int *xp, int *xn, int *yp, int *yn){
        int k = 0;
        //eje X positivo
            xp[L-1] = -(L-1);
            for (k=0; k<L-1; k++) xp[k] = 1;
        //eje X negativo
            xn[0] =(L-1);
            for (k=1; k<L; k++)  xn[k] = -1;
        //eje Y positivo
            yp[L-1] = -L*(L-1);
            for (k=0; k<L-1; k++) yp[k] = L;
        //eje Y negativo
            yn[0] = L*(L-1);
            for (k=1; k<L; k++)  yn[k] = -L;


    }
 int calcula_energia(int *s, int *xp, int *yp){
    int Energia=0;
    int j, i;
    for (j=0; j<L; j++){
        for (i=0; i<L; i++)
                Energia += (s[i+(j*L)]*(s[i+(j*L)+xp[i]] + s[i+(j*L)+yp[j]]));
    }
    Energia = -Energia;
    return Energia;
 }

 int calcula_magnetizacion(int *conifg){
     int i, magnet=0;
     for (i=0; i<(L*L); i++)
         magnet += conifg[i];
     return magnet;
 }

  void metropolis (int V, int *s, int *xp, int *yp, int *xn, int *yn, double *prob, int *energia, int *magnetizacion){
    int n=0, y, x, indice=0, DE=0;

    for (y=0; y<L; y++){ //Aplicamos metropolis a los puntos en orden normal (no aleatorio)
        for (x=0; x<L; x++){
            DE = s[n]*(s[n+xp[x]] + s[n+yp[y]] + s[n+xn[x]] + s[n+yn[y]])*2;
            indice = DE/4 + 2;

            if (Random()<prob[indice]){
                s[n] *= -1;
                *energia += DE;
                *magnetizacion += s[n] * 2;
            }
            n++;
        }
    }
}

double funcion (double Energia, double beta){
    return exp(-beta*Energia);
}


int main(){

    int config[(L*L)], xp[L], yp[L], xn[L], yn[L], energia=0, magnetizacion=0;
    for (int index=0; index<L;index++){
        xp[index]=0;
        yp[index]=0;
        xn[index]=0;
        yn[index]=0;
    }
    for (int index=0; index<L*L; index++) config[index]=0;
    double prob[5], beta = beta_inicial, t_0, t_f;
    int contador=0;

//inicializamos el generador Random
    for(int i=0; i<256; i++)  irr[i] = (rand()<<16)+rand();
    ini_ran(time(NULL));
    t_0 = time(NULL);
//


    FILE *conf;
    FILE *labels;
    FILE *temp;
    FILE *termal;
    FILE *energy;
    FILE *magnet;
    conf=fopen("Config_memoria.txt", "w");
    labels=fopen("Labels_memoria.txt", "w");
    temp=fopen("Temperature_memoria.txt", "w");
    energy=fopen("Energia_memoria.txt", "w");
    magnet = fopen("Magnet_memoria.txt", "w");
    Def_contornos(L, xp, xn, yp, yn);


    for (int j=0; j<Nbetas; j++){
        printf("Beta = %lf\n", beta);
        t_f=time(NULL);
        printf("Tiempo transcurrido desde inicio: %d segundos\n", (int)(t_f-t_0));

        for (int i=0; i<5; i++){//Calculamos los posibles valores de nuestra función, así evitamos muchos cálculos en el algoritmo de metrópolis
            prob[i] = funcion((-8+i*4), beta);
            //printf("Prob %d = %lf\n", i, prob[i]);
        }
        for (int k=0; k<Nconfig; k++){
            inicializa_config(config);
            energia = calcula_energia(config, xp, yp);
            magnetizacion = calcula_magnetizacion(config);
            //printf("Energia inicial: %lf   magnetizacion inicial: %lf\n", energia/(double)(2.0*L*L),magnetizacion/((double) (L*L)));
            int aleatorio = (int)(10000*(rand()/(double)(RAND_MAX+1)));
            //printf("Aleatorio: %d\n", aleatorio);

            for (int i=0; i<(Nter+aleatorio); i++) {
                metropolis((L*L), config, xp, yp, xn, yn, prob, &energia, &magnetizacion);
            }
            //printf("Acabado metropolis\n");
            for (int index=0; index<(L*L); index++){
                if (index==(L*L)-1) {
                    fprintf(conf, "%d", config[index]);
                }
                else{
                    fprintf(conf, "%d\t ", config[index]);
                }
                //printf("%d\t ", config[index]);
            //printf("Index: %d\n", index);
            }

            fprintf(conf, "\n");
            //printf("Escrito en el fichero\n");
            int ord_dis;
            if (beta<0.4406) ord_dis=1; //temp alta = desordenado (dis)
            else if (beta>0.4406) ord_dis=0; //temp baja = ordenado (ord)
           // printf("Ord_dis: %d\n", ord_dis);
            fprintf(labels, "%d\n ", ord_dis);
            fprintf(temp, "%lf\n", beta);
            fprintf(energy, "%g\n", energia/(2.0*(L*L)));
            fprintf(magnet, "%g\n", magnetizacion/((double) (L*L)));
            //printf("Contador: %d\n", contador);
            contador++;
        }
        beta+=paso_beta;
    }
    /*
    beta=0.5;

    for (int j=0; j<Nbetas; j++){
        printf("Beta = %lf\n", beta);
        t_f=time(NULL);
        printf("Tiempo transcurrido desde inicio: %d segundos\n", (int)(t_f-t_0));

        for (int i=0; i<5; i++){//Calculamos los posibles valores de nuestra función, así evitamos muchos cálculos en el algoritmo de metrópolis
            prob[i] = funcion((-8+i*4), beta);
            //printf("Prob %d = %lf\n", i, prob[i]);
        }
        for (int k=0; k<Nconfig; k++){
            inicializa_config(config);
            //for (int huh=0; huh<(L*L); huh++) printf("%d\t", config[huh]);
            //printf("Config inicializada: \n");
            /*for (int index=0; index<(L*L); index++){
                printf("%d\t ", config[index]);
            }

            energia = calcula_energia(config, xp, yp);
            magnetizacion = calcula_magnetizacion(config);
            //printf("Energia inicial: %lf   magnetizacion inicial: %lf\n", energia/(double)(2.0*L*L),magnetizacion/((double) (L*L)));
            int aleatorio = (int)(10000*(rand()/(double)(RAND_MAX+1)));
           // printf("Aleatorio: %d\n", aleatorio);

            for (int i=0; i<(Nter+aleatorio); i++) {
                metropolis((L*L), config, xp, yp, xn, yn, prob, &energia, &magnetizacion);
            }
            //printf("Acabado metropolis\n");
            for (int index=0; index<(L*L); index++){
                if (index==(L*L)-1) {
                    fprintf(conf, "%d", config[index]);
                }
                else{
                    fprintf(conf, "%d\t ", config[index]);
                }
                //printf("%d\t ", config[index]);
            //printf("Index: %d\n", index);

            }
            fprintf(conf, "\n");
            //printf("Escrito en el fichero\n");
            int ord_dis;
            if (beta<0.44) ord_dis=1; //beta peq, temp alta = desordenado (dis)
            else if (beta>0.44) ord_dis=0; //temp baja = ordenado (ord)
            //printf("Ord_dis: %d\n", ord_dis);
            fprintf(labels, "%d\n ", ord_dis);
            fprintf(temp, "%lf\n", beta);
            fprintf(en_mag, "%d %g %g\n", contador, energia/(2.0*(L*L)), magnetizacion/((double) (L*L)));
            //printf("Contador: %d\n", contador);
            contador++;
        }
        beta+=paso_beta;
    }
    */



    fclose(conf);
    fclose(temp);
    fclose(labels);
    fclose(energy);
    fclose(magnet);

return 0;
}
