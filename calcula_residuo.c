#include <stdio.h>
#include <math.h>

#define EPS 0.0000000000000001
#define MAXN 1000
#define MAXM 1000

double A[MAXN][MAXM];
double b[MAXN];
int n = -1;
int m = -1;

double calcula_Residuo(double x_aprox[MAXM]);
void leArquivo(char nomearq[]);
void le_X(char nomearq[], double x_aprox[MAXM]);
void imprimeMatriz(double M[][MAXM], int l, int c);
void imprimeVetor(double v[], int l);

/*Usa um x aproximado estimado a partir de um polinômio, que pode ter
grau menor ou igual ao original e verifica o resíduo. O intuito é checar
se a aproximação obtida foi boa.*/

double calcula_Residuo(double x_aprox[MAXM]) {
	int i, j;
	double residuo, temp[MAXN];
	residuo = 0;

	for (i = 0; i < n; i++)
		temp[i] = 0;
		for (j = 0; j < m; j++)
			temp[i] += A[i][j] * x_aprox[j]; 

	for (i = 0; i < n; i++)
		temp[i] = b[i] - temp[i];

	for (i = 0; i < n; i++)
		residuo += temp[i] * temp[i];

	return sqrt(residuo);
}

void leArquivo(char nomearq[]) {

	FILE *arq;
	int i, j, l, c;
	double valor, max = 0;

	arq = fopen(nomearq, "r");
    if (arq == NULL) {
        printf("Arquivo não encontrado\n");
        return;
    }

    fscanf(arq, "%d ", &n);
    fscanf(arq, "%d ", &m);
    
    for (i = 0; i < n * m; i++) {
        fscanf(arq, "%d", &l);
        fscanf(arq, "%d", &c);
        fscanf(arq, "%lf", &valor);
        A[l][c] = valor;
    }

    for (i = 0; i < n; i++) {
        fscanf(arq, "%d", &l);
        fscanf(arq, "%lf", &valor);
        b[l] = valor;
    }

    fclose(arq);

    printf("A: \n");
	imprimeMatriz(A, n, m);

	printf("b: \n");
	imprimeVetor(b, n);
}

void le_X(char nomearq[], double x_aprox[]) {

	FILE *arq;
	int i, l;
	double valor, max = 0;

	arq = fopen(nomearq, "r");
    if (arq == NULL) {
        printf("Arquivo não encontrado\n");
        return;
    }

    fscanf(arq, "%d ", &m);
    
    for (i = 0; i < m; i++) {
        fscanf(arq, "%d", &l);
        fscanf(arq, "%lf", &valor);
        b[l] = valor;
    }

    fclose(arq);

	printf("x: \n");
	imprimeVetor(x_aprox, m);
}

void imprimeMatriz(double M[][MAXM], int l, int c) {

	int i, j;

	for (i = 0; i < l; i++) {
		for (j = 0; j < c; j++)
			printf("%lf ", M[i][j]);
		printf("\n");
	}
	printf("\n");
}

void imprimeVetor(double v[], int l) {

	int i;

	for (i = 0; i < l; i++)
		printf("%lf ", v[i]);
	printf("\n\n");
}

int main(int argc, char* argv[]) {

	int i, j, k;
	double residuo, x_aprox[MAXM];

	if (argc <= 2) {
        printf("\nModo de usar:\n"); 
        printf("1º argumento: o nome do arquivo contendo A e b\n");
        printf("2º argumento: o nome do arquivo contendo x aproximado\n");
        return 0;
    }

    leArquivo(argv[1]);
    le_X(argv[2], x_aprox);
    residuo = calcula_Residuo(x_aprox);

    printf("Resíduo = %lf\n", residuo);

    return 0;
}