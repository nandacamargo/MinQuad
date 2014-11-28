#include <stdio.h>
#include <string.h>
#include <math.h>

#define EPS 0.0000000000000001
#define MAXN 1000
#define MAXM 1000
#define TRUE 1
#define FALSE 0

/* ------- Variáveis Globais ------ */

double A[MAXN][MAXM];
double b[MAXN];
int n = -1;
int m = -1;
int verbose = FALSE;

/* ----------- Funções ------------ */

void imprimeMatriz(double A[][MAXM], int n, int m);

void imprimeVetor(double v[], int n);

double leArquivo(char nomearq[]);

void reescalaMatriz(double max);

void calculaNormas(double normas[]);

void recalculaNormas(double normas[], int k);

void decomposicaoQR(int j, int k, double normas[]);

/* ------------ Main -------------- */

int main(int argc, char* argv[]) {

	int i, k;
	double max;
	double normas[MAXM];

	/* Lidando com os argumentos */

	if (argc < 2 || (argc >= 2 && strcmp(argv[1], "-v") == 0)) {
        printf("\nModo de usar:\n"); 
        printf("1º argumento: o nome do arquivo de entrada\n");
        printf("'-v' para ligar modo verbose\n");
        return 0;
    }

	i = 1;
    if (strcmp(argv[1], "-v") == 0) { 
    	i = 2;
    	verbose = TRUE;
    }

    else if (argc >= 3 && strcmp(argv[2], "-v") == 0)
    	verbose = TRUE;    	

    /* Leitura do Arquivo */

    max = leArquivo(argv[i]);
    if (max == -1) return 1;
    if (max <= EPS) {
    	printf("A matriz é nula!\n");
    	return 1;
    }

    /* Reescalando a matriz */

    reescalaMatriz(max);

    /* Construindo o vetor de normas */

    calculaNormas(normas);

    /* Decomposição QR */
    for (k = 0; k <= m; k++) {
    	recalculaNormas(normas, k);
    	decomposicaoQR(0, 0, normas);
    }

    return 0;
}

/* -- Implementação das Funções --- */

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

	int i, j;

	for (i = 0; i < l; i++)
		printf("%lf ", v[i]);
	printf("\n\n");
}


/* A função le e o arquivo e devolve o maior elemento (em módulo) da matriz */

double leArquivo(char nomearq[]) {

	FILE *arq;
	int i, j, l, c;
	double valor, max = 0;

	arq = fopen(nomearq, "r");
    if (arq == NULL) {
        printf("Arquivo não encontrado\n");
        return -1;
    }

    fscanf(arq, "%d", &n);
    fscanf(arq, "%d", &m);
    
    for (i = 0; i < n * m; i++) {
        fscanf(arq, "%d", &l);
        fscanf(arq, "%d", &c);
        fscanf(arq, "%lf", &valor);
        A[l][c] = valor;
        if (fabs(valor) > max) max = fabs(valor);
    }

    fscanf(arq, "%d", &n);
    for (i = 0; i < n; i++) {
        fscanf(arq, "%d", &l);
        fscanf(arq, "%lf", &valor);
        b[l] = valor;
    }

    fclose(arq);

    if (verbose) {
    	printf("A: \n");
		imprimeMatriz(A, n, m);

		printf("b: \n");
		imprimeVetor(b, n);
	}

    if (verbose) printf("Max: %lf\n", max);
    return max;

}


/* Essa função reescala a matriz e o vetor, encontrando o seu maior elemento (em módulo) e dividindo toda a matriz por ele. 
   Isso evita eventuais overflows e underflows que ocorreriam durante o processo. Os overflows não ocorrem porque 
   cada elemento a_ij da matriz é agora menor ou igual a 1 (em módulo). O underflow, por sua vez, só ocorreria no
   caso de existir algum a_ij muito menor do que o maior elemento da matriz (que agora será 1). Tais erros, no en-
   tanto, são extremamente pequenos, menores até do que a precisão do computador. Portanto, podemos ignorá-los. */

void reescalaMatriz(double max) {

	int i, j;

	for (i = 0; i < n; i++) {	
		for (j = 0; j < m; j++)
			A[i][j] /= max;
		b[i] /= max;
	}

	if (verbose) {
		printf("Matriz reescalada: \n");
		imprimeMatriz(A, n, m);
	}
}

/* Essa função calcula as normas das colunas da matriz A e devolve um vetor contendo as normas. */
void calculaNormas(double normas[]) {

	int i, j;

	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			normas[j] += A[i][j] * A[i][j];

	if (verbose) {
		printf("Vetor de normas:\n");
		imprimeVetor(normas, m);
	}
}

void recalculaNormas(double normas[], int k) {

	int j;
	if (k > 0){
		for (j = k; j < m; j++)
			normas[j] -= A[k - 1][j] * A[k - 1][j];

		if (k < m) {
			printf("Vetor de normas recalculadas: \n");
			imprimeVetor(normas, m);
		}
	}
}

/* Essa função calcula tau, gama e o vetor u (será guardado em A) tal que Q = I - gama * u * uT é um refletor 
   que zera tudo abaixo do primeiro elemento da coluna  */
void decomposicaoQR(int j, int k, double normas[]) {

	double tau, gama;

	tau = normas[j];
	if (A[k][k] < 0) tau *= -1;


}