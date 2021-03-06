#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

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

void imprimeColuna(double A[][MAXM], int j, int start, int end);

double leArquivo(char nomearq[]);

void reescalaMatriz(double max);

void calculaNormas(double normas[]);

int maiorNorma(double normas[], int k);

void recalculaNormas(double normas[], int k);

void permutacao(int j, int k, double normas[]);

double decomposicaoQR(int k, double normas[]);

void multiplicaQt(int k, double gama);

void back_subst(double x[MAXN], int k);

/* ------------ Main -------------- */

int main(int argc, char* argv[]) {

	int i, j, k, c;
    int p[MAXM];
	double max, gama, aux;
	double normas[MAXM], x[MAXM];
    clock_t inicio, fim;

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

    inicio = clock();

    for (k = 0; k < m; k++) {

        if (verbose) printf("----- Q%d -----\n", k);

        recalculaNormas(normas, k);
        c = maiorNorma(normas, k);
        p[k] = c;

        if (normas[c] <= EPS)
            break;

        permutacao(c, k, normas);

        gama = decomposicaoQR(k, normas);

        multiplicaQt(k, gama);

        if (verbose) printf("--------------\n");
    }

    printf("Fim do QR. Posto(A) = %d\n", k);

    fim = clock();
    printf("%lf mili-segundos\n", ((double)(fim- inicio) / ((double)CLOCKS_PER_SEC))*1000);
    
    /* Cálculo de x */

    back_subst(x, k);

    /* Permuta x */
    for (j = k - 1; j >= 0; j--)
        if (j != p[j]) {
            if (verbose) printf("Trocando as colunas %d e %d de x...\n", j, p[j]);
            aux = x[j];
            x[j] = x[p[j]];
            x[p[j]] = aux;           
        }
    printf("x: \n");
    imprimeVetor(x, k);

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

	int i;

	for (i = 0; i < l; i++)
		printf("%lf ", v[i]);
	printf("\n\n");
}

void imprimeColuna(double A[][MAXM], int j, int start, int end) {

    int i;

    for (i = start; i < end; i++)
        printf("%lf ", A[i][j]);
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

    fscanf(arq, "%d ", &n);
    fscanf(arq, "%d ", &m);
    
    for (i = 0; i < n * m; i++) {
        fscanf(arq, "%d", &l);
        fscanf(arq, "%d", &c);
        fscanf(arq, "%lf", &valor);
        A[l][c] = valor;
        if (fabs(valor) > max) max = fabs(valor);
    }

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

        printf("b reescalado: \n");
        imprimeVetor(b, n);
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

int maiorNorma(double normas[], int k) {

    int j, c = -1;
    double max = -1;

    for (j = k; j < m; j++)
        if (normas[j] > max){
            c = j;
            max = normas[j];
        }

    if (verbose && c != -1)
        printf("Maior norma: %lf (normas[%d])\n", max, c);

    return c;


}

void recalculaNormas(double normas[], int k) {

	int j;
	if (k > 0){
		for (j = k; j < m; j++)
			normas[j] -= A[k - 1][j] * A[k - 1][j];

		if (verbose && k < m) {
			printf("Vetor de normas recalculadas: \n");
			imprimeVetor(normas, m);
		}
	}
}

void permutacao(int j, int k, double normas[]) {

    double aux, auxNorma;
    int i;

    if (j != k) {
        if (verbose) {
            printf("Trocando as colunas %d e %d...\nAntes:\n", j, k);
            imprimeMatriz(A, n, m);
        }

        for (i = 0; i < n; i++) {
            aux = A[i][j];
            A[i][j] = A[i][k];
            A[i][k] = aux;
        }

        auxNorma = normas[j];
        normas[j] = normas[k];
        normas[k] = auxNorma;

        if (verbose) {
            printf("Depois:\n");
            imprimeMatriz(A, n, m);
        }
    }
}

/* Essa função calcula tau, gama e o vetor u (será guardado em A) tal que Q = I - gama * u * uT é um refletor 
   que zera tudo abaixo do primeiro elemento da coluna. Como o primeiro elemento de u é sempre 1, -tau é guar-
   dado em seu lugar.                                                                                         */
double decomposicaoQR(int k, double normas[]) {

    int j;
    double tau, gama;

    j = k;
    tau = sqrt(normas[j]);
    if (A[j][j] < 0) tau *= -1;
    A[j][j] += tau;
    if (fabs(tau) > EPS)
        gama = A[j][j] / tau;
    else gama = 0;

    for (k = j + 1; k < n; k++)
        A[k][j] /= A[j][j];

    A[j][j] = -tau;
    if (verbose && gama != 0) {
        printf("Tau: %lf\n", tau);
        printf("Gama: %lf\n", gama);
        printf("u: \n");
        printf("1 ");
        imprimeColuna(A, j, j + 1, n);
    }

    return gama;
}
/* Essa função multiplica o vetor b por um Qt = I - gama * uT * u */
void multiplicaQt(int k, double gama) {

    int i, j;
    double c, xk;
    double vT[MAXN], temp[MAXN];

    /* ---- Qt * b ---- */
    /* c = gama * uT * b */
    c = gama * b[k];

    for (i = k + 1; i < n; i++) 
        c += gama * A[i][k] * b[i]; /* A[i][k] é u[i] */

    /* b = b - u * c */
    b[k] -= c; 
    for (i = k + 1; i < n; i++)
        b[i] -= A[i][k] * c;

    /* ---- Qt * A ---- */ 

    /* vT = gama * uT Obs: vT existe em no intervalo [k, n[ */
    vT[k] = gama;
    temp[k] = 0;
    for (i = k + 1; i < n; i++){
        vT[i] = gama * A[i][k];
        temp[i] = 0;
    }

    /* temp = vT * A */

    /* primeira linha de A (pois não precisamos mexer em A[k][k]) */

    temp[k] = vT[k] * (1 - gama) * A[k][k];
    for (j = k + 1; j < m; j++)
        temp[j] += vT[k] * A[k][j];

    /* Resto das linhas de A */

    for (i = k + 1; i < n ; i++) 
        for (j = k + 1; j < m; j++)
            temp[j] += vT[i] * A[i][j];  

    /* A = A - u * temp */
    /* Primeira linha de A (pois não precisamos mexer em A[k][k])) */
    for (j = k + 1; j < m; j++)
        A[k][j] -= temp[j];

    /* Resto das linhas de A */

    for (i = k + 1; i < n; i++) 
        for (j = k + 1; j < m; j++)
            A[i][j] -= A[i][k] * temp[j];


    if (verbose) {
        printf("Q^t * b: \n");
        imprimeVetor(b, n);

        printf("A modificado (R + vetores u): \n");
        imprimeMatriz(A, n, m);
    }

}

void back_subst(double x[], int k) {

    int i, j;      

    
    for (i = k-1; i >= 0; i--) {
        x[i] = b[i];
        for (j = k-1; j > i; j--){
            x[i] -= A[i][j] * x[j];
        }
        x[i] = x[i] / A[i][i]; 
    }

    if (verbose) {
        printf("x:\n");
        imprimeVetor(x, n);
    }   
}
