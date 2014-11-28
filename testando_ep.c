#include <stdio.h>
#include <math.h>

#define eps 0.0000000000000001
#define nmax 1000
#define mmax 1000

void back_subst(double A[nmax][mmax], double b[nmax], double x[nmax], int k);

void back_subst(double A[][mmax], double b[], double x[], int k) {

	int i, j;      

    for (i = k-1; i >= 0; i--) {
        for (j = k-1; j >= i; j--)
            x[i] = b[i] - A[i][j] * x[j];
        x[i] = x[i] / A[i][i];  
    }

    for (i = 0; i < k; i++)
				printf("x[%d] = %lf\n", i, x[i]);
}

int main(int argc, char* argv[]) {

	double A[nmax][mmax], b[nmax], vT[nmax], temp[nmax];
	double gama, v[nmax], x[nmax], col_square[2];
	double max, aux, valor, tau;
	int i, j, n, m, k, l, c, p[nmax];
	FILE *arq;

	max = -10000000; 

	if (argc < 2) {
        printf("\nModo de usar:\n"); 
        printf("1º argumento: o nome do arquivo de entrada\n");
        return 0;
    }

	arq = fopen(argv[1], "r");
    if (arq == NULL) {
        printf("Arquivo não encontrado\n");
        return 0;
    }

    fscanf(arq, "%d", &n);
    fscanf(arq, "%d", &m);
    
    for (i = 0; i < n*m; i++) {
        fscanf(arq, "%d", &l);
        fscanf(arq, "%d", &c);
        fscanf(arq, "%lf", &valor);
        A[l][c] = valor;
        printf("l = %d, valor = %lf\n", valor);
    }

    fscanf(arq, "%d", &n);
    for (i = 0; i < n; i++) {
        fscanf(arq, "%d", &l);
        fscanf(arq, "%lf", &valor);
        b[l] = valor;
        printf("l = %d, valor = %lf\n", valor);
    }

    printf("Começo do b:\n");
    for (i = 0; i < n; i++) 
    	printf("b[%d] = %lf\n", i, b[i]);
/*    for (i = 0; i < n; i++) {
    	for (j = 0; j < m; j++)
    		printf("A[%d][%d] = %lf\n", i, j, A[i][j]);
    	printf("b[%d] = %lf\n", i, b[i]);
    }*/

    fclose(arq);

    /*Reescalando a matriz*/
	for (i = 0; i < n; i++)	{	/*Encontro o maior elemento da matriz*/
		for (j = 0; j < m; j++)		
	    	if (fabs(A[i][j]) > max) /*Essa parte é importante para evitar que ocorram overflows ou underflows*/
	        	max = A[i][j];	
	}

	if (max <= eps) {
		printf("A matriz é nula\n");
		return 0;
	}

	for (i = 0; i < n; i++) {	
		for (j = 0; j < m; j++)
			A[i][j] /= max;
		b[i] /= max;
	}

	printf("A:\n");
	for (i =0; i < n; i++) {
		for (j = 0; j<m; j++)
			printf("%lf ", i, j, A[i][j]);
		printf("\n");
	}	


	/*for (i = 0; i < n; i++) {	
		for (j = 0; j < m; j++)
			printf("A[%d][%d] = %lf\n", i, j, A[i][j]);
		printf("b[%d] = %lf\n", i, b[i]);
	}*/

	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			v[j] += A[i][j] * A[i][j];


	for (i =0; i < n; i++)
		printf("v[%d] = %lf\n", i, v[i]);

	for (k = 0; k <= m; k++) {
		
		/*recalculando as normas ao quadrado de cada coluna*/
		if (k > 0 && k < m) {
			for (j = k; j < m; j++)			
				v[j] = v[j] - A[k-1][j]*A[k-1][j];
		}

		col_square[0] = 0;
		col_square[1] = -1;

		for (j = k; j < m; j++) {
			if (v[j] > col_square[0])
				col_square[0] = v[j];
				col_square[1] = j;
		}

		p[k] = col_square[1]; /*vetor da permutação*/
		printf("p[%d] = %d\n", k, p[k]);
		printf("col_square[0] = %lf\n", col_square[0]);

		if (col_square[0] <= eps || k == m) {
			gama = 0;
			printf("n = %d\n", n);
			printf("Acabou o algoritmo QR. O posto da matriz é: %d.\n", k);
			for (i = 0; i < n; i++)
				for (j = 0; j < m; j++)
					printf("A[%d][%d] = %lf\n", i, j, A[i][j]);
			
			back_subst(A, b, x, k);

			for (j = n-1; j >= 0; j--) {
				aux = b[p[j]];
				b[p[j]] = b[j];
				b[j] = aux;
			}

			for (i = 0; i < n; i++)
				printf("b[%d] = %lf\n", i, b[i]);

			

			return 0;
		}
		
		if (col_square[1] != k) {			/*Feito o swap das colunas*/
			j = col_square[1];
			for (i = 0; i < n; i++) {
				aux = A[i][j];
				A[i][j] = A[i][k];
				A[i][k] = aux;
			}

			/*Falta permutar o b*/
			aux = v[j];
			v[j] = v[k];
			v[k] = aux;
		}

		/*for (i = 0; i < n; i++) {	
			for (j = 0; j < m; j++) {
				printf("A[%d][%d] = %lf\n", i, j, A[i][j]);
				printf("v[%d] = %lf\n", i, v[i]);
			}
		}*/
		tau = sqrt(v[k]);
	    if (A[k][k] < 0)
		   tau = -tau;
	
		A[k][k] += tau;
		gama = A[k][k]/tau;
		if (gama <= eps) {
			printf("A é singular\n");
			return 0;
		}
		printf("tau = %lf e gama = %lf\n", tau, gama);

		/*Checar se começa do k*/
		for (i = k+1; i < n; i++)   /*Desse modo, os vetores u's estão armazenados em A aonde seriam os zeros*/
			A[i][k] /= A[k][k];

		A[k][k] = tau;

		/*Aplicando o refletor Qk à matriz A*/
		vT[k] = gama;
		for (i = k+1; i < n; i++)
			vT[i] = gama*A[i][k];

		for (i = 0; i < n; i++)
			printf("vt[%d] = %lf\n", i, vT[i]);

		for (i = k; i < n; i++) {
			temp[i] = vT[i];
			vT[i] = 0;
		}

		printf("temp\n");
		for (i = 0; i < n; i++)
			printf("%lf\n", temp[i]);
				 
		  /*Checar se realmente começa do k ou k+1*/
		for (j = k+1; j < m; j++)
			for (i = k; i < n; i++)
				vT[j-1] += temp[i] * A[i][j];

		printf("vt\n");
		for (i = 0; i < n; i++)
			printf("%lf\n", vT[i]);
			/*Aqui obtenho vT que representa a matriz 1 X m-1, 1 X m-2 ... de acordo com a iteração*/

		printf("A antes:\n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++)
				printf(" %lf ", A[i][j]);	
			printf("\n");
		}
		
		for (i = k; i < n; i++)
			for (j = k+1; j < m; j++)	
				A[i][j] = A[i][j] - A[i][k]*vT[i];

		for (i = 0; i < n; i++)
			for (j = 0; j < m; j++)
				printf("A[%d][%d] = %lf\n", i, j, A[i][j]);		
/*-----------------------------------------------------------------------*/

		vT[k] = gama;
		for (i = k+1; i < n; i++)
			vT[i] = gama*A[i][k]; /*Lembrando-se de que A[i][k] contém u*/

		/*for (j = k; j < n; j++) {
			temp[j] = vT[j];
			vT[j] = 0;
		}*/

		aux = 0;
		for (i = k+1; i < n; i++)
			aux += A[i][k] * b[i];

		for (j = k; j < n; j++)
			vT[j] *=  aux;

		for (i = k; i < n; i++)
			b[i] = b[i] - vT[i];

		for (i = k; i < n; i++)
			printf("b[%d] = %lf\n", i, b[i]);		
	}

	return 0;
}