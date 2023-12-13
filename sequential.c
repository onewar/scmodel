#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define A_1 -4.0
#define A_2 -1.0
#define B_1 4.0
#define B_2 5.0
#define M 40
#define N 40
#define delta 1e-6

double h_1 = (B_1 - A_1) / M;
double h_2 = (B_2 - A_2) / N;
double eps;

typedef struct
{
	size_t height, width;
	double *data;
} Matrix;

typedef struct
{
	size_t height, width;
	int *data;
} IMatrix;

#define index(_m, _x, _y) (((_m)->data)[(_m)->width * (_y) + (_x)])

#define matrix_for_each_start(m)                \
	for (size_t i = 0; i < (m)->height; i++)    \
	{                                           \
		for (size_t j = 0; j < (m)->width; j++) \
		{                                       \
			double *elem = &index(m, i, j);

#define matrix_for_each_end \
	}                       \
	}

#define imatrix_for_each_start(m)               \
	for (size_t i = 0; i < (m)->height; i++)    \
	{                                           \
		for (size_t j = 0; j < (m)->width; j++) \
		{                                       \
			int *elem = &index(m, i, j);

Matrix *w;

/*
 *	Always return successfully, no need for assert :P
 */
Matrix *
matrix_new(size_t height, size_t width)
{
	Matrix *m = (Matrix *)malloc(sizeof(Matrix));
	assert(m);

	size_t size = width * height;
	*m = (Matrix){
		.height = height,
		.width = width,
		.data = (double *)malloc(sizeof(double) * size),
	};

	matrix_for_each_start(m)
		*elem = 0;
	matrix_for_each_end

		assert(m->data);

	return m;
}

void matrix_free(Matrix *m)
{
	free(m->data);
	free(m);
	return;
}

IMatrix *
imatrix_new(size_t height, size_t width)
{
	IMatrix *m = (IMatrix *)malloc(sizeof(IMatrix));
	assert(m);

	size_t size = width * height;
	*m = (IMatrix){
		.height = height,
		.width = width,
		.data = (int *)malloc(sizeof(int) * size),
	};

	imatrix_for_each_start(m)
		*elem = 0;
	matrix_for_each_end

		return m;
}

void imatrix_free(IMatrix *m)
{
	free(m->data);
	free(m);
	return;
}

void init(void)
{
	w = matrix_new(M + 1, N + 1);
	eps = pow(h_1 > h_2 ? h_1 : h_2, 2);
	return;
}

double
F(double x, double y)
{
	return (((x >= -3) || (x <= 3)) && y >= 0 &&
			((4 * x + 3 * y - 12 <= 0) || (4 * x - 3 * y + 12 >= 0)))
			   ? 1
			   : 0;
}

double
funVScalProd(Matrix *u, Matrix *v, double h_1, double h_2)
{
	double res = 0;
	for (int i = 1; i < M; i++)
	{
		for (int j = 1; j < N; j++)
			res += h_1 * h_2 * index(u, i, j) * index(v, i, j);
	}
	return res;
}

double
funVNorm(Matrix *u, double h_1, double h_2)
{
	return sqrt(funVScalProd(u, u, h_1, h_2));
}

Matrix *
funVSubtract(Matrix *u, Matrix *v)
{
	Matrix *res = matrix_new(M + 1, N + 1);
	for (int i = 1; i < M; i++)
	{
		for (int j = 1; j < N; j++)
			index(res, i, j) = index(u, i, j) - index(v, i, j);
	}
	return res;
}

Matrix *
funVmultConst(Matrix *u, double c)
{
	Matrix *res = matrix_new(M + 1, N + 1);
	matrix_for_each_start(u)
		index(res, i, j) = *elem * c;
	matrix_for_each_end return res;
}

void nodeTypeDef(IMatrix *nodeType)
{
	imatrix_for_each_start(nodeType) if (F(A_1 + (j + 0.5) * h_1, A_2 + (i + 0.5) * h_2) > 0)
		*elem = 1;
	matrix_for_each_end
}

void FijDef(Matrix *Fij, IMatrix *nodeType)
{
	for (int i = 1; i < M; i++)
	{
		for (int j = 1; j < N; j++)
		{
			int *elem = &index(nodeType, i, j);
			double yU = A_2 + (i + 0.5) * h_2;
			double yD = A_2 + (i - 0.5) * h_2;
			double xR = A_1 + (j + 0.5) * h_1;
			double xL = A_1 + (j - 0.5) * h_1;

			double *res = &index(Fij, i, j);

			if (xR <= -3 || yU <= 0 || yD >= 4 || xL >= 3 ||
				4 * xL + 3 * yD - 12 >= 0 ||
				4 * xR - 3 * yD + 12 <= 0)
			{
				continue;
			}
			else if (xR > -3 && xL < -3)
			{
				/*	1 - 4	*/
				int left = index(nodeType, i - 1, j);
				switch (*elem * 2 + left)
				{
				case 2: // 1
					*res = 1 * (((-3. / 4 * yU + 3 + xR) +
								 (3 + xR)) *
								(yU / 2) / (h_1 * h_2));
					break;
				case 0: // 2
					*res = 1 * (3 + xR) * (4. / 3 * xR + 4) / 2 /
						   (h_1 * h_2);
					break;
				case 3: // 3
					*res = 1 * (((-3. / 4 * yU + 3 + xR) + (-3. / 4 * yD + 3 + xR)) * (h_2 / 2)) /
						   (h_1 * h_2);
					break;
				case 1: // 4
					*res = 1 * (-3. / 4 * yD + 3 + xR) *
						   (4. / 3 * xR + 4 - yD) / 2 /
						   (h_1 * h_2);
					break;
				}
			}
			else if (xR > 3 && xL < 3)
			{
				/*	5 - 8	*/
				if (index(nodeType, i, j - 1) &&
					!index(nodeType, i - 1, j - 1))
				{
					// 5
					*res = 1 * ((3. / 4 * yU - 3 - xL) + (3 - xL)) * yU / 2 / (h_1 * h_2);
				}
				else if (!index(nodeType, i, j - 1) &&
						 !index(nodeType, i - 1, j - 1))
				{
					// 6
					*res = 1 * ((3 - xL) * (-4. / 3 * xL + 4) / 2) /
						   (h_1 * h_2);
				}
				else if (index(nodeType, i, j - 1) &&
						 index(nodeType, i - 1, j - 1))
				{
					// 7
					*res = 1 * (((-3. / 4 * yU + 3 - xL) + (-3. / 4 * yD + 3 - xL)) * h_2 / 2) /
						   (h_1 * h_2);
				}
				else if (!index(nodeType, i, j - 1) &&
						 index(nodeType, i - 1, j - 1))
				{
					// 8
					*res = 1 * (((-4. / 3 * xL + 4 - yD) * (-3. / 4 * yD + 3 - xL)) / 2) /
						   (h_1 * h_2);
				}
			}
			else if (xL < 0 && xR > 0)
			{
				/*	9 - 36	*/
				if (yD < 0)
				{
					// 9
					if (*elem && index(nodeType, i, j - 1))
						*res = yU / h_2;
				}
				else if (yD > 0)
				{
					// 10
					if (*elem && index(nodeType, i, j - 1))
						*res = 1;
				}
				else if (yD > 0 && yU < 4)
				{
					// 11
					if (!*elem && index(nodeType, i - 1, j) &&
						index(nodeType, i - 1, j - 1) &&
						!index(nodeType, i, j - 1))
					{
						*res = 1 * (h_1 * h_2 - ((xR + 3. / 4 * yU - 3) * (yU + 4. / 3 * xR - 4) / 2) - (-xL + 3. / 4 * yU - 3) * (yU - 4. / 3 * xL + 4) / 2) /
							   (h_1 * h_2);
					}
				}
				else if (yD > 0 && yU > 4) // 12
				{
					if (!*elem && index(nodeType, i - 1, j) &&
						index(nodeType, i - 1, j - 1) &&
						!index(nodeType, i, j - 1))
					{
						*res = 1 * ((((-4. / 3 * xR + 4 - yD) + (4 - yD)) * xR / 2) + ((4. / 3 * xL + (4 - yD)) * (-xL / 2)) / (h_1 * h_2));
					}
				}
				else if (yU < 4) // 13-27
				{
					// 13-16
					if (yD > 0 && 4 * xL - 3 * yU + 12 < 0 && 4 * xL - 3 * yD + 12 < 0) // 13
					{
						*res = 1 * (((-3. / 4 * yD + 3 + (-3. / 4 * yU + 3)) * h_2) / (h_1 * h_2));
					}
					else if (yD > 0 && 4 * xL - 3 * yU + 12 < 0 && 4 * xL - 3 * yD + 12 > 0) // 14
					{
						*res = 1 * (((-3. / 4 * yD + 3 + (-3. / 4 * yU + 3)) * h_2 / 2) + ((-xL * h_2) - ((-xL + 3. / 4 * yU - 3) * (yU - 4. / 3 * xL - 4) / 2)) / (h_1 * h_2));
					}
					else if (yD < 0 && 4 * xL - 3 * yU + 12 < 0 && xL < -3) // 15
					{
						*res = 1 * ((3 + (-3. / 4 * yU + 3)) * yU) / (h_1 * h_2);
					}
					else if (yD < 0 && 4 * xL - 3 * yU + 12 < 0 && xL > -3) // 16
					{
						*res = 1 * ((3 + (-3. / 4 * yU + 3) * yU / 2) + (-xL * yU) - ((-xL + 3. / 4 * yU - 3) * (yU - 4. / 3 * xL - 4) / 2) / (h_1 * h_2));
					}
					// 17-18
					else if (4 * xL - 3 * yD + 12 < 0 && xR < 3 && 4 * xR + 3 * yU - 12 > 0)
					{
						if (yD > 0) // 17
						{
							*res = 1 * (((-3. / 4 * yD + 3) + (-3. / 4 * yU + 3)) * h_2 / 2 + (xR * h_2) - (yU + 4. / 3 * xR - 4) * (xR + 3. / 4 * yU - 3) / 2) / (h_1 * h_2);
						}
						else if (yD < 0) // 18
						{
							*res = 1 * ((-xL + 3 - 3. / 4 * yU) * yU - (3 - xR * (-4. / 3 * xR + 4)) / 2) / (h_1 * h_2);
						}
					}
					// 19
					else if (yD < 0 && xL > -3 && xR < 3 && 4 * xL - 3 * yU + 12 < 0 && 4 * xR + 3 * yU - 12 > 0)
					{
						*res = (h_1 * yU - (xR + 3. / 4 * yU - 3) * (yU + 4. / 3 * xR - 4) / 2 - (-xL + 3. / 4 * xL - 3) * (yU - 3. / 4 * xL + 4) / 2) / (h_1 * h_2);
					}
					else if (4 * xR + 3 * yU - 12 < 0 && 4 * xL - 3 * yU + 12 < 0) // 20-23
					{
						if (yD < 0 && xL > -3) // 20
						{
							*res = (xR * yU + (3 + 3 - 3. / 4 * yU) * yU / 2) / h_1 * h_2;
						}
						else if (yD < 0 && xL > -3) // 21
						{
							*res = (yU * h_1 - (yU - 4. / 3 * xL - 4) * (-xL + 3. / 4 * yU - 3) / 2) / (h_1 * h_2);
						}
						else if (yD > 0 && 4 * xL - 3 * yD + 12 < 0) // 22
						{
							*res = 1 * ((xR * h_2 + ((-3. / 4 * yU + 3) + (-3. / 4 * yD + 3)) * h_2 / 2) / (h_1 * h_2));
						}
						else if (yD > 0 && 4 * xL - 3 * yD + 12 > 0) // 23
						{
							*res = (h_1 * h_2 - (yU - 4. / 3 * xL - 4) * (-xL + 3. / 4 * yU - 3) / 2) / (h_1 * h_2);
						}
					}
					else if (4 * xL - 3 * yU + 12 > 0) // 24-27
					{
						if (yD < 0 && xR > 3) // 24
						{
							*res = ((-xL * yU) + (3 - 3. / 4 * yU + 3) * yU / 2) / (h_1 * h_2);
						}
						else if (yD < 0 && xR > 3) // 25
						{
							*res = (h_1 * yU - (yU + 4. / 3 * xR - 4) * (xR + 3. / 4 * yU - 3) / 2) / (h_1 * h_2);
						}
						else if (yD > 0 && 4 * xR + 3 * yD - 12 > 0) // 26
						{
							*res = (-xL * h_2 + (-3. / 4 * yD + 3 + (-3. / 4 * yU + 3)) * (h_2 / 2)) / (h_2 * h_1);
						}
						else if (yD > 0 && 4 * xR + 3 * yD - 12 < 0) // 27
						{
							*res = (h_1 * h_2 - (yU + 4. / 3 * xR - 4) * (xR + 3. / 4 * yU - 3) / 2) / (h_1 * h_2);
						}
					}
				}
				else if (yU > 4) // 28-35
				{
					if (yD > 0) // 28-31
					{
						if (4 * xR + 3 * yD - 12 > 0 && 4 * xL - 3 * yD + 12 < 0) // 28
						{
							*res = (4 - yD) * (-3. / 4 * yD + 3) / (h_1 * h_2);
						}
						else if (4 * xR + 3 * yD - 12 > 0 && 4 * xL - 3 * yD + 12 > 0) // 29
						{
							*res = (4 - yD) * (-3. / 4 * yD + 3) / 2 + ((4. / 3 * xL + 4 - yD) + (4 - yD)) * (-xL / 2) / (h_1 * h_2);
						}
						else if (4 * xR + 3 * yD - 12 < 0 && 4 * xL - 3 * yD + 12 < 0) // 30
						{
							*res = (4 - yD) * (-3. / 4 * yD + 3) / 2 + (-4. / 3 * xR + 4 - yD) + (4 - yD) * (xR / 2) / (h_1 * h_2);
						}
						else if (4 * xR + 3 * yD - 12 < 0 && 4 * xL - 3 * yD + 12 > 0) // 31
						{
							*res = ((-4. / 3 * xR + 4 - yD) + (4 - yD) * (xR / 2) + ((4. / 3 * xL + 4 - yD) + (4 - yD)) * (-xL / 2)) / (h_1 * h_2);
						}
					}
					else if (yD < 0) // 32-35
					{
						if (4 * xR + 3 * yD - 12 > 0 && xL < -3) // 32
						{
							*res = 12 / (h_1 * h_2);
						}
						else if (4 * xR + 3 * yD - 12 > 0 && xL > -3) // 33
						{
							*res = (12 - (3 + xL) * (4. / 3 * xL + 4) / 2) / (h_1 * h_2);
						}
						else if (4 * xL - 3 * yD + 12 < 0 && xR > 3) // 34
						{
							*res = (12 - (3 - xR) * (-4. / 3 * xR + 4) / 2) / (h_1 * h_2);
						}
						else if (xR < 3 && xL > -3) // 35
						{
							*res = (12 - (3 - xR) * (-4. / 3 * xR + 4) / 2 - (3 + xL) * (4. / 3 * xL + 4) / 2) / (h_1 * h_2);
						}
					}
				}
			}
		}
	}
}

double
coefA(int i, int j, IMatrix *nodeType)
{
	double x = A_1 + (j - 0.5) * h_1;
	double yU = A_2 + (i + 0.5) * h_2;
	double yD = A_2 + (i - 0.5) * h_2;
	int U = index(nodeType, i, j - 1);
	int D = index(nodeType, i - 1, j - 1);

	if (!(D + U))
	{
		if (x < -3 || x > 3 || yU < 0 || yD > 4 || 4 * x + 3 * yD - 12 > 0 || 4 * x - 3 * yD + 12 < 0)
		{
			return 1 / eps;
		}
		else if (x > 0 && x < 3)
		{
			return (-4. / 3 * x + 4) / h_2 + (1 - (-4. / 3 * x + 4) / h_2) / eps;
		}
		else if (x > -3 && x < 0)
		{
			return (4. / 3 * x + 4) / h_2 + (1 - (4. / 3 * x + 4) / h_2) / eps;
		}
	}
	else if (D + U == 1)
	{
		if (yD < 0)
		{
			return (yU) / h_2 + (1 - yU / h_2) / eps;
		}
		else
		{
			if (x > 0 && x < 3)
			{
				return (-4. / 3 * x + 4 - yD) / h_2 + (1 - (-4. / 3 * x + 4 - yD) / h_2) / eps;
			}
			else if (x > -3 && x < 0)
			{
				return (4. / 3 * x + 4 - yD) / h_2 + (1 - (4. / 3 * x + 4 - yD)) / h_2 / eps;
			}
		}
	}
	else
	{
		return 1;
	}
}

double
coefB(int i, int j, IMatrix *nodeType)
{
	double y = A_2 + (i - 0.5) * h_2;
	double xR = A_1 + (j + 0.5) * h_1;
	double xL = A_1 + (j - 0.5) * h_1;
	int R = index(nodeType, i - 1, j);
	int L = index(nodeType, i - 1, j - 1);
	if (L + R == 0)
	{
		if (y < 0 || y > 4 || xL > 3 || xR < -3 || 4 * xL + 3 * y - 12 > 0 || 4 * xR - 3 * y + 12 < 0)
		{
			return 1 / eps;
		}
		else
		{
			return (-3. / 2 * y + 6) / h_1 + (1 - (-3. / 2 * y + 6) / h_1) / eps;
		}
	}
	// if one of the two endpoints of the segment is outside d
	// else if (nodeType[i - 1][j - 1] + nodeType[i][j] == 1) {
	else if (L + R == 1)
	{
		if (xL < 3 && xL > -3)
		{
			return (-3. / 4 * y + 3 - xL) / h_1 + (1 - (-3. / 4 * y + 3 - xL) / h_1) / eps;
		}
		else if (xR > -3 && xR < 3)
		{
			return (-3. / 4 * y + 3 + xR) / h_1 + (1 - (-3. / 4 * y + 3 + xR) / h_1) / eps;
		}
	}
	else
	{
		return 1;
	}
}

Matrix *operatorA(Matrix *w, IMatrix *nodeType)
{
	Matrix *res = matrix_new(M + 1, N + 1);
	for (int i = 1; i < M; i++)
	{
		for (int j = 1; j < N; j++)
		{
			index(res, i, j) = -(coefA(i, j + 1, nodeType) * (index(w, i, j + 1) - index(w, i, j)) / h_1 - coefA(i, j, nodeType) * (index(w, i, j) - index(w, i, j - 1)) / h_1) / (h_1) - (coefB(i + 1, j, nodeType) * (index(w, i + 1, j) - index(w, i, j)) / h_2 - coefB(i, j, nodeType) * (index(w, i, j) - index(w, i - 1, j)) / h_2) / (h_2);
		}
	}
	return res;
}

int main()
{
	init();
	IMatrix *nodeType = imatrix_new(M, N);
	Matrix *Fij = matrix_new(M + 1, N + 1);
	nodeTypeDef(nodeType);
	FijDef(Fij, nodeType);
	int iterNum = 0;
	Matrix *wNewer = w;
	double norm;
	Matrix *residual = funVSubtract(operatorA(w, nodeType), Fij);
	Matrix *Ar = operatorA(residual, nodeType);
	double tau = funVScalProd(Ar, residual, h_1, h_2) / pow(funVNorm(Ar, h_1, h_2), 2);
	// Generalized minimal residual method
	clock_t start = clock();
	do
	{
		w = wNewer;
		residual = funVSubtract(operatorA(w, nodeType), Fij);
		Ar = operatorA(residual, nodeType);
		tau = funVScalProd(Ar, residual, h_1, h_2) / pow(funVNorm(Ar, h_1, h_2), 2);
		wNewer = funVSubtract(w, funVmultConst(residual, tau));
		iterNum++;
	} while (funVNorm(funVSubtract(wNewer, w), h_1, h_2) >= delta);
	clock_t end = clock();
	double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Execution succeed!\n");
	printf("Total iteration: %d\n", iterNum);
	printf("Time elapsed: %f ms\n", elapsed_time * 1000);

	// for (int i = M; i >= 0; i--)
	// {
	// 	for (int j = 0; j < N + 1; j++)
	// 	{
	// 		printf("%f ", index(wNewer, i, j));
	// 		if (j <= N - 1)
	// 		{
	// 			printf("%f, ", index(wNewer, M - i, j));
	// 		}
	// 		else
	// 		{
	// 			printf("%f\n", index(wNewer, M - i, j));
	// 		}
	// 	}
	// }
	return 0;
}