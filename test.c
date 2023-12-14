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
#define INDEX(m, i, j) \
    ((m)[((i) * N) + (j)])
double h_1 = (B_1 - A_1) / M;
double h_2 = (B_2 - A_2) / N;
double eps;

double *w;
void init(void)
{
    w = (double *)malloc(sizeof(double) * ((M + 1) * (N + 1)));
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
funVScalProd(double *u, double *v, double h_1, double h_2)
{
    double res = 0;
    {
        for (int i = 1; i < M; i++)
        {
            for (int j = 1; j < N; j++)
                res += h_1 * h_2 * u[i * N + j] * v[i * N + j];
        }
        return res;
    }
}

double
funVNorm(double *u, double h_1, double h_2)
{
    return sqrt(funVScalProd(u, u, h_1, h_2));
}

double *
funVSubtract(double *u, double *v)
{
    double *res = (double *)malloc(sizeof(double) * ((M + 1) * (N + 1)));
    {

        for (int i = 1; i < M; i++)
        {
            for (int j = 1; j < N; j++)
                res[i * N + j] = u[i * N + j] - v[i * N + j];
        }
    }
    return res;
}

double *
funVmultConst(double *u, double c)
{
    double *res = (double *)malloc(sizeof(double) * ((M + 1) * (N + 1)));
    for (int i = 1; i < M; i++)
    {
        for (int j = 1; j < N; j++)
        {
            res[i * N + j] = u[i * N + j] * c;
        }
    }
    return res;
}

void nodeTypeDef(int *nodeType)
{
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (F(A_1 + (j + 0.5) * h_1, A_2 + (i + 0.5) * h_2) > 0)
            {
                nodeType[i * N + j] = 1;
                // cout << 1;
            }
        }
    }
}

void FijDef(double *Fij, int *nodeType)
{
    {
        for (int i = 1; i < M; i++)
        {
            for (int j = 1; j < N; j++)
            {
                int elem = nodeType[i * N + j];
                double yU = A_2 + (i + 0.5) * h_2;
                double yD = A_2 + (i - 0.5) * h_2;
                double xR = A_1 + (j + 0.5) * h_1;
                double xL = A_1 + (j - 0.5) * h_1;

                if (xR <= -3 || yU <= 0 || yD >= 4 || xL >= 3 ||
                    4 * xL + 3 * yD - 12 >= 0 ||
                    4 * xR - 3 * yD + 12 <= 0)
                {
                    continue;
                }
                else if (xR > -3 && xL < -3)
                {
                    /*	1 - 4	*/
                    int left = INDEX(nodeType, i - 1, j);
                    switch (elem * 2 + left)
                    {
                    case 2: // 1
                        Fij[i * N + j] = 1 * (((-3. / 4 * yU + 3 + xR) +
                                               (3 + xR)) *
                                              (yU / 2) / (h_1 * h_2));
                        break;
                    case 0: // 2
                        Fij[i * N + j] = 1 * (3 + xR) * (4. / 3 * xR + 4) / 2 /
                                         (h_1 * h_2);
                        break;
                    case 3: // 3
                        Fij[i * N + j] = 1 * (((-3. / 4 * yU + 3 + xR) + (-3. / 4 * yD + 3 + xR)) * (h_2 / 2)) /
                                         (h_1 * h_2);
                        break;
                    case 1: // 4
                        Fij[i * N + j] = 1 * (-3. / 4 * yD + 3 + xR) *
                                         (4. / 3 * xR + 4 - yD) / 2 /
                                         (h_1 * h_2);
                        break;
                    }
                }
                else if (xR > 3 && xL < 3)
                {
                    /*	5 - 8	*/
                    if (INDEX(nodeType, i, j - 1) &&
                        !INDEX(nodeType, i - 1, j - 1))
                    {
                        // 5
                        Fij[i * N + j] = 1 * ((3. / 4 * yU - 3 - xL) + (3 - xL)) * yU / 2 / (h_1 * h_2);
                    }
                    else if (!INDEX(nodeType, i, j - 1) &&
                             !INDEX(nodeType, i - 1, j - 1))
                    {
                        // 6
                        Fij[i * N + j] = 1 * ((3 - xL) * (-4. / 3 * xL + 4) / 2) /
                                         (h_1 * h_2);
                    }
                    else if (INDEX(nodeType, i, j - 1) &&
                             INDEX(nodeType, i - 1, j - 1))
                    {
                        // 7
                        Fij[i * N + j] = 1 * (((-3. / 4 * yU + 3 - xL) + (-3. / 4 * yD + 3 - xL)) * h_2 / 2) /
                                         (h_1 * h_2);
                    }
                    else if (!INDEX(nodeType, i, j - 1) &&
                             INDEX(nodeType, i - 1, j - 1))
                    {
                        // 8
                        Fij[i * N + j] = 1 * (((-4. / 3 * xL + 4 - yD) * (-3. / 4 * yD + 3 - xL)) / 2) /
                                         (h_1 * h_2);
                    }
                }
                else if (xL < 0 && xR > 0)
                {
                    /*	9 - 36	*/
                    if (yD < 0)
                    {
                        // 9
                        if (elem && INDEX(nodeType, i, j - 1))
                            Fij[i * N + j] = yU / h_2;
                    }
                    else if (yD > 0)
                    {
                        // 10
                        if (elem && INDEX(nodeType, i, j - 1))
                            Fij[i * N + j] = 1;
                    }
                    else if (yD > 0 && yU < 4)
                    {
                        // 11
                        if (!elem && INDEX(nodeType, i - 1, j) &&
                            INDEX(nodeType, i - 1, j - 1) &&
                            !INDEX(nodeType, i, j - 1))
                        {
                            Fij[i * N + j] = 1 * (h_1 * h_2 - ((xR + 3. / 4 * yU - 3) * (yU + 4. / 3 * xR - 4) / 2) - (-xL + 3. / 4 * yU - 3) * (yU - 4. / 3 * xL + 4) / 2) /
                                             (h_1 * h_2);
                        }
                    }
                    else if (yD > 0 && yU > 4) // 12
                    {
                        if (!elem && INDEX(nodeType, i - 1, j) &&
                            INDEX(nodeType, i - 1, j - 1) &&
                            !INDEX(nodeType, i, j - 1))
                        {
                            Fij[i * N + j] = 1 * ((((-4. / 3 * xR + 4 - yD) + (4 - yD)) * xR / 2) + ((4. / 3 * xL + (4 - yD)) * (-xL / 2)) / (h_1 * h_2));
                        }
                    }
                    else if (yU < 4) // 13-27
                    {
                        // 13-16
                        if (yD > 0 && 4 * xL - 3 * yU + 12 < 0 && 4 * xL - 3 * yD + 12 < 0) // 13
                        {
                            Fij[i * N + j] = 1 * (((-3. / 4 * yD + 3 + (-3. / 4 * yU + 3)) * h_2) / (h_1 * h_2));
                        }
                        else if (yD > 0 && 4 * xL - 3 * yU + 12 < 0 && 4 * xL - 3 * yD + 12 > 0) // 14
                        {
                            Fij[i * N + j] = 1 * (((-3. / 4 * yD + 3 + (-3. / 4 * yU + 3)) * h_2 / 2) + ((-xL * h_2) - ((-xL + 3. / 4 * yU - 3) * (yU - 4. / 3 * xL - 4) / 2)) / (h_1 * h_2));
                        }
                        else if (yD < 0 && 4 * xL - 3 * yU + 12 < 0 && xL < -3) // 15
                        {
                            Fij[i * N + j] = 1 * ((3 + (-3. / 4 * yU + 3)) * yU) / (h_1 * h_2);
                        }
                        else if (yD < 0 && 4 * xL - 3 * yU + 12 < 0 && xL > -3) // 16
                        {
                            Fij[i * N + j] = 1 * ((3 + (-3. / 4 * yU + 3) * yU / 2) + (-xL * yU) - ((-xL + 3. / 4 * yU - 3) * (yU - 4. / 3 * xL - 4) / 2) / (h_1 * h_2));
                        }
                        // 17-18
                        else if (4 * xL - 3 * yD + 12 < 0 && xR < 3 && 4 * xR + 3 * yU - 12 > 0)
                        {
                            if (yD > 0) // 17
                            {
                                Fij[i * N + j] = 1 * (((-3. / 4 * yD + 3) + (-3. / 4 * yU + 3)) * h_2 / 2 + (xR * h_2) - (yU + 4. / 3 * xR - 4) * (xR + 3. / 4 * yU - 3) / 2) / (h_1 * h_2);
                            }
                            else if (yD < 0) // 18
                            {
                                Fij[i * N + j] = 1 * ((-xL + 3 - 3. / 4 * yU) * yU - (3 - xR * (-4. / 3 * xR + 4)) / 2) / (h_1 * h_2);
                            }
                        }
                        // 19
                        else if (yD < 0 && xL > -3 && xR < 3 && 4 * xL - 3 * yU + 12 < 0 && 4 * xR + 3 * yU - 12 > 0)
                        {
                            Fij[i * N + j] = (h_1 * yU - (xR + 3. / 4 * yU - 3) * (yU + 4. / 3 * xR - 4) / 2 - (-xL + 3. / 4 * xL - 3) * (yU - 3. / 4 * xL + 4) / 2) / (h_1 * h_2);
                        }
                        else if (4 * xR + 3 * yU - 12 < 0 && 4 * xL - 3 * yU + 12 < 0) // 20-23
                        {
                            if (yD < 0 && xL > -3) // 20
                            {
                                Fij[i * N + j] = (xR * yU + (3 + 3 - 3. / 4 * yU) * yU / 2) / h_1 * h_2;
                            }
                            else if (yD < 0 && xL > -3) // 21
                            {
                                Fij[i * N + j] = (yU * h_1 - (yU - 4. / 3 * xL - 4) * (-xL + 3. / 4 * yU - 3) / 2) / (h_1 * h_2);
                            }
                            else if (yD > 0 && 4 * xL - 3 * yD + 12 < 0) // 22
                            {
                                Fij[i * N + j] = 1 * ((xR * h_2 + ((-3. / 4 * yU + 3) + (-3. / 4 * yD + 3)) * h_2 / 2) / (h_1 * h_2));
                            }
                            else if (yD > 0 && 4 * xL - 3 * yD + 12 > 0) // 23
                            {
                                Fij[i * N + j] = (h_1 * h_2 - (yU - 4. / 3 * xL - 4) * (-xL + 3. / 4 * yU - 3) / 2) / (h_1 * h_2);
                            }
                        }
                        else if (4 * xL - 3 * yU + 12 > 0) // 24-27
                        {
                            if (yD < 0 && xR > 3) // 24
                            {
                                Fij[i * N + j] = ((-xL * yU) + (3 - 3. / 4 * yU + 3) * yU / 2) / (h_1 * h_2);
                            }
                            else if (yD < 0 && xR > 3) // 25
                            {
                                Fij[i * N + j] = (h_1 * yU - (yU + 4. / 3 * xR - 4) * (xR + 3. / 4 * yU - 3) / 2) / (h_1 * h_2);
                            }
                            else if (yD > 0 && 4 * xR + 3 * yD - 12 > 0) // 26
                            {
                                Fij[i * N + j] = (-xL * h_2 + (-3. / 4 * yD + 3 + (-3. / 4 * yU + 3)) * (h_2 / 2)) / (h_2 * h_1);
                            }
                            else if (yD > 0 && 4 * xR + 3 * yD - 12 < 0) // 27
                            {
                                Fij[i * N + j] = (h_1 * h_2 - (yU + 4. / 3 * xR - 4) * (xR + 3. / 4 * yU - 3) / 2) / (h_1 * h_2);
                            }
                        }
                    }
                    else if (yU > 4) // 28-35
                    {
                        if (yD > 0) // 28-31
                        {
                            if (4 * xR + 3 * yD - 12 > 0 && 4 * xL - 3 * yD + 12 < 0) // 28
                            {
                                Fij[i * N + j] = (4 - yD) * (-3. / 4 * yD + 3) / (h_1 * h_2);
                            }
                            else if (4 * xR + 3 * yD - 12 > 0 && 4 * xL - 3 * yD + 12 > 0) // 29
                            {
                                Fij[i * N + j] = (4 - yD) * (-3. / 4 * yD + 3) / 2 + ((4. / 3 * xL + 4 - yD) + (4 - yD)) * (-xL / 2) / (h_1 * h_2);
                            }
                            else if (4 * xR + 3 * yD - 12 < 0 && 4 * xL - 3 * yD + 12 < 0) // 30
                            {
                                Fij[i * N + j] = (4 - yD) * (-3. / 4 * yD + 3) / 2 + (-4. / 3 * xR + 4 - yD) + (4 - yD) * (xR / 2) / (h_1 * h_2);
                            }
                            else if (4 * xR + 3 * yD - 12 < 0 && 4 * xL - 3 * yD + 12 > 0) // 31
                            {
                                Fij[i * N + j] = ((-4. / 3 * xR + 4 - yD) + (4 - yD) * (xR / 2) + ((4. / 3 * xL + 4 - yD) + (4 - yD)) * (-xL / 2)) / (h_1 * h_2);
                            }
                        }
                        else if (yD < 0) // 32-35
                        {
                            if (4 * xR + 3 * yD - 12 > 0 && xL < -3) // 32
                            {
                                Fij[i * N + j] = 12 / (h_1 * h_2);
                            }
                            else if (4 * xR + 3 * yD - 12 > 0 && xL > -3) // 33
                            {
                                Fij[i * N + j] = (12 - (3 + xL) * (4. / 3 * xL + 4) / 2) / (h_1 * h_2);
                            }
                            else if (4 * xL - 3 * yD + 12 < 0 && xR > 3) // 34
                            {
                                Fij[i * N + j] = (12 - (3 - xR) * (-4. / 3 * xR + 4) / 2) / (h_1 * h_2);
                            }
                            else if (xR < 3 && xL > -3) // 35
                            {
                                Fij[i * N + j] = (12 - (3 - xR) * (-4. / 3 * xR + 4) / 2 - (3 + xL) * (4. / 3 * xL + 4) / 2) / (h_1 * h_2);
                            }
                        }
                    }
                }
            }
        }
    }
}
double
coefA(int i, int j, int *nodeType)
{
    double x = A_1 + (j - 0.5) * h_1;
    double yU = A_2 + (i + 0.5) * h_2;
    double yD = A_2 + (i - 0.5) * h_2;
    int U = INDEX(nodeType, i, j - 1);
    int D = INDEX(nodeType, i - 1, j - 1);

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
coefB(int i, int j, int *nodeType)
{
    double y = A_2 + (i - 0.5) * h_2;
    double xR = A_1 + (j + 0.5) * h_1;
    double xL = A_1 + (j - 0.5) * h_1;
    int R = INDEX(nodeType, i - 1, j);
    int L = INDEX(nodeType, i - 1, j - 1);
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

double *operatorA(double *w, int *nodeType)
{
    double *res = (double *)malloc(sizeof(double) * ((M + 1) * (N + 1)));
    {
        for (int i = 1; i < M; i++)
        {
            for (int j = 1; j < N; j++)
            {
                INDEX(res, i, j) = -(coefA(i, j + 1, nodeType) * (INDEX(w, i, j + 1) - INDEX(w, i, j)) / h_1 - coefA(i, j, nodeType) * (INDEX(w, i, j) - INDEX(w, i, j - 1)) / h_1) / (h_1) - (coefB(i + 1, j, nodeType) * (INDEX(w, i + 1, j) - INDEX(w, i, j)) / h_2 - coefB(i, j, nodeType) * (INDEX(w, i, j) - INDEX(w, i - 1, j)) / h_2) / (h_2);
            }
        }
    }
    return res;
}

int main()
{
    init();
    int *nodeType = (int *)malloc(sizeof(int) * M * N);
    double *Fij = (double *)malloc(sizeof(double) * ((M + 1) * (N + 1)));

    nodeTypeDef(nodeType);
    FijDef(Fij, nodeType);
    int iterNum = 0;
    double *wNewer = (double *)malloc(sizeof(double) * ((M + 1) * (N + 1)));
    wNewer = w;
    double norm;
    double *residual = funVSubtract(operatorA(w, nodeType), Fij);
    double *Ar = operatorA(residual, nodeType);
    double tau = funVScalProd(Ar, residual, h_1, h_2) / pow(funVNorm(Ar, h_1, h_2), 2);
    // Generalized minimal residual method
    double DVMH_start = 0.0, DVMH_end = 0.0;
    // DVMH_start = dvmh_wtime();
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

    // dvmh_barrier();
    // DVMH_end = dvmh_wtime();
    // printf("DVMH TIME:%f\n", DVMH_end - DVMH_start);
    printf("Execution succeed!\n");
    printf("Total iteration: %d\n", iterNum);
    printf("Time elapsed: %f ms\n", elapsed_time * 1000);

    // for (int i = M; i >= 0; i--)
    // {
    //     for (int j = 0; j < N + 1; j++)
    //     {
    //         printf("%f ", INDEX(wNewer, i, j));
    //         if (j <= N - 1)
    //         {
    //             printf("%f, ", INDEX(wNewer, M - i, j));
    //         }
    //         else
    //         {
    //             printf("%f\n", INDEX(wNewer, M - i, j));
    //         }
    //     }
    // }
    return 0;
}