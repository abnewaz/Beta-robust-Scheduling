#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>



#define  n 15
#define m 2
#define row ((long long int)(pow(m, n + 1) - 1))
#define column (2*m+1)


float T;

float estimate_T(float mu, float var)
{
    return mu + (float)sqrt(var * m);
}

int main()
{
    int i, j, k, r, Job[n], counter = 1, rcnt = 0;
    float Mu[n], Var[n];
    float Mu_sum = 0, Var_sum = 0;

    srand(time(NULL));

    for (i = 0; i < n; i++) //  Initialization
    {
        Mu[i] = rand() % 15 + 10;
        Var[i] = rand() % 4 + 1;
        Job[i] = i;
        Mu_sum += Mu[i];
        Var_sum += Var[i];
    }

    T = estimate_T((float)(Mu_sum / m), (float)(Var_sum / m));
    j = 0;
    float S[row][column];


    for (j = 0; j < (long long int)(pow(m, n + 1) - 1); j++)
    {
        for (k = 0; k < 2 * m + 1; k++)
        {
            S[j][k] = 0.0;
        }

    }


    int index = 1;
    int c;
    for (j = 0; j < n; j++)
    {

        for (i = 0; i < m; i++)
        {

            for (k = (long long int)(pow(m, j) - 1); k < (long long int)(pow(m, j + 1) - 1); k++)
            {

                //printf("%d ", 2 * m + 1);
                for (c = 0; c < 2 * m + 1; c++)     //copying
                {

                    S[index][c] = S[k][c];
                }
                //printf("here\n");

                S[index][0] = j + 1;
                //for (c = 1; c < 2 * m + 1; c = c + 2)
                //{
                S[index][m*i + 1] += Mu[j];
                S[index][m*i + 2] += Var[j];
                //}
                index++;
            }
            printf("\n");
        }
    }
    for (j = 0; j < (long long int)(pow(m, n + 1) - 1); j++)
    {
        for (k = 0; k < 2 * m + 1; k++)
        {
            printf("%f ", S[j][k]);
        }
        printf("\n");
    }



    return 0;
}
