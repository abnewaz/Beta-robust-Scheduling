#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

int n =3;
int m = 2;
#define row  (int)(pow(m, n))
#define column  (int)(2 * m + 1)
double T;

double estimate_T(double mu, double var)
{
    return mu + (double)sqrt(var) * m;
}

double KL_divergence(double m1, double s1, double m2, double s2)
{
    double res, up, mid, log_s2_s1, s1_sqr, m1_m2, m1_m2_sqr, s2_sqr;
    log_s2_s1 = log(s2 / s1);
    s1_sqr = pow(s1, 2);
    m1_m2 = m1 - m2;
    m1_m2_sqr = pow(m1_m2, 2);
    s2_sqr = 2 * (pow(s2, 2));
    up = s1_sqr + m1_m2_sqr;
    mid = up / (s2_sqr);
    res = log_s2_s1 + mid - .5;

    return res;
}

double CDF(double z)
{

    return 0.5 * erfc(-z * (double)sqrt(0.5));
}

double likelihood(double M_load_mu[], double M_load_var[])
{
    int i;
    double ans = 1;
    for (i = 0; i < m; i++)
    {
        ans *= (double) CDF((T - M_load_mu[i]) / (double)sqrt(M_load_var[i]));
    }

    return ans;
}



int main()
{

    int i, j, k, Job[n];
    double Mu_sum = 0, Var_sum = 0;

    int temp2, min_load_index;
    double KL_distance[n], Mu[n], Var[n], C_max, temp, M_load[m], assignment[m], min_load;
    double M_load_mu[m], M_load_var[m];
    double sum_KL_dist = 0;
    srand(time(NULL));

    for (i = 0; i < n; i++) //  Initialization
    {
        Mu[i] = rand() % 180 + 90;
        Var[i] = rand() % 15 + 1;
        Job[i] = i;
        Mu_sum += Mu[i];
        Var_sum += Var[i];
        KL_distance[i] = KL_divergence(Mu[i], (double)sqrt(Var[i]), 0, 1);
        sum_KL_dist += KL_distance[i];

        printf("For job %d, N(%.2f,%.2f)\n", i,Mu[i],Var[i]);
    }

    T = estimate_T((double)(Mu_sum / m), (double)(Var_sum / m));

    // ******** ALG start here ********

    for (i = 0; i < m; i++)
    {
        M_load[i] = 0;
    }
    for (i = 0; i < n - 1; i++) // sorting based on KL_score
    {
        for (j = i + 1; j < n; j++)
        {

            if (KL_distance[i] < KL_distance[j])
            {
                temp = KL_distance[i];
                KL_distance[i] = KL_distance[j];
                KL_distance[j] = temp;

                temp2 = Job[i];
                Job[i] = Job[j];
                Job[j] = temp2;

                temp = Mu[i];
                Mu[i] = Mu[j];
                Mu[j] = temp;

                temp = Var[i];
                Var[i] = Var[j];
                Var[j] = temp;


            }
        }

    }



    for (j = 0; j < n; j++) // job assigning in minimum loaded machine
    {
        min_load = sum_KL_dist;

        for (i = 0; i < m; i++) // finding the minimum loaded machine
        {

            if (M_load[i] < min_load)
            {
                min_load = M_load[i];
                min_load_index = i;
            }

        }


        M_load[min_load_index] += KL_distance[j];
        M_load_mu[min_load_index] += Mu[j];
        M_load_var[min_load_index] += Var[j];



    }

    double result = likelihood(M_load_mu, M_load_var);

    printf("---------------------------------\n");


    // ******** ALG stop here ********
    maxl = 0;
    double S[1][5] = {0.0};
    // put job1 in mc 1
    S[1][0] = 1;
    S[1][1] += Mu[0];
    S[1][2] += Var[0];
    S[1][3] = 0;
    S[1][4] = 0;

    //put job 2 in mc 1  // mc 2

    S[1][0] = 2;
    S[1][1] += Mu[1];
    S[1][2] += Var[1];
    S[1][3] = 0;
    S[1][4] = 0;

    //put job 3 in mc 1 // mc2 branch

    S[1][0] = 3;
    S[1][1] += Mu[2];
    S[1][2] += Var[2];
    S[1][3] = 0;
    S[1][4] = 0;

    // calculate maxl
    //update maxl if necessary
    //free



    return 0;
}



