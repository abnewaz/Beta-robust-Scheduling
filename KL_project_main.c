/*
This code is prepared by
Mahabub Uz Zaman,
Anisa Bente Newaz

for the paper Titled
"Beta-Robust Stochastic Scheduling for identical parallel machines using kullback-leibler divergence and reinforcement learning"

*/
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
int m = 20;
int n = 100;
double T = 0;

int LSPT(double mu[], double var)
{

    return 0;
}

double KL_divergence(double m1, double s1, double m2, double s2)
{
    double res, m, n, log_s2_s1, s1_sqr, m1_m2, m1_m2_sqr, s2_sqr;
    log_s2_s1 = log(s2 / s1);
    s1_sqr = pow(s1, 2);
    m1_m2 = m1 - m2;
    m1_m2_sqr = pow(m1_m2, 2);
    s2_sqr = 2 * (pow(s2, 2));
    m = s1_sqr + m1_m2_sqr;
    n = m / (s2_sqr);
    res = log_s2_s1 + n - .5;
    return res;
}

double CDF(double z)
{
    // return (1 / 2) * (1 + erf(z / sqrt(2)));
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

double estimate_T(double mu, double var)
{
    return mu + (double)sqrt(var * m);
}

int main()
{
    int i, j, Job[n], temp2, min_load_index;
    double KL_distance[n], Mu[n], Var[n], C_max, temp, M_load[m], assignment[m], min_load;
    double M_load_mu[m], M_load_var[m];
    double Mu_sum = 0, Var_sum = 0;

    srand(time(NULL));

    for (i = 0; i < n; i++) //  Initialization
    {
        Mu[i] = rand() % 15 + 10;
        Var[i] = rand() % 4 + 1;
        Job[i] = i;
        Mu_sum += Mu[i];
        Var_sum += Var[i];
        KL_distance[i] = KL_divergence(Mu[i], (double)sqrt(Var[i]), 0, 1);
    }

    T = estimate_T((double)(Mu_sum / m), (double)(Var_sum / m));

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

                temp = Mu[i];
                Mu[i] = Mu[j];
                Mu[j] = temp;

                temp = Var[i];
                Var[i] = Var[j];
                Var[j] = temp;

                temp2 = Job[i];
                Job[i] = Job[j];
                Job[j] = temp2;
            }
        }
    }

    for (j = 0; j < n; j++) // job assigning in minimum loaded machine
    {
        min_load = 10000;
        for (i = 0; i < m; i++) // finding the minimum loaded machine
        {
            if (M_load[i] < min_load)
            {
                min_load = M_load[i];
                min_load_index = i;
            }
        }
        M_load[min_load_index] += KL_distance[j];
        // Job[j] = i;
        M_load_mu[min_load_index] += Mu[j];
        M_load_var[min_load_index] += Var[j];
    }

    double result = likelihood(M_load_mu, M_load_var);

    //printf("checkpoint 1");
    printf("%lf", result);

    return 0;
}
