#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

int n =6;
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
        Mu[i] = rand() % 250 + 250;
        Var[i] = rand() % 25 + 50;
        Job[i] = i;
        Mu_sum += Mu[i];
        Var_sum += Var[i];
        KL_distance[i] = KL_divergence(Mu[i], (double)sqrt(Var[i]), 0, 1);
        sum_KL_dist += KL_distance[i];
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

    int midp = (int)(n/2);

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

        if (j<=midp){
        M_load[min_load_index] += KL_distance[n-1-j];
        M_load_mu[min_load_index] += Mu[n-1-j];
        M_load_var[min_load_index] += Var[n-1-j];
        }
        else{
        M_load[min_load_index] += KL_distance[j-midp];
        M_load_mu[min_load_index] += Mu[j-midp];
        M_load_var[min_load_index] += Var[j-midp];

        }



    }

    double result = likelihood(M_load_mu, M_load_var);
//    printf("---------------------------------\n");
//    printf("Threshold T   : %.6f\n",T);
//    printf("---------------------------------\n");
//    for (i = 0; i < m; i++)
//    {
//        printf("For machine %d : N(%.2f,%.2f)\n", i, M_load_mu[i],M_load_var[i]);
//    }
//
//    printf("ALG result    : %f\n", result);
    printf("---------------------------------\n");


    // ******** ALG stop here ********



    //int* ptr = malloc((r * c) * sizeof(int));
    //int (*arr)[col] = calloc(row, sizeof *arr)
    float (*S)[column] = calloc(row, sizeof(*S));


    for (j = 0; j < row; j++)
    {
        for (k = 0; k < column; k++)
        {
            S[j][k] = 0;
        }

    }


    int index = 1;
    int start, stop;
    int c;
    for (j = 0; j < n; j++)
    {
        stop = index-1;

        for (i = 0; i < m; i++)
        {
            for (k = 0 ; k<(int)pow(m,j) ; k++)
            {

                for (c = 0; c < column; c++)     //copying
                {
                    S[(index%row)][c] = S[(stop - k)%row][c];
                }

                S[(index%row)] [0] = j + 1;
                S[(index%row)] [2*i+1] += Mu[j];
                S[(index%row)][2*i + 2] += Var[j];

                index++;
            }

        }
    }

//    for (j = 0; j < row; j++)
//    {
//        for (k = 0; k < column; k++)
//        {
//            printf("%f ", S[j][k]);
//        }
//        printf("\n");
//    }

    float OPT_result = -1;
    float OPT_temp_result;
    double M_load_mu_opt[m],M_load_var_opt[m];
    int best_sol;

    for (j = 0; j < row; j++)
    {
        //printf("%d", m);
        for (i = 0; i < m; i++)
        {
            M_load_mu_opt[i] = S[j][i*2 + 1];
            M_load_var_opt[i] = S[j][i*2 + 2];
            //printf("%f\n",M_load_mu[i]);

        }
        OPT_temp_result = likelihood(M_load_mu_opt,M_load_var_opt);
        //printf("   %f\n",ALG_temp_result);
        if (OPT_result < OPT_temp_result)
        {
            OPT_result = OPT_temp_result;
            best_sol = j;
//            printf("---------------------------------\n");
//            printf("Best solution %.2f found in row %d\n",OPT_result,j);



        }
//        else{printf("    Solution %.2f found in row %d\n", OPT_temp_result,j);}

    }


    printf("OPT result    : %.6f\n",OPT_result);
    printf("ALG result    : %.6f\n", result);
    printf("Total Gap     : %.6f\n\n",OPT_result - result);
//    printf("---------------------------------\n");
//    printf("Total Mu sum  : %.6f\n",Mu_sum);
//    printf("Total Var sum : %.6f\n",Var_sum);

    printf("---------------------------------\n");
    printf("Threshold T   : %.6f\n",T);
    printf("---------------------------------\n");
    printf("Best Solution :\n");
    for (i = 0; i < m; i++)
    {
        printf("For machine %d : N(%.2f,%.2f)\n", i, S[best_sol][2*i+1],S[best_sol][i*2+2]);
    }

    printf("---------------------------------\n");

    free(S);
    printf("ALG Solution :\n");
    for (i = 0; i < m; i++)
    {
        printf("For machine %d : N(%.2f,%.2f)\n", i, M_load_mu[i],M_load_var[i]);
    }

    printf("---------------------------------\n");

    return 0;
}


