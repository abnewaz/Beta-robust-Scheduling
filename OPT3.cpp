#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cmath>

using namespace std;

int n = 19;
int m = 3;
int row = pow(m, n);
int column =  2 * m + 1 ;
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

double likelihood(vector <double> M_load_mu, vector <double> M_load_var)
{
    int i;
    double ans = 1.0;
    for (i = 0; i < m ; i++)
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
    vector <double> KL_distance (n);
    vector <double> Mu (n);
    vector <double> Var (n);
    double C_max, temp, min_load;
    vector <double> M_load(m), assignment(m), M_load_mu (m), M_load_var (m);
    double sum_KL_dist = 0;
    //srand(time(NULL));

    for (i = 0; i < n; i++) //  Initialization
    {
        Mu[i] = rand() % 180 + 90;
        Var[i] = rand() % 15 + 1;
        Job[i] = i;
        Mu_sum += Mu[i];
        Var_sum += Var[i];
        KL_distance[i] = KL_divergence(Mu[i], (double)sqrt(Var[i]), 0, 1);
        sum_KL_dist += KL_distance[i];

        //cout <<("For job %d, N(%.2f,%.2f)\n", i,Mu[i],Var[i]);
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

    //cout <<("---------------------------------\n");


    // ******** ALG stop here ********



    //int* ptr = malloc((r * c) * sizeof(int));
    //int (*arr)[col] = calloc(row, sizeof *arr)
    //float (*S)[column] = calloc(row, sizeof(*S));
    vector <vector <double> > S(row, vector <double> (column));


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
//            cout <<("%f ", S[j][k]);
//        }
//        cout <<("\n");
//    }

    float OPT_result = -1;
    float OPT_temp_result;
    vector <double> M_load_mu_opt(m),M_load_var_opt(m);
    int best_sol;

    for (j = 0; j < row; j++)
    {
        //cout <<("%d", m);
        for (i = 0; i < m; i++)
        {
            M_load_mu_opt[i] = S[j][i*2 + 1];
            M_load_var_opt[i] = S[j][i*2 + 2];


        }
        OPT_temp_result = likelihood(M_load_mu_opt,M_load_var_opt);

        if (OPT_result < OPT_temp_result)
        {
            OPT_result = OPT_temp_result;
            best_sol = j;

        }


    }


    /* cout <<"OPT result    : %.6f" << OPT_result << endl;
    cout <<"ALG result    : %.6f" << result << endl;
    cout <<"Total Gap     : %.6f\n\n",OPT_result - result);


    cout <<"---------------------------------\n");
    cout <<"Threshold T   : %.2f\n",T);
    cout <<"---------------------------------\n");
    cout <<"Best Solution :\n");
    for (i = 0; i < m; i++)
    {
        cout <<"For machine %d : N(%.2f,%.2f)\n", i, S[best_sol][2*i+1],S[best_sol][i*2+2]);
    }

    cout <<"---------------------------------\n");

    cout <<("ALG Solution :\n");
    for (i = 0; i < m; i++)
    {
        cout <<("For machine %d : N(%.2f,%.2f)\n", i, M_load_mu[i],M_load_var[i]);
    }

    cout <<("---------------------------------\n");
 */

    cout << "successfull" << endl;
    return 0;
}


