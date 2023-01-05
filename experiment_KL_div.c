#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

double H_dist(double m1, double s1, double m2, double s2)
{
    double res, power;
    double a, b,c,d ;
    a = pow((m1-m2),2);
    b = pow(s1,2)+pow(s2,2);
    c = 2*s1*s2/b;

    power = -.25 * a;
    power = power / b;

    res = 1 - sqrt(c)* exp(power);

    return res;
}


double B_dist(double m1, double s1, double m2, double s2)
{
    double right, left;
    double a, b,c ;
    a = pow((m1-m2),2);
    b = pow(s1,2)+pow(s2,2);
    c = b/(2*s1*s2);

    left = .25 * a;
    left = left / b;
    right = .5 * log(c);

    return right+left;
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

int main()
{
    double instance[25][2];
    double dist;
    int i;
    srand(time(NULL));
    for (i = 0; i < 25; i++){
        instance[i][0] = rand() % 250 + 50;
        instance[i][1] = rand() % 15 + 10;
     dist = H_dist(instance[i][0],instance[i][1],275,12.5);
        printf("For mu = %.2f, var = %.2f, H_dist = %.2f\n",
               instance[i][0],instance[i][1],dist );
    }



    return 0;
}
