#include <stdio.h>

double kahan_sum(double *xs, size_t n) 
{
    double s = 0.0, c = 0.0;

    for (size_t i = 0; i < n; i++) {
        double y = xs[i] - c;
        double t = s + y;
        c = (t - s) - y;
        s = t;
    }

    return s;
}

int main()
{
    double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};

    double s = kahan_sum(x, sizeof(x)/sizeof(x[0]));
    printf("sum = %.2lf\n", s);

    return 0;
}

