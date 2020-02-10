#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "utils.h"

int main(int argc, char *argv[]) {
    /* defining constants */
    double kk1  = a/(kk * sqrt(KK));
    double kk2  = 1.;
    double vv11 = 1.;
    double RR   = r / (kk * sqrt(KK));
    double b11  = kk / (E * S);

    double k1   = kk1 * S * S * b11;
    double k2   = kk2 * b11 * E * S;
    double R    = RR * S * S * E * b11;
    /* double bb11 = 1 / S; */
    double v11  = vv11 * S * S * b11;

    GenReg genReg;
    genReg.R   = R;
    genReg.k1  = k1;
    genReg.k2  = k2;
    genReg.b11 = b11;
    genReg.v11 = v11;

    double tmax = 0.;
     
    if (argc<2
        || sscanf(argv[1], "%lf", &tmax)!=1
        ) {
        fprintf(stderr,"%s tmax \n", argv[0]);
        return -1;
    }

    /* Initializing ... */
    double x[N] = {0., 0.};
    double w[M];
    double u[M][N] = {
        {+1., +0.},
        {-1., +0.},
        {-2., +1.},
        {+2., -1.},
    };

    double t = 0.;
    fprintf(stdout, "# t x1 x2 \n");
    clock_t start_t = clock();
    while (t < tmax) {   /* Gillespie Algorithm */
        /* Compute the transitions rates */
        computeRatesGillespie(w, x, genReg);

        /* Generate alpha */
        double alpha = generateAlpha(w);

        /* Generating random values, tau and z */
        double tau = uniform();
        tau = - log(1 - tau) / alpha;
        double  z  = uniform();

        /* Selecting reaction */
        int rx = selectReaction(w, z, alpha);

        /* Updating chemical concentrations and time */
        updateConc(x, u, rx);
        t += tau;

        /* printing results */
        fprintf(stdout, "%lf %lf %lf\n", t, x[0], x[1]);
    }
    clock_t end_t = clock();
    fprintf(stderr, "# Total Gillespie execution time :\n%lf\n", (double)(end_t -
                start_t)/CLOCKS_PER_SEC);
    return 0;
}
