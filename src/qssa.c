#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "utils.h"

/* 
 * To compile:
 *  gcc -o qssa -g -Wall qssa.c utils.c -lgsl -lgslcblas -lm
 */

int main (int argc, char* argv[]) {
    const gsl_rng_type *T;
    gsl_rng *s;
    
    /* create generator chosen by the enviornment variable GSL_RNG_TYPE */
    gsl_rng_env_setup();

    T = gsl_rng_default;
    s = gsl_rng_alloc(T);

    double kk1  = a / (kk * sqrt(KK));
    double kk2  = 1.;
    double vv11 = 1.;
    double RR   = r / (kk * sqrt(KK));
    double b11  = 1. / S;

    GenReg genReg;
    genReg.R   = RR;
    genReg.k1  = kk1;
    genReg.k2  = kk2;
    genReg.b11 = b11;
    genReg.v11 = vv11;

    double tmax = 0.;
     
    if (argc<2
        || sscanf(argv[1], "%lf", &tmax)!=1
        ) {
        fprintf(stderr,"%s tmax \n", argv[0]);
        return -1;
    }
    
    /* Re-scaling the t_max */
    tmax = tmax * kk;

    double x[N] = {0., 0.};
    double w[M];
    double u[M][N] = {
        {+1., +0.},
        {-1., +0.},
        {-2., +0.},
        {+2., +0.},
    };
    
    /* Initializing ... */
    double t  = 0.;
    double ns = x[0] / S;
    double p  = probQSSA(ns, genReg); 
    x[1] = gsl_ran_binomial(s, p, E);
    clock_t start_t = clock();

    fprintf(stdout, "# t tau x1 x2 \n");
    while (t < tmax) {   /* QSSA Algorithm */
        /* Gillespie step considering slow variable */
        /* Compute the transitions rates */
        computeRatesQSSA(w, x, genReg);

        /* Generate alpha */
        double alpha = generateAlpha(w);

        /* Generating random values, tau and z */
        double tau = uniform();
        tau = - log(1. - tau) / alpha;
        double z = uniform();

        /* Selecting reaction */
        int rx = selectReaction(w, z, alpha);

        /* Updating chemical concentrations and time for the slow variable  */
        updateConc(x, u, rx);
        t += tau;

        /* Updating the fast variable */
        ns = x[0] / S; 
        p = probQSSA(ns, genReg);
        x[1] = gsl_ran_binomial(s, p, E);

        /* printing results */
        fprintf(stdout, "%lf %lf %lf %lf\n", t/kk, t, x[0], x[1]);
    }
    clock_t end_t = clock();
    fprintf(stderr, "Total execution time QSSA method:\n%lf\n\n", 
            (double)(end_t - start_t)/CLOCKS_PER_SEC);
    return 0;
}
