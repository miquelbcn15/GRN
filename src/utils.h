#ifndef _UTILS_
#define _UTILS_

/* defining the parameters */

#define S 500.
#define E 5.
#define r 0.4
#define kk 2.
#define KK 10.
#define a 5.

#define N 2   /* number chemical species */
#define M 4   /* number of reactions */

typedef struct gen {
    double R;
    double k1;
    double k2;
    double b11;
    double v11;
} GenReg;

double uniform();
void computeRatesGillespie(double w[], double x[], GenReg genReg);
double generateAlpha(double w[]);
int selectReaction(double w[], double z, double alpha);
void updateConc(double x[], double u[M][N], int rx);
double binomialQSSA(double p, double e);
void computeRatesQSSA(double w[], double x[], GenReg genReg);
double probQSSA(double ns, GenReg genReg);

#endif
