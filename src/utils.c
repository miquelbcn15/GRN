#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include "utils.h"

double uniform() {
    int randomData = open("/dev/urandom", O_RDONLY);
    if (randomData < 0) {
        fprintf(stderr, "\nERROR: Stopping in open file...\n\n");
        exit(-1);
    }
    unsigned char randomNum;
    size_t result = read(randomData, &randomNum, 1);
    if (result < 0) {
        fprintf(stderr, "\nERROR: Stopping in obtaining random...\n\n");
        exit(-1);
    }
    close(randomData);
    return (double)randomNum/(255.0 + 1.);
}

void computeRatesGillespie(double w[], double x[], GenReg genReg) {
    w[0] = genReg.R + genReg.k1 * x[1]; 
    w[1] = genReg.k2 * x[0];
    w[2] = genReg.b11 * x[0] * (x[0] - 1) * (E - x[1]);
    w[3] = genReg.v11 * x[1];
}

void computeRatesQSSA(double w[], double x[], GenReg genReg) {
    w[0] = genReg.R + genReg.k1 * x[1] / E;
    w[1] = genReg.k2 * x[0] / S;
    w[2] = x[0] / (S * S) * (1. - x[1] / E) * (x[0] - 1.);
    w[3] = genReg.v11 * x[1] / E;
}

double generateAlpha(double w[]) {
    double sum = 0;
    for (int i = 0; i < M; i++) {
        sum += w[i];
    }
    return sum;
}

int selectReaction(double w[], double z, double alpha) {
    double sum = 0;
    double val = alpha * z;
    int i;
    for ( i = 0; i < M; i++) {
        if ( (sum += w[i]) > val ) return i;
    }
    return i;
}

void updateConc(double x[], double u[M][N], int rx){
    for (int i = 0; i < N; i++) {
        x[i] += u[rx][i];
    }
}

double binomialQSSA(double p, double e) {
   double bin = 0;
   double z;
   for (int i = 0; i < e; i++) {
       if ((z = uniform()) < p ) bin++;
   }
   return bin;
}

double probQSSA(double ns, GenReg genReg) {
    return (ns * (ns - genReg.b11)) / (ns * (ns - genReg.b11) + genReg.v11);
}
