// headers for spikefit lib. Barnett 5/19/15

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

#define MAXSPIKESPEREVENT 200000         // static, and can use for big clips
// note making this bigger statically allocates more memory

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define ABS(a)     (((a) < 0) ? -(a) : (a))

void spikemod(double* W, int M, int T, int K, int fac, int Ns, int* l,
	      double* t, double* a, int Nt, double* F, int subF, int* iran);
double nll(double* W, int M, int T, int K, int fac, int Ns, int* l,
	   double* t, double* a, int Nt, double* Y, double* F, double eta,
	   int locflag, double* srt);
void fitonesp(double* W, int M, int T, int K, int fac, int* lb,
		double* tb, double* ab, int Nt, double* Y, double eta,
	      double tpad, double* Fb, double* Jbest, double* nlps,
	      int locflag, double* srt);
void fitgreedy(double* W, int M, int T, int K, int fac, int* Ns, int* lb,
		double* tb, double* ab, int Nt, double* Y, double eta,
	       double tpad, int maxNs, double* J, double* R, double* nlps,
	       int locflag);
void multifitgreedy(double* W, int M, int T, int K, int fac, int* Ns,
		    int* lb, double* tb, double* ab, int *Tc, int Nc,
		    double* Y, double eta,
		    double tpad, int maxNs, int wantR, double* J, double* R,
		    double* nlps, int locflag);
