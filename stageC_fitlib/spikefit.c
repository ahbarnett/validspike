/* C library for spike fitting, a re-implementation of Matlab fitting.
   See spikefit.mw for any missing documentation, and MEX calling setup.
   Barnett 2/12/15 start.
*/

#include "spikefit.h"

//----------------------------------------------------------------------------
void spikemod(double* W, int M, int T, int K, int fac, int Ns, int* l,
	      double* t, double* a, int Nt, double* F, int subF, int* iran)
/* Algorithm is synthesis/spikemodel with 1-indexed changed to 0-indexed.

 Inputs:   
   W - (M*T*K) waveforms 3d array stored contiguously,
               indexed by M fast, T med, K slowest.
   M,T,K - number of channels, waveform time points (upsampled), spike types
   fac - (integer) upsampling ratio
   Ns - number of spikes
   l,t,a - pointers to Ns-length arrays of model params (label, time, ampl)
   Nt - number of time points in signal
   subF - flag: 0 zeros F first and adds spikes; 1 subtracts spikes from
          whatever was in F already, 2 assumes Ns=1 overwriting only supp(F)

 Outputs:
   F  - (M*Nt) signal array M fast, Nt slow.
   iran - length-2 integer array of [ilo,ihi] for the most recent spike added

 Barnett 2/13/15. Better scaling for big clips (limiting iput) 3/19/15.
 subF=2 case where only writes to one-spike supp(F), 4/30/15
 */
{
  int i,s,iput,m,ilo,ihi, icenW = (T-1)/2;  // center, 0-indexed
  double al[MAXSPIKESPEREVENT];      // local copy of ampls a
  //printf("%d %d %d %d %d %d\n",M,T,K,fac,Ns,Nt); // debug
  int pad = 1;      // padding for good luck (integer rounding)
  
  if (subF==1)
    for (s=0; s<Ns; ++s) al[s] = -a[s];       // negate ampls
  else {
    for (s=0; s<Ns; ++s) al[s] = a[s];        // usual ampls
    if (subF==0)
      for (i=0; i<M*Nt; ++i) F[i] = 0.0;      // initialize all F=0
    else if (subF==2 && Ns==1) {              // init F=0 only in 1-spike support
      ilo = MAX(0,floor(t[0] - icenW/(double)fac - pad));
      ihi = MIN(Nt,ceil(t[0] + icenW/(double)fac + pad));
      for (i=ilo; i<ihi; ++i)
	for (m=0; m<M; ++m)
	  F[i*M+m] = 0.0;     // initialize F=0 only in one-spike support of F
    }
  }
    
  for (s=0; s<Ns; ++s) {               // loop over spikes
    // range of output indices to sweep...
    ilo = MAX(0,floor(t[s] - icenW/(double)fac - pad));
    ihi = MIN(Nt,ceil(t[s] + icenW/(double)fac + pad));
    for (iput=ilo; iput<ihi; ++iput) {
      int iget = round(icenW + fac*(iput - t[s]));   // t starts at zero
      if (iget>=0 && iget<T) {         // only read values lying in waveform
	int oput = iput*M;             // write offset for the channel loop
	int oget = (iget + (l[s]-1)*T)*M; // read offset (ie for this time & k)
	for (m=0; m<M; ++m)            // loop over channels, copy col vec
	  F[oput+m] = F[oput+m] + al[s] * W[oget+m]; // use local ampl
      }
    }
  }
  iran[0] = ilo; iran[1] = ihi;  // output the signal index range
}

//----------------------------------------------------------------------------
double nll(double* W, int M, int T, int K, int fac, int Ns, int* l,
	   double* t, double* a, int Nt, double* Y, double* F, double eta,
	   int locflag, double* srt)
/* NLL - return negative log likelihood for spike model w/ iid Gaussian noise
   Only one type of simple model for now.

 Inputs:
   W   - (M*T*K) waveforms 3d array stored contiguously
   M,T,K - number of channels, waveform time points (upsampled), spike types
   fac - (integer) upsampling ratio
   Ns  - number of spikes to read from model params
   l,t,a - pointers to Ns-length arrays of model params
   Nt  - number of time points in signal
   Y   - (M*Nt) signal array stored continuously
   F   - (M*Nt) needed workspace, returns signal for the model params
                (unless locflag=1 then F is garbage outside supp(F))
   eta - noise sqrt(variance), only aspect of noise model used for now
   locflag - (integer) if 0, compute NLL usual way; if 1, use local update.
                       If 2, compute NLL usual way & initialize srt
   srt - size Nt double array of sum squared Y per time pt. Is updated on output.

 Output: minus log (likelihood of signal Y given spike model F)

 Barnett 2/15/15
 */
{
  int i, iran[2], m;
  double x, J = 0.0;          // J=NLL
  if (locflag==0) { // standard NLL eval
    spikemod(W, M, T, K, fac, Ns, l, t, a, Nt, F, 0, iran); // subF init's F
    for (i=0; i<M*Nt; ++i) {    // l_2^2 norm over channels & times...
      x = Y[i]-F[i];
      J += x*x;
    }
  } else if (locflag==2) {   // same as locflag=0 but also initializes srt
    spikemod(W, M, T, K, fac, Ns, l, t, a, Nt, F, 0, iran); // subF init's F
    for (i=0; i<Nt; ++i) {    // l_2^2 norm. over times...
      double srti = 0.0;
      for (m=0; m<M; ++m) { // channels
	x = Y[i*M+m]-F[i*M+m];
	srti += x*x;
      }
      srt[i] = srti;
      J += srti;
    }
  } else if (locflag==1) {
    spikemod(W, M, T, K, fac, Ns, l, t, a, Nt, F, 2, iran); // subF=2 not init F
    for (i=0; i<iran[0]; ++i)  // sum srt in region before supp(F)
      J += srt[i];
    for (i=iran[1]; i<Nt; ++i) // .. and after
      J += srt[i];
    // only write to support of F...
    //printf("iran = %d %d\n",iran[0],iran[1]); // debug
    for (i=iran[0]; i<iran[1]; ++i) {    // l_2^2 norm over supp(F)
      double srti = 0.0;
      for (m=0; m<M; ++m) { // channels
	x = Y[i*M+m]-F[i*M+m];
	srti += x*x;
      }
      J += srti;
    }
  }  
  J /= 2.0*eta*eta;  // overall noise scaling
  return J;
}

//----------------------------------------------------------------------------
void fitonesp(double* W, int M, int T, int K, int fac, int* lb,
		double* tb, double* ab, int Nt, double* Y, double eta,
	      double tpad, double* Fb, double* Jbest, double* nlps,
	      int locflag, double* srt)
/* FITONESP - fit one spike, based on stageC_fitting/fitonespike.m

 Inputs:
   W   - (M*T*K) waveforms 3d array stored contiguously
   M,T,K - number of channels, waveform time points (upsampled), spike types
   fac - (integer) upsampling ratio
   Nt  - number of time points in signal
   Y   - (M*Nt) signal array stored continuously
   eta - noise sqrt(variance), only aspect of noise model used for now
   tpad - (double) padding in time points to stop spike bumping into the end
   nlps - (double size K) -log prior probabilities for each spike type 1..K
   locflag - (integer) if 0, compute NLL usual way; if 1, use local update.
   srt - size Nt double array of squared resid per time pt. Is updated.

 Outputs:
   lb (int), tb (double), ab (double) - 1-spike best model params (lb in 1..K)
   Fb (M*Nt double) - model sig output at best-fit params, stored contiguously
   Jbest - best-fit objective func = neg log lik

 Algorithm : exhaustive search of max posterior prob over time t and label l.
   Doesn't store J's for all param choices (as would be needed for
   marginalization), just the running best one.

 todo: check only on original sample grid, then upsample only near ones close
       to min J
 todo: ampl fitting, via (d/da)NLL explicit formula.
 todo: Bayesian evidence, or analytic amplitude marginalization.

 Barnett 2/13/15-2/15/15. Fix l is 1-indexed not 0-indexed bug, 2/18/15.
 log priors 3/12/15. locflag & srt 4/30/15
*/
{
  int l,i,m;
  double t;
  double a = 1.0;          // fix all ampls for now
  double* F = (double *) malloc(sizeof(double)*M*Nt);      // alloc for nll
  int iran[2];        // for index range of supp(F)
  
  *Jbest = INFINITY;
  for (l=1; l<=K; ++l) {          // waveforms (outer loop since spaced in RAM)
    for (t=tpad; t<=(double)Nt-1-tpad; t+=1.0/fac) {    // timeshifts
      double J = nll(W, M, T, K, fac, 1, &l, &t, &a, Nt, Y, F, eta, locflag, srt);
      J = J + nlps[l-1];      // -log posterior = -log lik - log prior
      // printf("l=%d\tt=%.1f\tJ=%.3f\n",l,t,J); // debug
      if (J<*Jbest) {
	*Jbest = J;
	*lb = l; *tb = t; *ab = a; // save the best params
      }
    }
  }
  // re-eval Fb best model, and maybe update srt in supp(F)...
  spikemod(W, M, T, K, fac, 1, lb, tb, ab, Nt, Fb, 0, iran);
  if (locflag) {
    for (i=iran[0]; i<iran[1]; ++i) {    // l_2^2 norms over supp(F)
      double srti = 0.0;
      for (m=0; m<M; ++m) { // channels
	double x = Y[i*M+m]-Fb[i*M+m];
	srti += x*x;
      }
      srt[i] = srti;
    }
  }
  free(F);
}



//----------------------------------------------------------------------------
void fitgreedy(double* W, int M, int T, int K, int fac, int* Ns, int* lb,
		double* tb, double* ab, int Nt, double* Y, double eta,
	       double tpad, int maxNs, double* J, double* R, double* nlps,
	       int locflag)
/* FITGREEDY - greedy fitting of multiple spikes in one signal window

 Inputs:
   W   - (M*T*K double) waveforms 3d array stored contiguously
   M,T,K - number of channels, waveform time points (upsampled), spike types
   fac - (int) upsampling ratio
   Nt  - (int) number of time points in signal
   Y   - (M*Nt double) signal array stored continuously
   eta - (double) noise sqrt(variance), only aspect of noise model used for now
   tpad - (double) padding in time points to stop spike bumping into the end
   maxNs - (int) max # spikes to add into model
   nlps - (double size K) -log prior probabilities for each spike type 1..K
   locflag - (integer) if 0, compute NLL usual way; if 1, use local update.

 Outputs:
   Ns (int) - number of spikes in best fit
   lb (int), tb (double), ab (double) - multi-spike best model param arrays
     (each must be allocated up to size maxspikes). Are padded with 0 or NANs.
   J (double, maxNs+1) - best-fit obj func after each spike = neg log lik
     Note: J[0...Ns] contains greedy fitting history.
           If Ns<maxNs, J[Ns+1] contains the best J for Ns+1th spike not added.
   R (M*Nt double) - residual at best-fit params, stored contiguously

 Algorithm: multiple calls to fitonesp.

 See also: ../fitgreedyspikes.m which is similar but outputs signal not resid.

 Barnett 2/16/15. fixed bug that if *Ns=0, was wrong 2/17/15. log priors 3/12/15
*/
{
  int i, s;
  double *srt;
  double* Fb1 = (double *) malloc(sizeof(double)*M*Nt);      // alloc for sig
  if (locflag) {
    srt = (double *) malloc(sizeof(double)*Nt);      // nll per time-pt
    // here lb, tb, ab, Fb1 are dummies since it's a 0-spike model...
    J[0] = nll(W, M, T, K, fac, 0, lb, tb, ab, Nt, Y, Fb1, eta, 2, srt);
    // here the locflag = 2 makes nll initialize srt
  } else
    J[0] = nll(W, M, T, K, fac, 0, lb, tb, ab, Nt, Y, Fb1, eta, locflag, srt);

  for (i=0; i<M*Nt; ++i) R[i] = Y[i];   // init residual R (is the fitted thing)
  *Ns = 0;                              // init output # spikes
  for (s=0; s<maxNs; ++s) {     // loop over adding spikes...
    // note writes best params into correct place in lb,tb,ab, and J arrays...
    fitonesp(W, M, T, K, fac, lb+s, tb+s, ab+s, Nt, R, eta, tpad, Fb1, J+s+1,
	     nlps, locflag, srt);
    if (J[s+1]<J[s]) {        // accept the new spike
      *Ns = s+1;              // update number of spikes
      for (i=0; i<M*Nt; ++i) R[i] -= Fb1[i];    // subtract model from resid
    } else {              // clean up and stop
      lb[s] = 0; tb[s] = NAN; ab[s] = NAN; // leave the not-accepted J
      break;
    }
  }
  free(Fb1);
}


//----------------------------------------------------------------------------
void multifitgreedy(double* W, int M, int T, int K, int fac, int* Ns,
		    int* lb, double* tb, double* ab, int *Tc, int Nc,
		    double* Y, double eta,
		    double tpad, int maxNs, int wantR, double* J, double* R,
		    double* nlps, int locflag)
/* MULTIFITGREEDY - greedy fitting many spikes in many clip windows, OpenMP

 Inputs:
   W   - (M*T*K) waveforms 3d array stored contiguously
   M,T,K - number of channels, waveform time points (upsampled), spike types
   fac - (integer) upsampling ratio
   Tc (int size Nc) - number of time points in each signal clip
   Nc - number of signal clips
   Y   - (M*sum(Tc)) signal array stored contiguously
   eta - noise sqrt(variance), only aspect of noise model used for now
   tpad - (double) padding in time points to stop spike bumping into the end
   maxNs - max # spikes to add into model for fitting each clip
   wantR - 1 to write to residual array (size of Y); 0 not to.
   nlps - (double size K) -log prior probabilities for each spike type 1..K
   locflag - (integer) if 0, compute NLL usual way; if 1, use local update.
                       if 2 (default), switch based on which is fastest.

 Outputs:
   Ns (int size Nc) - number of spikes in best fit model for each clip
   lb (maxNs*Nc int), tb (maxNs*Nc double), ab (maxNs*Nc double) -
     multi-spike best model param arrays (only first Ns[c] valid in col c).
     maxNs is the fast (row) variable.
   J (double, (maxNs+1)*Nc) - history of greedy best-fit objective func
      ( = neg log lik) for each clip. maxNs+1 is the fast (row) axis.
   R (double, M*sum(Tc)) optional residual array for all clips contiguously

 Algorithm: multiple calls to fitonesp, distributed across cores.
 Notes:
 if want fast, think re flops & RAM mvmt. Here I'm doing a bit of local cacheing
 todo: test differing Tc's

 todo: make this choose locflag = (Nt > 100) on a per-clip basis if locflag==2.

 See also: test_multifitgreedy.m, ../fitgreedyspikes.m (on which based)

 Barnett 2/17/15. neg log priors 3/12/15
*/
{
  int Nt, ns, i, c, o, *l;  // c will loop over clips
  double *t,*a,*Jloc,*Yloc,*Rloc;
  
 // set up p array (since omp not sequential) = indexes along time axis (not M)
  int* p = (int*)malloc(Nc*sizeof(int));
  p[0] = 0; for (c=0;c<Nc-1;++c) p[c+1] = p[c] + Tc[c];

  int maxNt = 0;
  for (c=0; c<Nc; ++c) if (Tc[c]>maxNt) maxNt=Tc[c]; // compute max over Tc

  // note shared is default. The idea of splitting the parallel and for comes
  // from: http://stackoverflow.com/questions/2352895/how-to-ensure-a-dynamically-allocated-array-is-private-in-openmp
  //omp_set_num_threads(8);
#pragma omp parallel private(c,Yloc,Rloc,Jloc,l,t,a,i,ns,Nt,o)
  {
    // alloc stuff once for each thread... (saves re-mallocing inside loop)
    Yloc = (double*) malloc(M*maxNt*sizeof(double));
    Rloc = (double*) malloc(M*maxNt*sizeof(double));
    Jloc = (double*) malloc((maxNs+1)*sizeof(double));
    l = (int*) malloc(maxNs*sizeof(int));
    t = (double*) malloc(maxNs*sizeof(double));
    a = (double*) malloc(maxNs*sizeof(double));

#pragma omp for schedule(dynamic) // dynamic important: effort varies as O(Nt^2)
    for (c=0; c<Nc; ++c) {
      Nt = Tc[c];            // Nt for this clip
      for (i=0; i<maxNs; ++i) {   // reset w/ NAN the local arrays
	l[i] = 0; t[i] = NAN; a[i] = NAN; // there's no int NAN; use 0
      }
      for (i=0; i<maxNs+1; ++i) Jloc[i] = NAN;
      for (i=0; i<M*Nt; ++i) Yloc[i] = Y[i+p[c]*M]; // get from distant Y array
      // we just made local Y copy since fitgreedy accesses Yloc a lot...
      int locflagc = locflag;
      if (locflag==2) locflagc = (Nt>100); // choose on a per-clip basis
      fitgreedy(W, M, T, K, fac, &ns, l,t,a, Nt, Yloc, eta, tpad, maxNs,
		Jloc, Rloc, nlps, locflagc);
      // now uncaching: copy local arrays into the big distant arrays...
      Ns[c] = ns;
      //printf("    ns=%d Jbest=%.3f\n",ns,Jloc[ns]);
      //for (i=0;i<maxNs+1;++i) printf("Jloc[%d] = %.3f\n",i,Jloc[i]);
      o = c*maxNs;         // offset in the maxNs*something arrays
      for (i=0; i<maxNs; ++i) {
	lb[i+o] = l[i]; tb[i+o] = t[i]; ab[i+o] = a[i];
      }
      o = c*(maxNs+1);         // J has 1 more row!
      for (i=0; i<maxNs+1; ++i) J[i+o] = Jloc[i]; // copy the J history into col
      if (wantR)       // copy local resid into big output array
	for (i=0; i<M*Nt; ++i) R[i+p[c]*M] = Rloc[i];
    }
    free(Yloc); // free for each thread
    free(Rloc);
    free(Jloc);
    free(l);
    free(a);
    free(t);
  } // end omp parallel
  free(p);
}

