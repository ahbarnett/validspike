/* C library for "glutton" spike fitting (greedy on whole dataset)
   Barnett 5/1/15
*/

#include "spikefit.h"

// local declarations
void waveformnorms(double* W, int M, int T, int K, int fac, double* nrms);

#define MAXK 1000   // # spike types
#define MAXNTHREADS 1000

//----------------------------------------------------------------------------
void fillscore(double* W, int M, int T, int K, int fac, double* Y, int Nt,
	       double* tsh, int Nsh, double* S, double eta, int* shflags)
/* fill S score matrix given data Y, waveform W, and eta noise level

 Arguments: 
   W - (M*T*K double input) waveforms 3d array stored contiguously,
               indexed by M fast, T med, K slowest.
   M,T,K - (integer inputs) number of channels, waveform time points
          (upsampled), spike types
   fac - (integer input) upsampling ratio beta
   Y - (M*Nt double input) signal data
   Nt - (integer input) length of signal
   tsh - (size-Nsh double input) array of time shift values to compute S at
   Nsh - (integer input) number of time shift values
   S - (K*Nsh double output) score matrix, ordered K fast, Nsh slow
   eta - sqrt of variance in iid noise model
   shflags - (size-Nsh int input). If jth entry is 0, skip that tsh index, else
              if 1, compute jth col of S as usual. If -1 write a Nan

   Cost is O(Nsh.K.T.M) and Nsh = fac.N in applications

   Barnett 5/1/15. JFM's single norm test 5/14/15 (turned out not to help)
*/
{
  int i,j,k,m,iW,ilo,ihi,icenW = (T-1)/2;  // center, 0-indexed
  double pad = 2.0;       // time spillover padding either end in samples
  double* Wnrms; Wnrms = (double*)malloc(sizeof(double)*K);
  
  waveformnorms(W, M, T, K, fac, Wnrms); // get ||W_k|| for k=1..K
  
  for (j=0;j<Nsh;++j) {    // loop over time shifts in S
    if (shflags[j]==1) {      // only compute desired cols of S
      ilo = MAX(0,floor(tsh[j] - icenW/(double)fac - pad)); // index range in Y
      ihi = MIN(Nt,ceil(tsh[j] + icenW/(double)fac + pad));
      double S0t = 0.0;         // compute -(l2 norm)^2 of local signal w/o spike
      for (i=ilo; i<ihi; ++i) {    // loop over needed nearby signal
	iW = round(icenW + fac*(i - tsh[j]));   // t offset; t starts at zero
	if (iW>=0 && iW<T) {         // only read values lying in waveform
	  for (m=0; m<M; ++m)            // loop over channels
	    S0t -= Y[i*M + m]*Y[i*M + m];
	}
      }
      double vnrm = sqrt(-S0t);  // v is vector of local Ymt signal vals
      for (k=0; k<K; ++k) {
	//printf("vnrm=%g, Wnrmk=%g\n",vnrm,Wnrms[k]); // debug
	double tmp = -2*vnrm + Wnrms[k];
	if (tmp>0.0) {  // S must be >0, Jeremy's norm test... hardly ever pass
	  S[k + j*K] = tmp*Wnrms[k];    // write lower bnd for S_kt
	} else {
	  S[k + j*K] = S0t;        // initialize S. note k 0-indexed
	  for (i=ilo; i<ihi; ++i) {
	    iW = round(icenW + fac*(i - tsh[j])); // time index offset
	    if (iW>=0 && iW<T) {         // only read values lying in waveform
	      int oW = (iW + k*T)*M; // W read offset (ie for this time & k)
	      for (m=0; m<M; ++m) {           // loop over channels
		double x = Y[i*M + m] - W[oW+m];
		S[k + j*K] += x*x;
	      }
	    }
	  }
	}
	S[k + j*K] *= 1/(2.0*eta*eta);       // scale for likelihood
      }
    } else if (shflags[j]==-1) {
      for (k=0; k<K; ++k) S[k + j*K] = NAN;  // write Nan col of S
    }                                        // elseif shflags[j]=0 don't write S
  }
  free(Wnrms);
}


//----------------------------------------------------------------------------
void smartscore(double* W, int M, int T, int K, int fac, double* Y, int Nt,
		double* tsh, int Nsh, double* S, double eta, int skip)
/* Faster way to do fillscore using initial coarse grid
   
   doc see fillscore: todo
   skip - (integer input) 1,2,.. gives factor to skip in tsh list
   
   called by fillscore.m created by gluttonfit.mw

   Barnett 5/14/15
*/
{
  int j,k;
  int* shflags; shflags = (int*)malloc(sizeof(int)*Nsh);
  
  for (j=0;j<Nsh;++j) shflags[j] = -1;  // fill S with nans off the coarse grid
  for (j=0;j<Nsh;j+=skip) shflags[j] = 1; // initial coarse grid
  fillscore(W, M, T, K, fac, Y, Nt, tsh, Nsh, S, eta, shflags);

  if (skip>1) { // setup for local fill-ins
    int* shflagsloc; shflagsloc = (int*)malloc(sizeof(int)*(skip-1));
    for (j=0;j<skip-1;++j) shflagsloc[j] = 1;
    
    // now fill in S values around a t at which coarse Skt<0, for any k
    for (j=0;j<Nsh-skip;j+=skip) {
      int fillme = 0;           // if coarse S dips, trigger a fill-in of S
      for (k=0;k<K;++k)
	if (S[j*K +k]<0.0 || S[(j+skip)*K+k]<0.0) fillme = 1;
      if (fillme) {
	fillscore(W, M, T, K, fac, Y, Nt, tsh+j+1, skip-1, S+K*(j+1), eta,
		  shflagsloc);        // only write into local bit of S
      }
    }
  }
  free(shflags);
}
 

//----------------------------------------------------------------------------
void waveformnorms(double* W, int M, int T, int K, int fac, double* nrms)
/* Compute norms of downsampled waveforms 1...K over mt axes.
   The min over 1/fac fractional time-shifts is taken
   
 Arguments: 
   W - (M*T*K double input) waveforms 3d array stored contiguously,
               indexed by M fast, T med, K slowest.
   M,T,K - (integer inputs) number of channels, waveform time points
          (upsampled), spike types
   fac - (integer input) upsampling ratio beta
   nrms - (length-K double output) desired norms ||W_k||_2, note not squared

   Barnett 5/14/15
*/
{
  int m,k,t,off,j;
  for (k=0;k<K;++k) {
    nrms[k] = INFINITY;
    int ind = k*M*T;    // index offset for this type k
    for (off=0;off<fac;++off) {  // loop over fac fractional timeshifts
      double thisnrm = 0.0;
      for (t=off;t<T;t+=fac)
	for (m=0; m<M; ++m) {
	  double w = W[ind+t*M+m];
	  thisnrm += w*w;
	}
      //printf("%g\n",thisnrm);
      if (thisnrm<nrms[k]) nrms[k]=thisnrm;       // keep lowest over offsets
    }
    nrms[k] = sqrt(nrms[k]);  // turn from sum of squares into a norm
  }
}


//----------------------------------------------------------------------------
void locvalidmins(double* S, int K, int Nsh, double* nlps, int* jt,
		  int* l, double* val, int* Nmin)
/* find local valid minima in S score matrix, using local (t+-1) properties only

   jt (int), l (int), val (double) outputs, must be allocated before calling.
        time-shift indices, labels in 1...K, and S-values of valid local
	minima found. output S-vals include the neg log prior (nlps).
	Index jt is 0-indexed, but l runs from 1..K.
   Nmin (int output), # minima found

 */
{
  int j,k;
  int c = 0;     // counter find "valid" local minima
  for (j=1;j<Nsh-1;++j) {
    double minS = INFINITY;
    for (k=0;k<K;++k) if (S[j*K+k]<minS) minS = S[j*K+k];   // S col min 
    for (k=0;k<K;++k) {
      int i = j*K+k;    // index
      if (S[i]<-nlps[k] && S[i]<=minS && S[i]<=S[i-K] && S[i]<=S[i+K]) {
	jt[c] = j; l[c] = k+1;     // convert type to 1-indexed
	val[c] = S[i]+nlps[k];       // note includes neg log prior for type
	++c;
      }
    }
  }
  *Nmin = c;
}


//---------------------------------------------------------------------------
int minisgammaloc(int* jt, int* l, double* val, double* tsh, int N, int j,
		  double gamma)
/*
  Returns true if the jth LVM is global with respect to the window +-gamma either
  side of the center time tsh(jt(j)).

Input arguments:
  jt - list of indices of LVM times in the tsh array
  l - list of spike types of LVMs
  val - list of S-values of LVMs
  N - (int) number of LVMs
  tsh - (double array) time shifts
  j - index of LVM to test (0-indexed)
  gamma - time padding either side (in time samples) for global min test

  see keepgamlocmins.m, gluttonfit.mw/minisgammaloc
*/
{
  int i,keep=1;
  double t = tsh[jt[j]];   // current time shift of LVM
  if (t-tsh[jt[MAX(0,j-1)]] <= gamma | tsh[jt[MIN(N-1,j+1)]]-t <= gamma) {
    for (i=0;i<N;++i)
      if (fabs(tsh[jt[i]]-t)<=gamma & val[i] < val[j]) keep = 0;
  }
  return keep;
}


//---------------------------------------------------------------------------
void gluttonstuffme(double* W, int M, int T, int K, int fac, double* Y, int Nt,
		    double tpad, double eta, int skip, double gamma,
		    double* nlps, double* t, int* l, double* a, int* Ns,
		    int maxNs, int verb)
/* repeated glutton fit rounds on same chunk, single thread version

 Arguments: 
   W - (M*T*K double input) waveforms 3d array stored contiguously,
               indexed by M fast, T med, K slowest.
   M,T,K - (integer inputs) number of channels, waveform time points
          (upsampled), spike types
   fac - (integer input) upsampling ratio beta
   Y - (M*Nt double input/output) signal data; on exit, model residual
   Nt - (integer input) length of signal
   tpad - (double input) padding in sample units at either end of time-shifts
          list
   eta - (double input) sqrt of variance in iid noise model
   skip - (integer input) coarsening factor for skipping through score
          evaluation
   gamma - (double input) time padding either side (in time samples) for
          global score min test
   nlps - (K double input) negative log priors (lambda_l) on firing rates for
          spike types
   t,l,a - output arrays of times (double), labels (int in 1...K), amplitude
          (double). Each must be allocated to at least length maxNs.
   Ns - pointer to number of spikes found, ie length of the t,l,a arrays output.
   maxNs - maximum number of spikes possible.
   verb - (int input) 0,1,... verbosity

   Cost is O(NKTM), memory workspace O(KN)

   Barnett 5/19/15
*/
{
  int j, s, Nsh = round((Nt-1-2*tpad)*fac+2);  // set up time-shifts
  // (note the last tsh overlaps by 2 points the next in multiglutton)
  if (verb>1)
    printf("gluttonstuffme: Nt=%d, Nsh=%d, gamma=%.3g, skip=%d\n",Nt,Nsh,gamma,
	   skip);
  double* tsh; tsh = (double*)malloc(sizeof(double)*Nsh);
  for (j=0;j<Nsh;++j) tsh[j] = tpad + j/(double)fac;  // formula needed later
  double* S; S = (double*)malloc(sizeof(double)*K*Nsh);
  // do initial expensive S fill...
  smartscore(W, M, T, K, fac, Y, Nt, tsh, Nsh, S, eta, skip);

  int Nm, maxNm = Nsh;  // worst-case number of LVMs
  int* jm; jm = (int*)malloc(sizeof(int)*maxNm);   // alloc LVM arrays: tsh ind
  int* lm; lm = (int*)malloc(sizeof(int)*maxNm);   // type
  int* gm; gm = (int*)malloc(sizeof(int)*maxNm);   // whether gamma global min
  double* sm; sm = (double*)malloc(sizeof(double)*maxNm); // s-value
  
  int r = 1, maxr = 20, newspikes = 1;             // todo: opts for maxr
  int c=0; // c is counter for collected spikes
  int cst = 0;            // c where this round started
  while (newspikes && r<=maxr) {     // loop over rounds
    locvalidmins(S, K, Nsh, nlps, jm, lm, sm, &Nm);  // get LVMs
    newspikes = 0;  // count how many new spikes found (valid LVMs)
    for (j=0;j<Nm;++j)  // generate once the array of valid flags
      if (gm[j] = minisgammaloc(jm, lm, sm, tsh, Nm, j, gamma)) // assign gm
	++newspikes;
    if (verb>2)
      printf("round %d: %d valid mins found, %d new spikes\n",r,Nm,newspikes);

    // for gamma-valid LVMs, collect spike and update 
    for (j=0;j<Nm;++j)
      if (gm[j]) {
	  t[c] = tsh[jm[j]];  // add spike to output list (jm is 1-indexed)
	  l[c] = lm[j];
	  a[c] = 1.0;        // dummy ampl (spikemod needs)
	  if (c<maxNs) c++;
	}

    if (c>=maxNs)
      fprintf(stderr,"gluttonstuffme exceeded maxNs=%d number of spikes!\n",maxNs);
    
    // run fwd model with only newest spikes (since cst), subtracting from Y
    int iran[2];  // unused
    spikemod(W,M,T,K,fac,c-cst,l+cst,t+cst,a+cst,Nt,Y,1,iran);

    double jpad = ceil(2.0*T);  // assumes 1/fac in W equals tsh spacing
    if (2*jpad*(c-cst) > Nsh) {  // cheaper to just recompute S afresh
      smartscore(W, M, T, K, fac, Y, Nt, tsh, Nsh, S, eta, skip);
    } else {              
    for (j=0;j<Nm;++j) // overwrite S only in windows around valid minima
      if (gm[j]) {
	int joff = MAX(jm[j]-jpad,0);  // don't fall off bottom of arrays
	int thisNsh = 2*jpad+1;
	if (joff+thisNsh>Nsh) thisNsh = Nsh-joff; // don't fall off top
	smartscore(W, M, T, K, fac, Y, Nt, tsh+joff, thisNsh, S+K*joff,
		   eta, skip);  // note joff offset to the tsh index and in S
      }
    }
    cst = c;   // update counters
    ++r;
  }
  if (verb && r>maxr) printf("gluttonstuffme reached max # rounds, %d\n",maxr);
  *Ns = c;  // total # spikes found
  free(tsh);
  free(S);
  free(jm);
  free(lm);
  free(gm);
  free(sm);
}


//---------------------------------------------------------------------------
void multiglutton(double* W, int M, int T, int K, int fac, double* Y, int N,
		    double tpad, double eta, int skip, double gamma,
		  double* nlps, double* t, int* l, double* a, int* Ns,
		  double* R, int wantR, int verb)
/* repeated glutton fit rounds on entire dataset, multi-thread version

 Arguments: 
   W - (M*T*K double input) waveforms 3d array stored contiguously,
               indexed by M fast, T med, K slowest.
   M,T,K - (integer inputs) number of channels, waveform time points
          (upsampled), spike types
   fac - (integer input) upsampling ratio beta
   Y - (M*Nt double input) signal data
   N - (integer input) length of signal
   tpad - (double input) padding in sample units at either end of time-shifts
          list
   eta - (double input) sqrt of variance in iid noise model
   skip - (integer input) coarsening factor for skipping through score
          evaluation
   gamma - (double input) time padding either side (in time samples) for
          global score min test
   nlps - (K double input) negative log priors (lambda_l) on firing rates for
          spike types 1..K
   t,l,a - output arrays of times (double), labels (int in 1...K), amplitude
          (double). Each must be allocated to length at least N
   Ns - pointer to number of spikes found, ie length of the t,l,a arrays output
   R - (M*Nt double output) model residual. *** note has errors at chunk ends
   wantR - (int input) if 0 don't write to R, otherwise do (not recommended)
   verb - (int input) 0,1,... verbosity

   Cost is O(NKTM/P) for P cores, memory workspace O((K+M)N)

   Barnett 5/22/15
*/
{
  int i,*irng,m,c,o,j;
  //omp_set_num_threads(8);
  int Nc = omp_get_max_threads();             // # chunks, one per thread
  int toff[MAXNTHREADS], Nt[MAXNTHREADS];     // time offsets & size each chunk
  int maxNt = ceil((N-1)/(double)Nc + 2*tpad+1);  // max chunk time length
  toff[0] = 0; Nt[0] = maxNt;                 // set up chunk locations...
  for (c=1;c<Nc;++c) {
    Nt[c] = maxNt; toff[c] = toff[c-1]+maxNt-1-2*tpad;
  }
  Nt[Nc-1] = N-toff[Nc-1];        // last chunk potentially bit shorter
  int maxNs = maxNt-2*tpad;       // at most ~1 spike per time sample!
  int Nss[MAXNTHREADS];           // spikes found per chunk
  if (verb) printf("Nc=%d threads, maxNt=%d, maxNs=%d, tpad=%.3g, eta=%.3g, skip=%d, gamma=%.3g\n",Nc,maxNt,maxNs,tpad,eta,skip,gamma);
  double *Yc;
#pragma omp parallel private(c,Yc,i,irng,o,m)
  {
    Yc = (double*) malloc(M*maxNt*sizeof(double)); // alloc for each thread
    irng = (int*) malloc(2*sizeof(int));

#pragma omp for schedule(dynamic) // dynamic not sure if needed, all sizes same
    for (c=0; c<Nc; ++c) {        // loop over chunks
      o = c*maxNs;                // offset for t,l,a arrays for this chunk

      for (i=0;i<Nt[c];++i)       // copy Y into local Yc: loop over times
	for (m=0;m<M;++m)         // channels
	  Yc[m+M*i] = Y[m+M*(i+toff[c])];
      
      gluttonstuffme(W, M, T, K, fac, Yc, Nt[c], tpad, eta, skip, gamma,
		     nlps, t+o, l+o, a+o, Nss+c, maxNs, verb);  // glutton on Yc
      if (verb>1) printf("done, Nss[%d]=%d\n",c,Nss[c]);

      if (wantR) {                // write middle of chunk Yc to R, todo fix
	irng[0] = (c>0)*tpad;     // time ind range in Yc
	irng[1] = Nt[c]-1 - (c<Nc-1)*tpad; // (ends special, to fill R to ends)
	if (verb>1) printf("irng+toff[%d]=[%d,%d)\n",c,irng[0]+toff[c],irng[1]+toff[c]);
	for (i=irng[0];i<irng[1];++i)  // times
	  for (m=0;m<M;++m)         // channels
	    R[m+M*(i+toff[c])] = Yc[m+M*i];
      }
    }
    free(Yc);  // free for each thread
    free(irng);
  }            // end omp parallel block
  i = 0; // counter to make t,l,a arrays contiguous...
  for (c=0;c<Nc;++c) {
    o = c*maxNs;                  // offset for this chunk
    for (j=0;j<Nss[c];++j) { 
      t[i] = t[o+j] + toff[c]; l[i] = l[o+j]; a[i] = a[o+j];  // time offset
      ++i;
    }
  }
  *Ns = i;  // total spikes
}

