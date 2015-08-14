function [p Ns Jbest info R] = multifitgreedy(wf,Y,Tc,noi,maxNs,opts)
% MULTIFITGREEDY  greedy fit times, labels in variable-duration clips (MEX/OpenMP)
%
% [p Ns Jbest info R] = multifitgreedy(wf,Y,Tc,noi,maxNs,opts)
%
% INPUTS:
%  wf :    waveform struct with fields: W   - M*T*K array of waveforms
%                                       fac - upsampling factor
%  Y :     M*Ttot signal clips concatenated on time axes, where Ttot = sum(Tc)
%          or, M*Nt*Nc 3d array (equal-Nt case) in which case Tc is ignored.
%  Tc :    lengths in time points of each clip (defines the number of clips)
%  noi :   noise model struct with field noi.eta for std error of iid Gaussian
%  maxNs : (optional, default 5) maximum number spikes to fit per clip
%  opts :  (optional) struct controlling parameters:
%          opts.tpad (default 3.0) clip edge padding in sample time units
%          opts.nlps (default 20 for all types) neg log prior for each spike
%                         type 1..K
%          opts.locflag (default 2): 0 - use std NLL calc, 1 - local update,
%                                    2 - choose fastest on per-clip basis.
%
% OUTPUTS:
%  p :     1*Nc struct array of best-fit params, containing labels l in 1..K,
%          times t (doubles), amplitudes a (doubles).
%  Ns :    1*Nc int number of spikes found for each clip
%  Jbest : 1*Nc best-fit obj func (negative log likelihood) for each clip
%  info :  struct containing diagnostic info:
%          info.Jhist : (maxNs+1)*Nc history of obj func J in s-spike
%          best fits, for s=0 to Ns(c) for this clip, in the c'th column.
%          If Ns(c)<maxNs, the entry Jhist(Ns(c)+1,c) gives not-accepted next
%          spike's best J.
%  R :     M*Ttot complete best-fit residual signal clips (as Y).
%
% See for usage: test_multifitgreedy

Nc = numel(Tc);  % number of time clips
Ttot = sum(Tc);
nrJ = maxNs+1;    % # rows in Jhist
W = wf.W; fac = wf.fac; % enforcing int32 for ints in C breaks MEX interface!
if fac~=floor(fac), error('wf.fac must be positive integer!'); end
[M T K] = size(W);
TK = T*K; W = reshape(W,[M TK]);  % mwrap can't do 3d arrays
if mod(T,2)==0, warning('waveform T should be odd'); end
if numel(size(Y))==3                 % incoming 3d array, so build Tc
  [MY Nt Nc] = size(Y);
  Y = reshape(Y,[MY, Nt*Nc]);        % make it 2d
  Tc = Nt*ones(1,Nc); Ttot = sum(Tc);
else
  [MY Nt] = size(Y);
  if Nt~=Ttot, error('Y must have same # time points as sum of Tc'); end
end
if MY~=M, error('Y must have same number of channels M as W'); end
eta = noi.eta;
if isfield(noi,'covar'), warning('covar noise model not avail in MEX!'); end
if nargin<5, maxNs=5; end
if maxNs>100, maxNs=100; warning('maxNs too big, setting to 100 (a clip is too long!)'); end
if nargin<6, opts = []; end
if isfield(opts,'tpad'), tpad = opts.tpad;
else, tpad = 3.0; end         % default edge padding
if isfield(opts,'nlps'), nlps = opts.nlps;
else nlps = 20*ones(1,K); end           % default neg log priors per spike
if numel(nlps)==1, nlps = nlps*ones(1,K); end

if isfield(opts,'locflag'), locflag = opts.locflag; else locflag = 2; end % def
if locflag<0 | locflag>2, error('locflag must be 0, 1, or 2\n'); end

% watch out: "invalid argument" problem was wantR being logical not double type!

wantR = nargout>4;
if wantR
  mex_id_ = 'multifitgreedy(i double[], i int, i int, i int, i int, o int[xx], o int[xx], o double[xx], o double[xx], i int[], i int, i double[], i double, i double, i int, i int, o double[xx], o double[xx], i double[], i int)';
[Ns, l, t, a, Jhist, R] = sf(mex_id_, W, M, T, K, fac, Tc, Nc, Y, eta, tpad, maxNs, 1, nlps, locflag, 1, Nc, maxNs, Nc, maxNs, Nc, maxNs, Nc, nrJ, Nc, M, Ttot);
else
  R = 0;   % is dummy, might be slightly faster (or save RAM) not to write out.
  mex_id_ = 'multifitgreedy(i double[], i int, i int, i int, i int, o int[xx], o int[xx], o double[xx], o double[xx], i int[], i int, i double[], i double, i double, i int, i int, o double[xx], io double[], i double[], i int)';
[Ns, l, t, a, Jhist, R] = sf(mex_id_, W, M, T, K, fac, Tc, Nc, Y, eta, tpad, maxNs, 0, R, nlps, locflag, 1, Nc, maxNs, Nc, maxNs, Nc, maxNs, Nc, nrJ, Nc);
end

Ns = Ns(:)';

p = [];          % build struct array of params - speed around 1.5e5 clips/sec:
for c=1:Nc, j=1:Ns(c); p(c).l = l(j,c)'; p(c).t = t(j,c)'; p(c).a = a(j,c)'; end

Jbest = Jhist((maxNs+1)*(0:Nc-1) + Ns+1);   % get Ns+1'th entry from each col
info.Jhist = Jhist;
