function [C taus] = crosscorr(l,t,a,o)
% CROSSCORR - estimate cross-correlation vs time shift given times and labels
%
% [C taus] = crosscorr(l,t)
% [C taus] = crosscorr(l,t,a)
% [C taus] = crosscorr(l,t,a,opts)
%  computes K-by-K cross-correlation matrix dependent on time separation tau.
%  The mean is not subtracted, so C_{ij}(tau) has non-negative entries.
%
% Inputs:
%  l - 1D array of labels
%  t - 1D array of firing times
%  a - 1D array of firing amplitudes (set to 1 if empty or absent)
%  opts controls various options:
%   opts.dtau, taumax = tau time bin width and maximum tau, in sample units.
%   dtau should be odd
% Outputs:
%  C - K*K*Nt matrix of cross-correlations (K = max(l), Nt = # timeshift bins)
%  tau - list of timeshift bins
%
% Uses quantization onto sample time grid of spacing 1, in order to get O(N)
% scaling, where N is number of spikes.
% (a pre-sort and binary search could get O(N log N) with more bookkeeping)
%
% todo: rewrite meth=d in C as w/ mda i/o
%
% todo: another meth using K^2 FFTs and inner prods?
% or K FFTs and (ntau)*K^2 inner prods w/ exps (ie slow FT)?
%
% Barnett 3/1/15. docs & default ampls, 4/7/15. new meth=d 4/8/15

if nargin<1, show_crosscorr; return; end
if nargin<3, a = []; end
if isempty(a), a = 0*t+1; end % default ampls
if nargin<4, o = []; end
if ~isfield(o,'dtau'), o.dtau = 20; end
if ~isfield(o,'taumax'), o.taumax = 500; end   % ie 50 bins
taus = 0:o.dtau:o.taumax; taus = [-taus(end:-1:2) taus]; % tau bin centers
taue = [taus-o.dtau/2 taus(end)+o.dtau/2];  % tau bin edges
tend=taus(end)+o.dtau/2; % farthest end of tau bins
ntau = numel(taus); n = (ntau-1)/2;  % # rel time bins, max integer bin #

K = max(l);  % # types
i = l>0; l=l(i); t=t(i); a=a(i); % kill nonpositive labels
[t,i] = sort(t-min(t)+1); l = l(i); a = a(i);  % shift and sort t values
C = zeros(K,K,ntau);
shs = round((1.5-o.dtau)/2:(o.dtau-0.5)/2); % sample shifts for central bin
T = round(max(t)); N = numel(t);  % now make signal array (allow overlaps!)...

meth='d';

if meth=='a'  % version via making discrete A signal array
  fprintf('creating quantized array, size %.3gGB...\n',K*T*8/1e9)
  LA = zeros(K,T); for j=1:N, LA(l(j),round(t(j)))=LA(l(j),round(t(j)))+a(j); end
  %LA = bsxfun(@minus,LA,mean(LA,2));  % make zero-mean - decided against
  fprintf('computing C... ')
  for j=1:N, if mod(j,round(N/10))==0, fprintf('%d%% ',round(100*j/N)); end
    tj = t(j); lj = l(j);
    for d=-n:n
      i = round(tj)+d*o.dtau+shs;
      if min(i)>0 & max(i)<=T  % haven't fallen off end of LA
        a = sum(LA(:,i),2);   % col vec
        if d==0, a(lj) = a(lj) - 1; end    % remove self spike
        C(:,lj,d+n+1) = C(:,lj,d+n+1) + a;
      end
    end
  end
  
elseif meth=='d'  % work directly w/ t,l: faster, 1 min 50bins N=1e6 K=10, 1core
  fprintf('computing C... ')
  j0=1; while t(j0)<tend, j0=j0+1; end % get range 
  j1=N; while t(j1)>T-tend, j1=j1-1; end
  jm = 1;        % start index of events to include rel to jth event
  jp = j0; while t(jp)<2*tend, jp=jp+1; end % end index of same
  for j=j0:j1, if mod(j,round(N/10))==0, fprintf('%d%% ',round(100*j/N)); end
    tj = t(j); lj = l(j);
    while t(jm)<tj-tend, jm=jm+1; end  % update [jm,jp] index range to include...
    while t(jp)<tj+tend, jp=jp+1; end
    i = [jm:j-1, j+1:jp];  % index list in window, omit the current spike j
    tir = t(i)-tj; li = l(i); % times and labels in time window around tj    
    for k=lj:K  % only do lower-tri part since C tau-rev anti-symm
      lik = li==k; if sum(lik)>0
        c = histc(tir(lik),taue); c=c(:);  % last entry from histc is zero
        C(k,lj,:) = squeeze(C(k,lj,:)) + a(j)*c(1:end-1); % annoying sizes
      end
    end
    %ii=(li>=lj);
    %if numel(i)>0  % lower-tri part of C, more vectorized than above
    %  L = zeros(numel(i),K-lj+1); L(1:numel(i),
    %  c = histc(L,taue); c = c(1:end-1,:)'; % exploit histc works on all cols
    %  C(lj:k,lj,:) = C(lj:k,lj,:) + a(j)*c;
    %end
  end
end
fprintf('\n')
