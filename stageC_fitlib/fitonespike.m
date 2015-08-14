function [p Jbest f] = fitonespike(wf,Y,noi,opts)
% FITONESPIKE - maximum likelihood fit for single spike to signal event window
%
% [p Jbest f] = fitonespike(wf,Y,noi,opts) is an obsolete slow matlab code for
%  brute-force exploration of the time-shift and identity of a single spike.
%  Since it uses negloglik, it can handle a general noise model.
%
% Inputs:
%  wf - waveforms struct containing: W = M*T*K waveforms
%                                    d = raw EC data struct (d.A not needed)
%                                    fac = upsampling ratio
%                                    freqs = 1*K firing freq estimates
%  Y - M*Nt data array to fit (M channels, Nt time points)
%  noi - noise model struct, minimally noi.eta = sqrt(variance)
%  opts - optional, controls:  opts.verb : if >1, graph J(t,k) and J0
%
% Outputs:
%  p - best-fit param struct, eg p.t time, p.a ampl, p.m = identity
%  Jbest - best-fit objective func =neg log lik
%  f - model signal output at best-fit params p
%
% Note: has been superceded by fitonesp MEX interface for the iid Gauss noide.

% todo: allow log priors input on each spike type.
%
% Barnett Feb 2015; 2/18/15 plot diagnostics.
if nargin<1, test_fitonespike; return; end
if nargin<4, opts = []; end
if ~isfield(opts,'verb'), opts.verb = 0; end

[M Nt] = size(Y);
K = size(wf.W,3);         % # unit types
dt = 1/wf.fac; tpad = 1.0;   % time spacing, edge padding, to check
ts = 1+tpad:dt:Nt-1-tpad;
J = nan(numel(ts),K);        % store all objective funcs
for i=1:numel(ts) % loop over times
  p.t = ts(i); p.a = 1.0; % fixed ampl for now
  for k=1:K
    p.l = k;   % identity
    J(i,k) = negloglik(wf,p,Y,noi);
  end
end
[Jbest,indbest] = min(J(:));
[ibest,kbest] = ind2sub(size(J),indbest);
p.t = ts(ibest); p.l = kbest;  % copy best into output param struct
if opts.verb>1, figure; plot(ts, J,'+-'); hline(negloglik(wf,[],Y,noi));
  legnum(1:K); xlabel('t'); ylabel('J=NLL'); hold on;
  plot(p.t,Jbest,'k.','markersize',20); title('J(t,k) and J0 in fitonespike')
end
[~,~,f] = negloglik(wf,p,Y,noi);       % get best signal f (wastes small time)
%%%%

function test_fitonespike
%load data/waveforms_e_fac4  % chooses upsampling fac
%load data/waveforms_e_fac1
parallel = 0;                   % tests par toolbox
wf = loaddefaultwaveforms;
wf.W = wf.W/max(abs(wf.W(:)));  % L-infty normalize
[M,T,K] = size(wf.W); d = wf.d;   % setup wf, # channels
Nt = 30;
pe.t = 13.35; pe.l = 2;       % pick time and identity
Y = spikemodel(wf, pe, Nt);
noi = setup_noisemodel(d,Nt,0.1,0.0002);
[p Jbest f] = fitonespike(wf,Y,noi);
fprintf('no noise:\n'); pe, p

Y = Y + noisesample(noi);
[p Jbest f] = fitonespike(wf,Y,noi);
fprintf('iid noise:\n'); p
figure; plot(Y', '-','color', .5*[1 1 1], 'linewidth',2);
hold on; plot(f', '-'); drawnow;

% RAW TIMING
noi = setup_noisemodel(d,Nt,0.1);  % iid Gaussian test (no inv covar matvec)
n=1e2; tic, for i=1:n
  [p Jbest f] = fitonespike(wf,Y,noi);
end,
t=toc; fprintf('noi=iid: %.3g fits/sec (around %.3g Gflmas)\n',n/t,n*T*Nt*d.M*K/t/1e9)
if parallel, n=1e3; tic, parfor i=1:n
  [p Jbest f] = fitonespike(wf,Y,noi);
end,
t=toc; fprintf('par noi=iid: %.3g fits/sec (around %.3g Gflmas)\n',n/t,n*T*Nt*d.M*K/t/1e9)
end

noi = setup_noisemodel(d,Nt,0.1,.0005);  % cov Gaussian test (matvec)
n=1e2; tic, for i=1:n
  [p Jbest f] = fitonespike(wf,Y,noi);
end,
t=toc; fprintf('noi=cov: %.3g fits/sec (around %.3g Gflmas)\n',n/t,n*T*Nt^2*d.M*K/t/1e9)
if parallel, n=1e3; tic, parfor i=1:n
  [p Jbest f] = fitonespike(wf,Y,noi);
end,
t=toc; fprintf('par noi=cov: %.3g fits/sec (around %.3g Gflmas)\n',n/t,n*T*Nt^2*d.M*K/t/1e9)
end

% RELIABILITY & TIMING
%Ns = 1e3; X = nan(M,Nt,Ns); pe = struct([]);
%for i=1:Ns, [ = spikemodel(wf, pe, Nt);
%...
