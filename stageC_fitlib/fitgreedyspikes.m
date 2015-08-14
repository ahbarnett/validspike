function [p J0 Jbest f] = fitgreedyspikes(wf,Y,noi,maxspikes,opts)
% FITGREEDYSPIKES - max likelihood multi-spike greedy fit to an event window
%
% [p J0 Jbest f] = fitgreedyspikes(wf,Y,noi,maxspikes,opts)
%
% Inputs: (first three same as fitonespike)
%  wf - waveforms struct containing: W = M*T*K waveforms
%                                    d = raw EC data struct (d.A not needed)
%                                    fac = upsampling ratio
%                                    freqs = 1*K firing freq estimates
%  Y - M*Nt data array to fit (M channels, Nt time points)
%  noi - noise model struct, minimally noi.eta = sqrt(variance)
%  maxspikes - max # spikes (optional)
%  opts - (optional) controls:
%         opts.verb = 0 (silent, default), 1 (diagnostic output), 2 (&plots)
%
% Outputs:
%  p - best-fit param struct, eg p.t times, p.a ampls, p.m = identities
%  J0 - 0-spike obj func, ie MSE going in.
%  Jbest - best-fit objective func =neg log lik
%  f - model signal output at best-fit params p
%
% Note: made OBSOLETE by spikefitlib/multifitgreedy.m produced from spikefit.mw
%
% Barnett 2/11/15
if nargin<1, test_fitgreedyspikes; return; end
if nargin<4 || isempty(maxspikes), maxspikes = 5; end
if nargin<5, opts = []; end
if ~isfield(opts,'verb'), opts.verb = 0; end

[J0,~,f] = negloglik(wf,[],Y,noi); Jbest=J0; % NLL & signal for 0-spike model
p.l = []; p.a = []; p.t = [];
for s=1:maxspikes
  [p1 Jbest1 f1] = fitonespike(wf,Y,noi,opts);
  if opts.verb
    fprintf('s=%d: p1 = %.3g %.3g %d, \tJbest1=%.3f\n',s,p1.t,p1.a,p1.l,Jbest1)
  end
  if Jbest1<Jbest      % accept the new spike
    f = f + f1; Y = Y - f1;  % "move" spike from data to model signal
    Jbest = Jbest1;    % keep new best NLL
    p.t = [p.t p1.t]; p.a = [p.a p1.a]; p.l = [p.l p1.l]; % append to params
  else, return          % reject and exit with current p, Jbest, f
  end
end

%%%%
function test_fitgreedyspikes   % just single run, no failure stats, for now.
% tweaked from test_fitonespike
load data/waveforms_e_fac4  % chooses upsampling fac
W = W/max(abs(W(:)));  % L-infty normalize
wf.W = W; wf.fac = fac; [M,T,K] = size(W); d.M = M;   % setup wf, # channels
Nt = 30;
pe.t = [13.35 10.5 20]; pe.l = [2 1 3]; pe.a=[1 1 1]; % pick times & identities
Y = spikemodel(wf, pe, Nt);
figure; plot(Y', '-','color', .5*[1 1 1], 'linewidth',2);
noi = setup_noisemodel(d,Nt,0.1);
[p J0 Jbest f] = fitgreedyspikes(wf,Y,noi);
hold on; plot(f', '-'); drawnow;
fprintf('no noise: J0=%.3f, Jbest=%.3f. Params:\n',J0,Jbest); pe, p
fprintf('success = %d\n',spikesetmatch(pe,p))

if 1
Y = Y + noisesample(noi);
[p J0 Jbest f] = fitgreedyspikes(wf,Y,noi);
fprintf('noise: J0=%.3f, Jbest=%.3f. Params:\n',J0,Jbest); p
fprintf('success = %d\n',spikesetmatch(pe,p))
fprintf('expected residual in noi model J=%.3f\n',numel(Y)/2)
figure; plot(Y', '-','color', .5*[1 1 1], 'linewidth',2);
hold on; plot(f', '-'); drawnow;
end
