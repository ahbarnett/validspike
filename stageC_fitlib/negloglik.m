function [J DJ f] = negloglik(wf,x,Y,noi)
% NEGLOGLIK  negative log likelihood (obj func) for spike model in clip.
%
% J = negloglik(wf,p,Y,noi) returns J = -log prob(y | x), ie negative log
%  likelihood of signal window data Y given the model parameters x, under
%  noise model noi. This is where the noise model is implemented.
%
% [J DJ f fc] = negloglik(wf,x,y,noi) also returns diagnostic outputs
%
% Inputs:
%  wf  = waveform object (contains stuff spikemodel needs)
%  x   = parameter object (same as p, see spikemodel)
%  Y   = M*Nt signal window
%  noi = noise model object, containing:
%      noi.eta = noise level (for each y component) for iid Gaussian noise model
%
% Outputs:
%  J    = negative log likelihood
%  [DJ   = derivative of J wrt some parameters, not implemented]
%  f    = the model signal vector f(x)
%
% Calling without arguments does self-test
%
% Notes: 1) This is slow, calling the MATLAB spikemodel, and general noise model
% 2) Gaussian iid is the only noise model for now

% Barnett 1/29/15 based on fitwindow01 3/13/14
if nargin<1, test_negloglik; return; end

[M Nt] = size(Y);
f = spikemodel(wf,x,Nt);  % f is M*Nt
if isfield(noi,'invcovar')             % apply the noise model
  d = Y - f;
  J = 0; % loop over channels, have indep contrib...
  for m=1:M, J = J + (1/2)*d(m,:)*(noi.invcovar * d(m,:)'); % quadratic form
  end
else
  J = sum(abs(Y(:)-f(:)).^2)/(2*noi.eta^2);  % -upstairs in exp for Gaussian
end
DJ = [];  % dummy: todo: get deriv of spikemodel wrt all params?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_negloglik  % test noise generation and likelihood
wf = loaddefaultwaveforms; d = wf.d;
wf.W = wf.W/max(abs(wf.W(:)));  % L-infty normalize
[M,T,K] = size(wf.W);
Nt = 30;
p.t = 14.5; p.l = 1;       % pick time and identity
Wused = wf.W(:,:,p.l); sc = max(abs(Wused(:))); % plot y-scale
F = spikemodel(wf, p, Nt);
Y = F; plot_spike_shapes(Y, 'no noise',sc);
noi = setup_noisemodel(d,Nt,0.1); % data object d needed for d.samplefreq, d.M
[J DJ f] = negloglik(wf,p,Y,noi);
fprintf('no noise:\t J = %.3g, should be zero\n',J)

Y = F + noisesample(noi);
plot_spike_shapes(Y, 'iid noise',sc);
J = negloglik(wf,p,Y,noi);
fprintf('iid noise:\t J = %.3g, should be around %.3g\n',J, numel(Y)/2)

tau = .001; noi = setup_noisemodel(d,Nt,0.1,tau);    % correlated noise
Y = F + noisesample(noi);
plot_spike_shapes(Y, 't-corr noise',sc);
J = negloglik(wf,p,Y,noi);
fprintf('t-corr noise:\t J = %.3g, should be around %.3g\n',J, numel(Y)/2)
