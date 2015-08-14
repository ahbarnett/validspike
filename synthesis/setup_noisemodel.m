function noi = setup_noisemodel(d,Nt,eta,tau)
% SETUP_NOISEMODEL - initialize a noise model struct
%
% noi = setup_noisemodel(d,Nt,eta) sets up noi as a struct for a iid Gaussian
%  noise model with variance eta^2, using raw data info struct d, and signal
%  window length Nt
%
% noi = setup_noisemodel(d,Nt,eta,tau) sets up noi as a struct for a
%  (within-channel) time-auto-correlated Gaussian noise nodel with variance
%  eta^2, and decay time tau (in secs), using raw data info struct d, and signal
%  window length Nt.
%
% The noi struct is currectly Gaussian only, with properties:
%   noi.eta : std deviation at each signal index
%   noi.sqrtcovar : C^{1/2} where C is covariance matrix between time samples
%                   for each channel (no correlation between channels)
%   noi.invcovar : C^{-1}
%     If sqrtcover or invcovar not present, iid Gaussian asssumed.
%   noi.M, noi.Nt : # channels and # time-points
%
% Without input arguments, setup_noisemodel runs test, plotting matrices.
%
% See also: NOISESAMPLE for sampling from this noise model
%           NEGLOGLIK for computing likelihoods using this noise model

% Barnett 1/29/15, test matrix plots 2/10/15

if nargin<1, test_setup_noisemodel; return; end

noi.M = d.M;
noi.Nt = Nt;
noi.eta = eta;
if nargin==4 && noi.eta>0  % time-correlated within channel
  tau = tau*d.samplefreq;    % convert auto-corr decay time to index units
  covar = noi.eta^2 * toeplitz(exp(-(0:Nt-1)/tau));    % exp(-|t|/tau) autocorr
  noi.sqrtcovar = sqrtm(covar);              % for sampling; O(Nt^3) effort
  noi.invcovar = inv(covar);                 % for likelihood eval
  noi.sqrtinvcovar = inv(noi.sqrtcovar);     % for pre-whitening
  noi.covar = covar;                         % for kicks
end
%%%

function test_setup_noisemodel % just showing covar etc matrices
d.samplefreq = 2e4;
d.M = 1;
noi = setup_noisemodel(d,30,0.1,0.001);
figure;
subplot(2,2,1); imagesc(noi.covar); colorbar; axis equal tight; title('C');
subplot(2,2,2); imagesc(noi.sqrtinvcovar); colorbar; axis equal tight; title('C^{-1/2}');
subplot(2,2,3); imagesc(noi.invcovar); colorbar; axis equal tight; title('C^{-1}');
subplot(2,2,4); imagesc(noi.sqrtcovar); colorbar; axis equal tight; title('C^{1/2}');
% Note C^{-1} is v close to 2nd-deriv matrix, and C^{-1/2} also similar.
% So, prewhitening here means mult by |omega| in freq space.
