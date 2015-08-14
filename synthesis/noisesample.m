function Y = noisesample(noi)
% NOISESAMPLE - draw a sample from the additive noise model
%
% Y = noisesample(noi) where noi is noise model struct returns Y a sample from
%  noise model.
%
% See also SETUP_NOISEMODEL, NEGLOGLIK

% Barnett 1/29/15. Note if noi empty, cannot return a zero vector since don't
% know what size it is.
if nargin<1, negloglik; return; end

if ~isfield(noi,'sqrtcovar')  % iid
  Y = noi.eta * randn(noi.M, noi.Nt);
else                        % t-corr within channels
  Y = randn(noi.M, noi.Nt) * noi.sqrtcovar;   % xforms each row indep 
end
