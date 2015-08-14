function [t l p wf R] = spikesort_timeseries(Y,samplefreq,copts,smeth,sopts)
% SPIKESORT_TIMESERIES  spike sort a (filtered) multi-channel time-series.
%
% This is the main interface to our few-channel spike sorting algorithms.
%
% [t l p wf R] = spikesort_timeseries(Y,samplefreq)
% [t l p wf R] = spikesort_timeseries(Y,samplefreq,copts,smeth,sopts)
%
% Inputs:
%  Y - M*N real-valued data array (M = # channels, N = # time samples)
%  samplefreq - time sampling rate in samples/sec
%  copts - struct setting options for classifier method;
%          see WAVEFORMS_FROM_TIMESERIES
%  smeth - string setting sorting method; see FIT_TIMESERIES
%  copts - struct setting options for sorting method; see FIT_TIMESERIES
%
% Outputs:
%  t,l - lists of firing times (0-indexed, in sample units) and labels in 1...K
%  p - struct array of all spike parameters
%       p(:).t times, p(:).l labels, p(:).a amplitudes, etc
%  wf - waveform struct containing waveforms used for fitting
%  R - residual signal (same size as Y) after subtraction of spikes
%      (see comments about R in FIT_TIMESERIES)
%
% For example use: DRIVER_TIMESERIES
%
% See also: WAVEFORMS_FROM_TIMESERIES, FIT_TIMESERIES

% Barnett 6/11/15
if nargin==0, driver_timeseries; return; end
if nargin<3 || isempty(copts), copts = []; end
if nargin<4 || isempty(smeth), smeth = 'Gl'; end     % default
if nargin<5 || isempty(sopts), sopts = []; end
if ~isfield(sopts,'verb'), sopts.verb = 1; end

wf = waveforms_from_timeseries(Y,samplefreq,copts);   % build classifier

if nargout>3  % only compute R if asked for
  [t l R] = fit_timeseries(Y,wf,smeth,sopts);     % do fitting
  if sopts.verb, fprintf('\trms resid = %.3g\n',sqrt(mean(R(:).^2))), end
else
  [t l] = fit_timeseries(Y,wf,smeth,sopts);
end
if nargout>2, p = []; p.t = t; p.l = l; end   % package the only params we have
