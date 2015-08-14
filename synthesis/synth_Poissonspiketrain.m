function [Y p] = synth_Poissonspiketrain(wf,N,rates,noi,tpad,seed,o)
% SYNTH_POISSONSPIKETRAIN - make possibly noisy Poisson-firing spike time-series
%
% [Y p] = synth_Poissonspiketrain(wf,N,rates,noi,tpad,seed,opts) generates time
%  series using the waveform set in wf, with given firing rates.
%  Note: the rates determine an exact number of spikes for each type, with random
%  times, rather than a true Poisson process with variable number of each type.
%
% Inputs:
%  wf - waveform (classifier) object; defines K = size(wf.W,3);
%  N - number of time-points (overall time = N*wf.d.dt)
%  rates - list of K firing rates (spikes/sec); if <0 then whole list is negated
%          and interpreted as enforcing strict number of spikes of each type
%  noi - (optional) noise model struct (no noise added if empty)
%  tpad - (optional, real) duration (in time samples) to avoid at each end
%  seed - (optional, int) random number seed (unchanged if empty
%         randomizes if 'shuffle')
%  opts - (optional) struct controls options such as:
%         opts.ampl - controls relative iid Gaussian amplitude stddev (default 0)
%         opts.truePois - if true (default false), use true Poisson process with
%                     intensities given by rates (# spikes is variable)
%         opts.rateunits - 's' firing rates are in per-sample, 'Hz' in per second
%                     (default is 'Hz', which demands field wf.d.dt exists)
%
% Outputs:
%  Y - M*N multichannel time-series
%  p - true spikes parameter object: p.t = times (0-indexed times in
%        sample units), p.l = labels in 1...K, p.a = amplitudes
%
% See also: LOADDEMODATA, SETUP_NOISEMODEL

% todo: make opt for a true Poisson process w/ variable # of spikes
% todo: make opt to respect a refractory violation
% Barnett 6/2/15. 6/8/15: truePois & centered on (N-1)/2 not N/2

if nargin<1, test_synth_Poissonspiketrain; return; end
if nargin<4, noi=[]; end
if nargin<5 || isempty(tpad), tpad=0; end
if nargin>=6 && ~isempty(seed), rng(seed); end
if nargin<7, o = []; end
if ~isfield(o,'ampl'); o.ampl = 0; end
if ~isfield(o,'truePois'); o.truePois = 0; end
if ~isfield(o,'rateunits'); o.rateunits = 'Hz'; end

tmin = tpad; tmax = N-1-tpad; tlen = tmax-tmin; % real time interval containing t
if min(rates)<0, Ne = -round(rates);  % Ne = # of each type (specified)
  o.truePois = 0;                     % override Poisson
else
  if strcmp(o.rateunits,'Hz')
    rates = rates*wf.d.dt;            % convert to firings per timesample
  end
  if ~o.truePois                      % fixed Ne deduced from rates
    Ne = round(tlen*rates);
  end
end
K = size(wf.W,3);
p.l = []; p.t = []; p.a = [];  % param struct: labels, times (unsorted), ampls
for k=1:K                      % append each spike type in turn  
  if o.truePois                % make each Ne a Poisson random variable
    tk = []; t = tmin - log(rand(1))/rates(k);  % generate first firing
    while (t<tmax)
      tk = [tk t]; t = t - log(rand(1))/rates(k); % subsequent firings
    end
    Ne(k) = numel(tk);
  else                      % fixed # unif rand times (not strict Poisson rate)
    tk = tmin + rand(1,Ne(k))*tlen;
  end
  p.t = [p.t tk];
  p.l = [p.l k*ones(1,Ne(k))];
  p.a = [p.a ones(1,Ne(k))+o.ampl*randn(1,Ne(k))];
end
[p.t,i] = sort(p.t); p.l=p.l(i); p.a=p.a(i); % time-order the spikes

Y = spikemod(wf, p, N);  % run fwd model

if ~isempty(noi), Y = Y + noisesample(noi); end  % possibly add noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_synth_Poissonspiketrain
wf = loaddefaultwaveforms;
N = round(wf.d.samplefreq);         % 1 second of time
noi = setup_noisemodel(wf.d,N,25);  % choose noise std dev eta
K = size(wf.W,3);
'fake Poisson (fixed # events per type):'
[Y p] = synth_Poissonspiketrain(wf,N,100*ones(1,K),noi,2,0);
fprintf('pops: (should be fixed) '), histc(p.l,1:K)
if exist('spikespy'), spikespy({Y,p.t,p.l,'synth Poisson firing'}); end;

o.truePois = 1; tpad = 100; 'true Poisson:'
[Y p] = synth_Poissonspiketrain(wf,N,100*ones(1,K),noi,tpad,[],o);
fprintf('pops: (should vary) '), histc(p.l,1:K)
fprintf('t range %.1f to %.1f (should be in %d to %d)\n',min(p.t), max(p.t),tpad, N-tpad)
