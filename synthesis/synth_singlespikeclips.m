function [X L p] = synth_singlespikeclips(wf, Nt, pops, noi, tpad, seed, o)
% SYNTH_SINGLESPIKECLIPS - make noisy random synthetic single-spike event clips
%
% [X L] = synth_singlespikeclips(wf, Nt, pops, noi, tpad, seed, o) returns a 3d
%  array X with clips of equal duration Nt samples, containing synthesized
%  firing events from the waveforms in wf, with numbers of each type given in
%  pops. The ordering of events is randomized, and their identities (labels)
%  returned in L.
%
% [X L p] = ... also returns firing times and amplitudes in p{:}.t and p{:}.a
%
% Inputs:
%  wf    = waveform struct containing M*NT*K array W, integer fac giving
%          upsampling,  etc.
%          note that NT assumed odd.
%  Nt    = number of output time samples per clip.
%          Output grid will be taken from the center of the waveform grid.
%  pops  = 1*K integer, # events to make of each type.
%  noi   = noise model struct (if not present or [], no noise is added).
%          See SETUP_NOISEMODEL.
%  tpad  = buffer in time samples to keep firing time away from clip ends.
%          (default 10)
%          If inf, all firing times are centered on clips (ie no jitter),
%          ie t = (Nt-1)/2 for Nt odd (central sample), or
%             t = Nt/2 - 1 for Nt even (1/2 to left of center)
%  seed  = sets random number generator seed, or randomizes if empty.
%  o (optional struct) sets various options such as:
%      o.ampl - controls relative iid Gaussian amplitude stddev (default 0)
%
% Outputs:
%  X     = M*Nt*Ns data arry where Nt = number of data-rate time samples
%          Ns = sum(pops) total # spikes
%  L     = 1*Ns true labels of events, each an integer in 1,...,K
%  p     = 1*Ns cell array of true parameter structs for each event window,
%          each with fields: l (label), t (time), and a (amplitude).
%          Times are relative to 1st sample output point, in sample units.
%
% Calling without input arguments performs a self-test.
%
% See also: SYNTH_POISSONSPIKETRAIN

% Barnett/Magland 12/11/14. time jitter ta & downsampling 1/5/15-1/6/15
% ahb: noise filtering & autocorr 1/26/15. Use spikemodel, noi, 1/29/15
% 6/9/15: simplified to use Poisson & renamed, retested

if nargin<1, test_synth_singlespikeclips; return; end
if nargin<4, noi=[]; end
if nargin<5 || isempty(tpad), tpad=10; end
if isinf(tpad), tpad = (Nt-1)/2; end % sets no jitter
if nargin>=6 && ~isempty(seed), rng(seed); seed = []; else seed = []; end % only reseed once
if nargin<7, o = []; end
if Nt<1, error('Nt must be positive!'); end
[M NT K] = size(wf.W);
if numel(pops)~=K, error('numel(pops) doesnt match K in wf object'); end
if min(pops)<0, error('at least one entry of pops is negative'); end
Ns = sum(pops);       % # types, total # spikes
L = []; for k=1:K, L = [L, k*ones(1,pops(k))]; end % requested # w/ each label
j = randperm(Ns); L = L(j); % randomize label order
X = nan(M,Nt,Ns);
p = cell(1,Ns);
for j=1:Ns       % for each event window, given params, generate timeseries...
  rates = zeros(1,K); rates(L(j)) = -1;  % tell it exactly one spike of type L_j
  [X(:,:,j) p{j}] = synth_Poissonspiketrain(wf,Nt,rates,noi,tpad,seed,o);
end
%%%%%

function test_synth_singlespikeclips
wf = loaddefaultwaveforms; [M T K] = size(wf.W);
Nt = 30;
eta = 25; tau = .0002; % noise autocorr decay (Prentice '11, around 0.2ms)
noi = setup_noisemodel(wf.d,Nt,eta,tau); %noi = []; % no noise
pops = 1500 - 200*(1:K);   % numbers per type
tpad = 10; %tpad = inf; % no jitter
o = []; o.ampl = 0.1;
[X L p] = synth_singlespikeclips(wf, Nt, pops, noi, tpad, 0, o);
image_spike_shapes(X);
t = nan*L; for j=1:sum(pops), t(j) = p{j}.t; end %cell array struct field -> vec
j=find(L==1); j=j(1:50); plot_spike_shapes(X(:,:,j),[],[],t(j));
