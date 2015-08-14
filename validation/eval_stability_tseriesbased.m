function [fhat,fsam,info] = eval_stability_tseriesbased(alg, Y, o)
% EVAL_STABILITY_TSERIESBASED   stability metrics for a time-series spike sorter
%
% [fhat,fsam,info] = eval_stability_tseriesbased(alg, opts)
%
% Inputs:
%  alg - function handle of time-series spike sorter to be tested, with
%        interface:
%          [T L] = alg(Y)   where Y is input timeseries (Nchannels x Nsamples)
%          returning T (list of firing times in zero-indexed sample units) and
%          L (list of corresponding labels in 1...K).
%
%  opts.meth : chooses stability metric method:
%
%   'rev' - noise reversal (JFM method, default)
%
%   'add' - spike addition (AHB method). Produces info.Qs as induced confusion
%           matrices.
%           Method options:
%           opts.num_runs : # indep perturbed runs (Nr, default 20)
%           opts.ratescale : fraction of original spiking rates to add spikes at
%                           (default 0.2)
%           opts.ampl : std deviation of relative spike amplitudes about 1
%                      (passed to synth_Poissonspiketrain).
%  General options:
%  opts.Nt : # samples window to use for waveform estimation and add/subtraction
%  opts.upsampfac : upsampling factor (default 3)
%  opts.verb : verbosity (0,1,2,3..)  >=3 is for paper.
%  opts.T, opts.L : override the unpert run of alg using times and labels T,L,
%      for consistent figure generation (experts only)
%
% Outputs:
%    fhat - (1-by-K) best estimates of stability metrics for each neuron type
%    fsam - (Nr-by-K) complete set of samples of metrics for each neuron type
%    info - struct with:
%           Qs - confusion matrix (or cell array of induced confusion matrices)
%           pops, permL2 - populations and best permutation from first full run
%           ... and more
%
% See also: DRIVER_TIMESERIES_STABILITY for example usage, showing use of inline
%           function for alg

% Barnett 7/15/15, add 7/22/15


if nargin<3, o=[]; end
if ~isfield(o,'meth'), o.meth = 'rev'; end
if ~isfield(o,'verb'), o.verb = 1; end
if ~isfield(o,'Nt'), o.Nt = 60; end   % W pull window width in timeseries samples
if isfield(o,'num_runs'), Nr = o.num_runs; else Nr=20; end  % > dozen (Mackay96)
if ~isfield(o,'ratescale'), o.ratescale = 0.2; end    % relative density to add
[M N] = size(Y);
fsam = []; info.Qs = {};

if strcmp(o.meth,'rev')       % ================== noise reversal
  if ~isfield(o,'T')
    fprintf('ss noise-rev unperturbed run\n')
    [T{1} L{1}] = alg(Y);
  else
    fprintf('noise-rev overriding unpert ss run with input data\n')
    T{1} = o.T; L{1} = o.L;   % skip the unpert run (advanced)
  end
  K = max(L{1});
  wf = pull_waveforms_from_tseries(Y,T{1},L{1},o.Nt,o);
  F = spikemod(wf, struct('t',T{1},'l',L{1}), N);      % run ampl-1 fwd model
  Ynr = 2*F - Y;                                       % the noise reversal
  fprintf('ss noise-rev reversed run\n')
  [T{2} L{2}] = alg(Ynr);
  if o.verb>2, info.Ynr = Ynr; info.Tnr = T{2}; info.Lnr = L{2}; end  % for figs
  [info.permL2 info.Qs info.acc info.times] = times_labels_accuracy(T{1},L{1},T{2},L{2},o);  % dump all out to info for fun
  if o.verb==2 & exist('spikespy'), spikespy({Y,T{1},L{1},'unpert sort'},{F,'regen fwd model'},{Ynr,T{2},info.permL2(L{2}),'rev sort (best perm)'}); end    % show all signals
  if o.verb==3 & exist('spikespy'), spikespy({Y,T{1},L{1},'Y, unperturbed run'},{Ynr,T{2},info.permL2(L{2}),'Ynr, noise-reversed run'}); end    % show signals for paper
  % ... now could pull W's from 2nd run & compare re first?
  fsam = [fsam; info.acc.p];    % only works if K fixed. todo: fix

elseif strcmp(o.meth,'add')  % ================== self spike addition
  if ~isfield(o,'T')
    fprintf('ss spike-addition unperturbed run\n')
    [T L] = alg(Y);
  else
    fprintf('spike-addition overriding unpert ss run with input data\n')
    T = o.T; L = o.L;   % skip the unpert run (advanced)
  end
  K = max(L);
  wf = pull_waveforms_from_tseries(Y,T,L,o.Nt,o);
  pops = histc(L,1:K); rates = pops/N;    % found rates per time sample unit
  for r=1:Nr, fprintf('ss add run %d\n',r)
    o.rateunits = 's';     % since don't have known wf.d.dt or wf.d.samplefreq
    o.truePois = 1;        % allow spike # variation, according to fixed rates
    [F p] = synth_Poissonspiketrain(wf,N,o.ratescale*rates,[],o.Nt/2,[],o);
    newpops = histc(p.l,1:K), info.newpops = newpops;
    Ya = Y + F;            % linearly superpose new spikes
    [Ta La] = alg(Ya);     % rerun
    Te = [T p.t]; Le = [L p.l];   % expected times and labels for the rerun
    [Q tmiss tfals twrng info.permLa{r}] = times_labels_confusion_matrix(Te,Le,Ta,La,o);
    Qind = Q - diag([pops 0])    % induced confusion
    if 0                   % one version of metric
      tI = Qind(1:K,end)./newpops(:);  % col
      tII = Qind(end,1:K)./newpops;    % row
      disp('type-I error rates:'); fprintf('%.2g\t',tI); fprintf('\n')
      disp('type-II error rates:'); fprintf('%.2g\t',tII); fprintf('\n')
      f = 1 - tI(:)' - tII(:)'        % metric f_k = sum of error rates, kth unit
    else
      f = diag(Qind(1:K,1:K))./newpops'; f = f(:)';  % simple frac of new correct
    end
    if o.verb>3 & exist('spikespy'), spikespy({Y,T,L,'Y, unperturbed run'},{Ya,Ta,info.permLa{r}(La),'spike-addition run'}); end    % show signals for paper
    fsam = [fsam; f];    % only works if K fixed. todo: fix
    info.Qs = {info.Qs{:} Qind};
  end

end                           % =========== done validation algs

fhat = mean(fsam,1);  % best estimate of stability is mean
if iscell(L), L = L{1}; T = T{1}; end
info.pops = histc(L,1:K); info.L = L; info.T = T; % info out for 1st run only
info.wf = wf;

% show report card...
if o.verb
  disp('mean stabilities per spike type:'), fprintf('%6.3f ',fhat), fprintf('\n')
end
if o.verb>1, show_stabilities(fhat,fsam,o); end
