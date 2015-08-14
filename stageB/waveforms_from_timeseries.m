function wf = waveforms_from_timeseries(A,samplefreq,o)
% WAVEFORMS_FROM_TIMESERIES  estimate waveforms from time series.
%
% wf = waveforms_from_timeseries(A,opts) returns a waveform object wf
%  (essentially a classifier) containing K waveforms, by detection of
%  candidate single spike events, upsampling, and clip-based spike-sorting
%  (not all clips need be classified).
%
% A = M-by-N filtered signal timeseries data
% samplefreq - sample rate in samples/sec
% opts (optional) controls parameters:
%  opts.thresh - size of negative detection threshold
%  opts.maxNclus - max # to send to clustering, noting time scales O(N^2) in this
%  opts.Kkeep - number of waveforms to keep from clustering (first Kkeep kept)
%  opts.eps - DBSCAN radius (sensitive)
%  opts.verb - verbosity 0,1,2...
%  opts.maxwid - collect random sample of size maxNclus only up to this width
%             (in sec). If zero, use narrowest; if Inf use random sample.
%  opts also passed to DETECTEVENTS, SPIKESORT_CLIPS
%
% See also: SPIKESORT_TIMESERIES which calls this

% Barnett 3/12/15. maxwid 4/20/15. changed name from class_filtereddata 6/11/15

if nargin<2, o = []; end
if ~isfield(o,'verb'), o.verb = 1; end
if ~isfield(o,'fmethod'), o.fmethod = 'pca'; end
if ~isfield(o,'num_fea'), o.num_fea = 10; end

% clus params
if ~isfield(o,'maxNclus'), o.maxNclus = 2e3; end  % since dbk is O(N^2)
if ~isfield(o,'maxwid'), o.maxwid = .0007; end  % spike width in sec
if ~isfield(o,'upsampfac'), o.upsampfac = 3; end
if ~isfield(o,'cmethod'), o.cmethod = 'dbpe'; end % fastest DBSCAN, auto eps
if ~isfield(o,'minpop'), o.minpop = 6; end

[M N] = size(A);
d.A = A; d.samplefreq = samplefreq; d.dt = 1/d.samplefreq; % fake an EC obj d
d.T = d.dt*size(d.A,2); clear A;

if o.verb, noi = empiricalnoise(d); fprintf('\tnoi.eta = %.3g\n',noi.eta), end
[X t m dinfo] = detectevents(d,o);  % threshold (no param overrides)
if ~isfield(o,'eps'), o.eps = 3.5*noi.eta*sqrt(o.num_fea); end  % try auto eps
[M Nt Ns] = size(X);
[tm tv] = signal_moments(X.^2); % means & variances about mean in time
maxtv = (o.maxwid*samplefreq)^2;
if Ns>o.maxNclus
  [~,i] = sort(tv);  % explore width filtering:
  if maxtv<tv(i(o.maxNclus))    % just use the narrowest spikes
    i = i(1:o.maxNclus);
    fprintf('using just %d narrowest spikes\n',o.maxNclus)
  elseif maxtv>tv(i(end))       % randomly sample from all spikes
    i = randperm(Ns, o.maxNclus);
  else            % random sample from spikes up to width maxwid
    Nnarrow = sum(tv(i)<=maxtv); i = i(randperm(Nnarrow,o.maxNclus));
    fprintf('using %d random of %d narrow spikes\n',o.maxNclus,Nnarrow)
  end
else, i = 1:Ns;  % use all
end
X = X(:,:,i);  % subset of clips to cluster
%o.ta = tm(i); % using 1st moment to align - was bad. Stick to min peak align

[L W z cinfo ta finfo] = spikesort_clips(X,o);           % the meat

if isfield(o,'ta'), o = rmfield(o,'ta'); end
K = size(W,3);  % initial # clusters found
pops = histc(L,1:K)
if o.verb
  fprintf('frac classified: %.3g (should be 0.5 to 0.9)\n',sum(pops)/size(X,3))
end

% decide which waveforms to keep...
if isfield(o,'Kkeep'), pok = 1:min(o.Kkeep,K); % force K if big enough (overrides minpop)
else, pok = pops>=o.minpop; end
wf.W = W(:,:,pok); % waveforms to keep
%ikeep = remove_shifted_duplicates(wf.W); % just for now
ikeep = pok;   % no filtering
if o.verb, fprintf('from K=%d, keep types: ',sum(pok>0)), fprintf('%3d ',ikeep), fprintf('\n'), end
%plot_spike_shapes(wf.W,'waveforms before ikeep'); drawnow;
wf.W = wf.W(:,:,ikeep);
wf.freqs = pops(pok(ikeep)); % these freqs are biased due width-selection
fprintf('making waveforms done, keeping K=%d\n',size(wf.W,3))
%if o.verb>1, plot_spike_shapes(wf.W,'waveforms from tseries: final W'); end

% package classifier struct (as a wf waveform struct, but with extra fields)
d = []; d.dt = 1/samplefreq; d.samplefreq=samplefreq;
d.name = sprintf('wf from unknown t-series of duration %d',N);
wf.d = d; wf.o = o; wf.fac = o.upsampfac;
wf.dinfo = dinfo; wf.finfo=finfo; wf.cinfo = cinfo;
