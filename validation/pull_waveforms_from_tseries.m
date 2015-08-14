function wf = pull_waveforms_from_tseries(Y,T,L,Nt,o)
% PULL_WAVEFORMS_FROM_TSERIES  estimate upsampled W given times and labels and Y
%
%  wf = pull_waveforms_from_tseries(Y,T,L,Nt,opts)
%   returns wf a waveform struct with fields W (waveforms), fac (upsampling
%   factor).
%
%  Nt = window width to extract in samples (before upsampling which loses
%          some of the ends). Optional (default 60)
%  opts.upsampfac = (optional) upsampling factor (default 3)
%  opts.verb = verbosity (0 = no output, 1 = timing)
%
% Notes: 1) Plain averaging for now.
%  Could better be a linear system w/ MTK unknowns (1e4, too big!) and lots of
%  rows, but only matters if dominated by overlapping events.
% 2) the indexing of time alignments has been set to match that in spikemod
%  and alignspikes and upsample. All of these offsets matter, and are tricky.

% Barnett 7/16/15 based on old fac=1 extractmeanwaveforms

if nargin==0, test_pull_waveforms_from_tseries; return; end
if nargin<5, o = []; end
if ~isfield(o,'verb'), o.verb = 0; end
if ~isfield(o,'upsampfac'), o.upsampfac = 3; end  % needed so wf.fac always set
if ~isfield(o,'kerpars'), o.kerpars = []; end
if ~isfield(o.kerpars,'Tf'), o.kerpars.Tf=5*(o.upsampfac>1); end   % default

if o.verb, t1 = tic; fprintf('pull_W_from_tseries...'); end
K = max(L);
[M N] = size(Y);
keep = find(T>=Nt/2 & T<N-Nt/2); T = T(keep); L = L(keep);   % no fall off ends
Ns = numel(T);
intT = round(T);            % best sample to center on (0-indexed)
sh = -floor(Nt/2);          % shift, adjust to center window (to nearest sample)
X = zeros(M,Nt,Ns);         % allocate (not-upsampled-yet) clips
for n=1:Ns                  % loop over spikes
  j = intT(n) + sh + (1:Nt); % sample indices in time series
  X(:,:,n) = Y(:,j);        % copy event into output array
end
X = upsample(X,o.upsampfac,o.kerpars);
ush = floor(Nt/2) - o.kerpars.Tf+1;            % time shift of central index 
X = alignspikes(X,o.upsampfac,[],T-intT+ush);  % do fractional alignment
wf.W = meanwaveforms(X,L);                     % clip-based means
wf.fac = o.upsampfac;
if o.verb, fprintf(' done in %.3g s\n',toc(t1)); end
%%%%%


function test_pull_waveforms_from_tseries    % tests (used to check alignment)

if 1, disp('small test...')
  wf = loaddefaultwaveforms; K = size(wf.W,3);
  Ns = 100; p.l = randi(K,1,Ns); p.t = 100*(1:Ns) + 5*rand(1,Ns);
  N = max(p.t+100);  % total time
  Y = spikemod(wf,p,N);
  Nt = 60;
  wf2 = pull_waveforms_from_tseries(Y,p.t,p.l,Nt);
  % size(wf.W), size(wf2.W) % will match if Nt=60 since that made default wf
  F = spikemod(wf2,p,N);
  if exist('spikespy'), spikespy({Y,p.t,p.l,'tseries'},{F,p.t,p.l,'regen fwd model'},{F-Y,'diff'}); end
  % note interp errs at edges of waveforms
  fprintf('rel l2 err norm recon fwd model = %.3g\n',norm(F(:)-Y(:))/norm(F(:)))
  fprintf('rel l2 err norm W = %.3g\n',norm(wf.W(:)-wf2.W(:))/norm(wf.W(:)))
  %keyboard
end

if 1, disp('larger test using default noisy time series...')
  wf = loaddefaultwaveforms; % the waveforms that were used for synth
  d = loaddemodata;          % synthetic time series
  h = fileparts(mfilename('fullpath')); % current dir
  gndfile = fullfile(h,'../data/EC_default_synth_groundtruth.mat'); % "
  load(gndfile,'p');         % its ground truth
  Nt = 60;
  wf2 = pull_waveforms_from_tseries(d.A,p.t,p.l,Nt);
  K = size(wf2.W,3); sc = 400;
  ta = ones(1,K) * (size(wf2.W,2)-1)/2; % supposed gnd-truth alignments
  plot_spike_shapes(wf.W,'true',sc,ta); plot_spike_shapes(wf2.W,'pulled',sc,ta);
  fprintf('rel l2 err norm W = %.3g\n',norm(wf.W(:)-wf2.W(:))/norm(wf.W(:)))
  % expect this to have noise at 1/sqrt(# spikes) level
  %keyboard
end
