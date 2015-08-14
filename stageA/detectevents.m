function [X t m info] = detectevents(d,opts)
% DETECTEVENTS - detect single and multple spiking events on a small # channels
%
% [X t m] = detectevents(d) where d is an extracellular data struct, extracts
%  believed single spiking events in s, with timeshifts t, and
%  multiple-spiking events in a cell array m.
%  Small number of channels only, since detects across all channels.
%
% [X t m info] = detectevents(d,opts) controls certain options and outputs info.
%   opts.verb = 0 (silent), 1 diagnostic output (default)
%   opts.Twin = single-event window time in sec (default 0.003)
%   opts.thresh = threshold for detection signal (default auto-set)
%   opts.T1spike = time (in sec) that single spike can be above-threshold.
%   opts.sig = detection signal : 'a' = "energy" u^+(tau.u_t)^2
%                                 'm' = minimum over channels
%                                 'M' = max abs over channels
%
% Outputs:
% X - 3d array (M*Nt*Ns) of purported single-spike events (same as Jeremy)
%      Nt is number of window samples, same for all events. Ns = # spike events
% t - vector (1*Ns) time offsets of starts of single-spike windows (in 0 to N-1)
% m - struct for "multiple" events:
%     m.t - (1*Nm) start time offset (0-indexed) in samples, 0 to N-1
%     m.Ts - (1*Nm) lengths of each window in samples
%     m.X - concatenated M*m.Ttot signal array (m.Ttot = sum(m.Ts))
%     m.tptr - (1*Nm) column pointers to (1-indexed) start of each event in m.X
%     m.Ns - stores Nm

% Barnett 11/17/14
% 12/4/14: tidied wj0,wj1 detection, fixed spill-over end of 1:N samples
% 2/10/15: opts for plain thresh, channelwise-min. Changed multiple m format
% 2/19/15: fixed offset by 1 bug in m.X
% 3/12/15: o.sig default now 'm'
% 6/4/15: 'M' added. 6/8/15: changed t outputs to 0-indexed
% 6/11/15: thresh auto-set, quan killed, fixed short-clip 1-offset bug

if nargin<1, test_detectevents; return; end

if nargin<2, opts = []; end
if ~isfield(opts,'verb'), opts.verb = 1; end
if ~isfield(opts,'Twin'), opts.Twin = 0.003; end % width of a single-event window (s), controls Nt the clip lengths
if ~isfield(opts,'sig'), opts.sig = 'm'; end
if ~isfield(opts,'T1spike'), opts.T1spike = 0.001; end  % max time a single spike can be above-threshold (sec)

Nt = ceil(opts.Twin * d.samplefreq);  % # samples in single-event window
[M N] = size(d.A); 

% make one-channel signal for detection...
if strcmp(opts.sig,'a')  % Alex's method: sqrt(energy)
  Ap = [diff(d.A,[],2) zeros(M,1)]; % upwind stencil derivative
  tau = 0.0002;    % spike rise timescale (sec); hopefully dataset-independent
  sig = sqrt(sum(d.A.^2 + (tau*d.samplefreq)^2*Ap.^2,1)); % energy: u^2 + (tau.u')^2
  %figure; hist(sig,50); figure; semilogy((1:N)*d.dt,sig)
elseif strcmp(opts.sig,'m')
 sig = abs(min(d.A,[],1));     % abs of minimum over channels
elseif strcmp(opts.sig,'M')
 sig = max(abs(d.A,[],1));     % max of abs over channels
end

% Choose threshold for this one-channel signal... (ok except for 'a')
if isfield(opts,'thresh'), thresh = opts.thresh;
else, thresh = autothreshold(d); end
if opts.verb, fprintf('threshold = %.3g\n',thresh); end

% fatten the detection peaks to include Nt/2 either side, choose centering...
b = [0, conv(double(sig>thresh), ones(1,Nt),'full'), 0]; % 0-padded spread peaks
sh = -floor(Nt/2);       % shift, adjust to center window (to nearest sample)
b = diff(b>0); wj0 = find(b==1)+sh; wj1 = find(b==-1)+sh; % inds starts & stops
t = round((wj0+wj1)/2 - Nt/2)-1; % time offsets (in samples 0...N-1) if single
falloff = (t<0 | t+Nt-1>=N);  % true for windows who'll fall off even if short
t = t(~falloff);          % kill them
wj0 = wj0(~falloff); wj1 = wj1(~falloff); % (for wj0 could leave <1 or wj1 >N)
info.b=b; info.sig=sig; info.thresh=thresh; info.wj0=wj0; info.wj1=wj1; % diagn
Ne = numel(wj0); % total # events

% use event durations to split into single- vs multiple-spike events
sing = (wj1-wj0)*d.dt < (opts.Twin+opts.T1spike);     % duration criterion
js = find(sing); t = t(js); % "single" inds in window list, their start times
Ns = numel(js);  % # purported single-spike events
Nm = Ne-Ns;    % # multi-spike events (all are kept)
if opts.verb, fprintf('detectevents: Nt=%d, %d single-events & %d multi-events\n',Nt,Ns,Nm); end
info.Nt = Nt;

X = nan(M,Nt,Ns);
for n=1:Ns                   % loop over single events
  j = t(n) + (1:Nt);         % sample indices in window (t is 0-indexed)
  X(:,:,n) = d.A(:,j);       % copy event into output array
end

js = find(~sing); % probably multiple-spike events...
m.t = max(0,wj0(js)-1);            % 0-indexed time offset of each window
m.tlast = min(N-1,wj1(js)-2);      % last time offset in each window
m.Ts = m.tlast - m.t;              % lengths in signals
m.Ttot = sum(m.Ts);
m.Ns = numel(m.Ts);
m.tptr = cumsum([1 m.Ts(1:Nm-1)]); % pointers to start col indices in m.X
m.X = nan(M,m.Ttot);
for n=1:Nm, j=0:m.Ts(n)-1; m.X(:,m.tptr(n)+j) = d.A(:,m.t(n)+1+j); end % t 0-offset
%%%%%%%%%

function test_detectevents % makes a plot showing event windows
d = loaddemodata;
d.A = freqfilter(d.A,d.samplefreq,300,[]);
[s t m info] = detectevents(d);
show_detect(d,t,m,info);
set(gca,'xlim',[0 1]); % only 1st 1 sec for sanity
