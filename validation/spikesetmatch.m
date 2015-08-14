function [b n info] = spikesetmatch(p,q,o)
% SPIKESETMATCH - decide if two spike parameter sets are close enough to match
%
% [b n info] = spikesetmatch(p,q) returns b=true if p and q match all spike
%  parameters sufficiently closely, and n the number of spikes matching
%  sufficiently closely.  p and q both have fields t, l, a, each a
%  row-vector of time-offsets, identities and amplitudes respectively.
%  A match means number of spikes and identities match exactly, and times
%  and amplitudes are within terr and aerr respectively. A greedy algorithm
%  is used, O(I^2) where I = number of spikes.
%
% Outputs:
%  b    - true if complete match, false otherwise
%  n    - integer number of spikes matching
%  info - struct with info about type of (mis)match:
%         info.i - which index in q matches each index in p (ie permutation)
%         info.pjmiss - indices of p spikes unable to match
%         info.qjmiss - indices of q spikes unable to match
%
% [b,...] = spikesetmatch(p,q,o) controls options such as:
%  o.terr - set allowable time error (in sample units, default 1)
%  o.aerr - amplitude error (default Inf, a wide berth)
%
% [fails] = spikesetmatch; tests the routine, returns an array of
%  structs of failure situations
%
% Notes:
% 1) This is a legacy code used only by accuracy testers in stageC_fitlib.
%  It has been superceded by TIMES_LABELS_ACCURACY.
% 2) The greedy algorithm can occasionally fail to find a match that does
%  exist, if there are same-identity spikes with nearby times. This is likely
%  to be rare, esp. if o.terr is set at least a couple of times larger than the
%  time fitting error.

% Barnett 2/11/15 tweaked (m->l) from fitwindow01 of 3/13/14. 2/18/15: n output

if nargin<1, b = test_spikesetmatch; return; end
if nargin<3, o = []; end
if ~isfield(o,'terr'), o.terr = 1; end
if ~isfield(o,'aerr'), o.aerr = Inf; end
Ip = numel(p.t); Iq = numel(q.t);
P = true(1,Ip); Q = true(1,Iq);  % true if each p, q spike is still available
info.i = nan(1,Ip);
n=0;
for i=1:Ip      % loop through each p spike and remove a matching q spike
  t = q.t(Q);  % available times
  % j = indices of q spikes in available list only, that match...
  j = find(q.l(Q)==p.l(i) & abs(t-p.t(i))<=o.terr & abs(q.a(Q)-p.a(i))<=o.aerr);
  if ~isempty(j)         % we have at least one match
    [~,k] = min(abs(t(j)-p.t(i))); % k is index of matches matching time best
    fQ = find(Q); jQ = fQ(j(k));   % match only this best one
    info.i(i) = jQ; Q(jQ) = false; P(i) = false; % update perm, remove pair
    n = n+1;      % count the matches
  end
end
b = sum(P)+sum(Q)==0;  % Boolean output is: are all spikes accounted for?
info.pjmiss = find(P); info.qjmiss = find(Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fails] = test_spikesetmatch
% test driver for spikesetmatch. Output is struct array of failure situations;
% see above for field names

if 1, % warm-up tests
  p.t = [-5 .4 13]; p.l = [1 3 4]; p.a = [.9 1 1.1]; I = numel(p.t);
q = p;
[b n i] = spikesetmatch(p,q); if b,fprintf('ok\n'), else fprintf('bad\n'), end
q.t = q.t + 2*rand(size(q.t))-1; q.a = q.a + 0.1*(2*rand(size(q.a))-1); % jiggle
if spikesetmatch(p,q), fprintf('ok\n'), else, fprintf('bad\n'), end
a = randperm(I); q.t = q.t(a); q.l = q.l(a); q.a = q.a(a); % permute spikes
[b n i] = spikesetmatch(p,q); if b,fprintf('ok\n'), else fprintf('bad\n'), end
q.t= q.t(1:2); q.l = q.l(1:2); q.a = q.a(1:2); % wrong number of spikes
[b n i] = spikesetmatch(p,q); if ~b,fprintf('ok\n'), else fprintf('bad\n'), end
end

fprintf('please wait...\n')
tic, n0=1e4; M=10; Imax = 10; % how many tests to do, # neuron types, # spikes
dt=0.3; da=0.5; % rand jiggle sizes
nfail = 0; fails.p = []; fails.q = []; fails.pjmiss = []; fails.qjmiss = []; fails.i = [];
for n=1:n0
  I = randi(Imax);
  p.t = 10*rand(1,I); p.l = randi(M,1,I); p.a = rand(1,I); % p
  a = randperm(I); q = p; q.t = q.t(a); q.l = q.l(a); q.a = q.a(a); % q = perm p
  q.t = q.t + dt*(2*rand(size(q.t))-1); q.a = q.a + da*(2*rand(size(q.a))-1);
  [b n i] = spikesetmatch(p,q);  % test the code
  fail = ~b || max(abs(q.t(i.i)-p.t))>1; % short-circuit logical or
  if fail, nfail=nfail+1; fails.p = [fails.p p]; fails.q = [fails.q q];
    fails.pjmiss = [fails.pjmiss i.pjmiss];
    fails.qjmiss = [fails.qjmiss i.qjmiss]; fails.i = [fails.i i.i]; end
end
fprintf('%d fails out of %d random trials, %.3g s\n',nfail,n0, toc)

% Notes from March 2014:

% eg number of missed spikes could be summing numel(i.pjmiss) over trials

% NB greedy can fail occasionally, if ratio dt/o.terr is too large, due to
% misassociating based on the nearest times (for subset of same identity).
% Here we have ratio = 0.3.
% Expect failure rate v. small in real world due to refractory period.
