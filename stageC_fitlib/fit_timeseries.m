function [t l R] = fit_timeseries(A,C,meth,o)
% FIT_TIMESERIES  spike times and labels from timeseries given waveforms (MEX)
%
% [t l] = fit_timeseries(A,wf,meth,o) returns times and labels given filtered
%  timeseries data and classifier by fitting waveforms via one of several
%  methods controlled by string meth.
%
% [t l R] = fit_timeseries(A,wf,meth,o) also returns the model residual of
%   same dimensions as A.
%
% Inputs:
%  A = M-by-N filtered signal data
%  wf = waveform object, a classifier from class_filtereddata()
%  meth - text string choosing fitting method (all involve MEX interfaces):
%         'g' greedy on clips
%         'gs'/'gl'/'Gl' glutton on whole dataset (fastest is 'Gl')
%  o (optional struct) controls options. General options:
%   o.noi - provides noise model object (otherwise it is fit from data)
%   o.nlps - neg log priors (lambda_l) for each type
%
%  Options specific for each method:
%  meth = 'g':
%   o.thresh - detection threshold, overrides C.o.thresh
%   o.maxspikesperclip - max # spikes fit per clip
%   o.Twin - fitting window length in sec
%  meth = 'gs', 'gl', 'Gl':
%   o.gamma - safe time interval for gamma-local min (in samples)
%   o.skip - integer coarsening factor for S score eval 
%
% Outputs:
%  t,l = 1-by-Ns arrays of real-valued times (in sample units starting at 0),
%        and integer labels in 1...K where K=size(C.W,3)
%  R = (optional) residual signal (same size as A) after subtraction of found
%      spikes
%
% Called without arguments, a self-test is performed on synthetic data
%
% Notes: 1) no amplitude fitting is done
%        2) residual R returned is computed in various ways that may not exactly
%           equal the residual computed from A-F where F is forward model, eg
%           F = spikemod(wf,p,size(A,1)), where p.t=t and p.l=l are fitted params
%
%  See also: SPIKESORT_TIMESERIES

% Barnett 3/12/15. 3/19/15, meth interface & glutton 5/20/15
% 6/2/15 name change & self-test; 6/8/15 0-offset for times
% 7/15/15 renamed from sort_timeseries

if nargin<1, test_fit_timeseries; return; end
if nargin<4, o = []; end
if ~isfield(o,'verb'), o.verb = 1; end  % some method defaults...
fac = C.fac; d = C.d; d.A = A; clear A; [M N]=size(d.A); % fake EC data struct
[M T K] = size(C.W);
if ~isfield(o,'nlps'), o.nlps=10*ones(1,K); end
if ~isfield(o,'gamma'), o.gamma = 0.0005*d.samplefreq; end
if ~isfield(o,'skip'), o.skip = 5; end
if isfield(o,'noi'), noi=o.noi; else noi = empiricalnoise(d); end  % noise model
if o.verb, fprintf('fit_timeseries, meth=%s...\n',meth), end
t1 = tic;

if strcmp(meth,'g')   % ........................ greedy clip-based fit

  if ~isfield(o,'thresh'), o.thresh = C.o.thresh; end
  if ~isfield(o,'maxspikesperclip'), o.maxspikesperclip = inf; end
  if ~isfield(o,'Twin'), o.Twin = .002; end  % fitting window
  
  [X t m dinfo] = detectevents(d,o); [M Nt Ns] = size(X);
  if o.verb>1, show_detect(d,t,m,dinfo); end
  if numel(t)==0, error('no events detected; stopping!'); end

  t = [t m.t]; m = mergeclips(X,m); % concatenate short & long clips
  Nshort = Ns; Ns = m.Ns;
  maxIpert = 0.1;  % max # spikes per unit time sample (not important)
  maxI=min(ceil(max(m.Ts)*maxIpert),o.maxspikesperclip);   % max spikes per clip (sets alloc size)

  o.tpad = 5.0;
  if nargout>2   % compute full signal resid
    r = m; r.X = [];    % copy time info to resid clip struct
    [p I Jbest info r.X] = multifitgreedy(C,m.X,m.Ts,noi,maxI,o); % clip resids
    R = d.A;     % paste resid from clips into full timeseries...
    for c=1:Ns, j = 0:r.Ts(c)-1; R(:,t(c)+j+1) = r.X(:,r.tptr(c)+j); end
  else           % no resid needed
    [p I Jbest info] = multifitgreedy(C,m.X,m.Ts,noi,maxI,o);
  end
  fprintf('# clips fit with 0 spikes: %d\n',sum(I==0))
  
  l = [p(:).l];     % stack overall l,t arrays
  for c=1:Ns, p(c).tabs = p(c).t + t(c); end  % convert p.t (rel) to abs time
  t = [p(:).tabs];

  
elseif strcmp(meth,'gs') % .....glutton whole dataset, slow (1-core, fresh S's)
  
  tpad = 2; tsh = tpad + (0:floor((N-2*tpad)*fac))/fac;  % time shifts for S
  newspikes = 1; nfit = 1; pp.t = []; pp.l = []; % store all spikes
  R{1} = d.A;  % keep residuals
  while newspikes   % ....... repeated fitting
    tic; S = fillscore(C,R{nfit},tsh,noi,o);
    if o.verb, fprintf('fit %d: S done in %.3g s = %.3g G pts/sec\n', nfit, ...
                       toc, 1e-9*numel(tsh)*M*K*T/fac/toc), end
    tic; [jt l s] = locvalidmins(S,o.nlps); Nm = numel(jt); % get all LVMs
    if o.verb, fprintf('%d LVMs found in %.3g s\n',Nm,toc), end
    [jt l s Nm] = keepgamlocmins(jt,l,s,tsh,o.gamma);
    p.l = l; p.t = tsh(jt);  % time shifts of LVMs (jt 1-indexed)
    if o.verb>1, figure;plot(tsh,S,'-');hold on;plot(p.t,s,'.','markersize',20);
      for i=1:numel(jt), text(p.t(i),s(i)-50,sprintf('%d',l(i))); end
      for i=1:numel(pe.t), text(pe.t(i),0,sprintf('%d',pe.l(i))); end
    end
    R{nfit+1} = R{nfit} - spikemod(C,p,N); % for now - subtracts a lot of zeros
    pp.t = [pp.t p.t]; pp.l = [pp.l p.l]; % append
    pf{nfit} = p;  % store each fit
    nfit = nfit+1; newspikes = numel(p.t)>0;
  end
  R = R{end}; t = pp.t; l = pp.l;
  
  
elseif strcmp(meth,'gl')   % ........ glutton whole dataset (1 core, faster)

  [t l R] = gluttonstuffme(C,d.A,noi,o.nlps,o);
  
  
elseif strcmp(meth,'Gl')   % ........ glutton whole dataset (multicore, fastest)

  [t l] = multiglutton(C,d.A,noi,o.nlps,o);
  if nargout>2, 
    F = spikemod(C,struct('t',t,'l',l),N); R = d.A - F;  % brute-force residual
  end
  
end

[t,i] = sort(t); l=l(i); % time-order the output spikes
if o.verb
  fprintf('fit_timeseries found %d spikes in %.3g s\npops: ',numel(t),toc(t1))
  fprintf('%6d ',histc(l,1:max(l))), fprintf('\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function test_fit_timeseries  % synthesis-based accuracy and speed test

wf = loaddefaultwaveforms; [M T K] = size(wf.W); fs = wf.d.samplefreq;
N = round(1.0*fs);         % 1 second of time
noi = setup_noisemodel(wf.d,N,25);  % choose noise std dev eta
firingrate = 100; rates = firingrate*ones(1,K);            % mean rates in Hz
[Y pe] = synth_Poissonspiketrain(wf,N,rates,noi,[],0);
fprintf('Y range = [%.3g,%.3g]\n',min(Y(:)),max(Y(:)))
o.nlps = log(wf.fac*fs./rates);  % prior = the true rates

meths = {'g','gs','gl','Gl'};
o.verb = 1; o.skip = 6; o.gamma = 0.0006*fs; % some stageC alg params
for i=1:numel(meths)
  fprintf('------ meth=%s :\n',meths{i}); t1=tic;
  [p.t p.l R] = fit_timeseries(Y,wf,meths{i},o);
  fprintf('------ done: %d spikes, %.3g s\n',numel(p.t),toc(t1))
  F = spikemod(wf,p,N); Re = Y - F;          % naive compute residual
  fprintf('             resid rms %.4g\n',sqrt(sum(Re(:).^2)/M/N))
  fprintf('             rms err vs fwd-model resid %.3g\n',norm(R(:)-Re(:))/sqrt(M*N))
  [Q tfals tmiss twrng]=times_labels_confusion_matrix(p.t,p.l,pe.t,pe.l);
  fprintf('\n')
  % todo: don't forget could 3000-Hz low-pass filter R here to find bad spots...
end

if exist('spikespy')    % Just show the last method's sorting results...
  spikespy({Y,pe.t,pe.l,'synth data Y & true spikes'},{F,p.t,p.l,'F & fitted spikes'},{R,[tfals,tmiss,twrng],[1+0*tfals,2+0*tmiss,3+0*twrng],'R + falsepos(1),miss(2),wrong(3)'});
end
