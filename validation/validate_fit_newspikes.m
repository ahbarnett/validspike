% script for inserting spikes to validate idea. Barnett 3/18/15

clear; d = loaddata('bb');
cm = 1;

A = freqfilter(d.A,d.samplefreq,300,[]); d.A = [];
% A is now filtered data; d only kept around for its dataset info

if ~cm
  if 0 % w/o mean subtraction
    %o.thresh = 150; o.eps = 850; o.maxNclus = 1000; o.K = 5; % 'b'
    o.thresh = 150; o.eps = 800; o.maxNclus = 1000; o.K = 5; % 'bb'
    C = class_filtereddata(A(:,1:1200000),o); % better since allow some wider spikes
    %set(gcf,'paperposition',[0 0 6 8]); print -depsc2 K5_clus.eps  
  else
    C = load('data/waveforms_2.5ms_b_fac3_tv_dbk_K5.mat');
  end
  sc = 500;
else    % common mode
  Acm = freqfilter(mean(A,1),d.samplefreq,[],3000); % common-mode signal
  A = bsxfun(@minus, A, Acm);  % kill common-mode
  if 0
    %o.thresh = 120; o.eps = 450; o.maxNclus = 1000; o.K = 7; % 'b'
    o.thresh = 120; o.eps = 450; o.maxNclus = 2000; o.K = 7; %'bb'
    %o.thresh = 140; o.eps = 400; o.maxNclus = 2500; o.K = 7; %'bb' messing
    C = class_filtereddata(A(:,1:1200000),o);
    %set(gcf,'paperposition',[0 0 6 6]); print -depsc2 K7_cm_clus.eps 
  else
    C = load('data/waveforms_2.5ms_b_fac3_tv_dbk_cm_K7.mat');
    %C = load('data/waveforms_2.5ms_b1_fac3_tv_dbk_cm_K7.mat'); % b1
  end
  sc = 400;
end
plot_spike_shapes(C.W,'W',sc); set(gcf,'position',[100 900 200 200]); drawnow

o  = []; tic; [t l R] = sort_filtereddata(A,C,o);
fprintf('raw sort found %d spikes in %.3f sec (%.3g spikes/s)\n',numel(t),toc,numel(t)/toc)
K = size(C.W,3); [M N] = size(A);
pops = histc(l, 1:K), Ns = numel(l);
F = A-R; ssview(A,F,R,{t l});

if 0  % finding overlaps - Break this expt out:
i = (l==1 | l==2 | l==5 | l==6); % find overlaps only within certain types
i = (l==5 | l==2); % find overlaps only within certain types
%i=l>0; % all firings
Tclose = .002; ll=l(i); tt = t(i); Ni = numel(ll); j=(diff(tt)<Tclose/d.dt);
sum(j)/Ni  % 208 of 4454 from types 1,2,5,6 are within 2 ms of next spike
% so 9% are within 2ms of another spike
Tclose/(d.T/Ni) % Poisson expect 0.06, so around same (bit under)

j=find(j);

j=(tt(3:end)-tt(1:end-2)<Tclose/d.dt); % 3 within Tclose?
sum(j)/Ni  % 1 from 4454 types 1,2,5,6 are 3 firings within 2 ms
(Tclose/(d.T/Ni))^2 % Poisson
end

if 1 % ADDING ALL TYPES AT ONCE...
if 0 % same # of each
  popsf = 100*ones(1,K);  % # fake spikes to add of each type
else % in proportion to existing
  popfrac = 0.25;  % what frac of original # spikes to add as new spikes
  popsf = ceil(popfrac*pops);
end
Ns = sum(popsf);
p.l = []; for k=1:K, p.l = [p.l, k*ones(1,popsf(k))]; end % build labels
p.t = N*rand(1,Ns); p.a = ones(1,Ns) + 0.05*randn(1,Ns); % t, ampl variation
%F = spikemod(C, p, N); ssview(F,{p.t+1,p.l}) % debug. must Ns<MAXSPIKESPEREVENT
Af = A + spikemod(C, p, N); % one big clip
te = [t p.t+1]; le = [l p.l]; % note time-shift of 1 sample since t 1-indexed
tic; [tf lf] = sort_filtereddata(Af,C,o);
fprintf('raw+(%d new spikes) sort found %d spikes in %.3f sec (%.3g spikes/s)\n',Ns,numel(tf),toc,numel(tf)/toc)
opts.max_matching_offset = 10; % +-0.5 ms
P = times_labels_confusion_matrix(te,le,tf,lf,opts);
Pnew = P - diag([pops 0])  % change in confusion mat
pmiss = Pnew(1:K,end)'./popsf
pfalspos = Pnew(end,1:K)./popsf
%reliperc = 100*(1 - (pmiss+pfalspos)/2);
reliperc = 100*(1 - (pmiss+pfalspos));
[sprintf('stabilities: \t'), sprintf('%.1f\t',reliperc)]
end  

if 0 % ADDING TYPES SEPARATELY...
nf = 100;   % # fake spikes to add of one type
P = zeros(K+1,K+1,K);
for k=1:K
  p.t = N*rand(1,nf); p.a = ones(1,nf); p.l = ones(1,nf)*k;
  Af = A + spikemod(C, p, N); % one big clip
  te = [t p.t]; le = [l p.l];
  [tf lf] = sort_filtereddata(Af,C,o);
  fprintf('k=%d sort found %d spikes\n',k,numel(tf))
  opts.max_matching_offset = 10; % 0.5 ms
  P(:,:,k) = times_labels_confusion_matrix(te,le,tf,lf,opts);
end
P
end
