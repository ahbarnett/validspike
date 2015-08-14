% test multifitgreedy on synthetic multi-spike clips
% Tweaked from ../test_fitgreedyspikes.m
% Barnett 2/16/15. defaultwaveforms 6/8/15

clear; wf=loaddefaultwaveforms; [M,T,K] = size(wf.W); d = wf.d; % setup wf
Nt = 40;       %30; overall clip length
sc = 1.5*max(abs(wf.W(:)));
pe.t = [13.35 10.5 20]; pe.l = [2 1 3]; pe.a=[1 1 1]; % pick times & identities
%pe.t = [13.35]; pe.l = [2]; pe.a=[1]; % pick times & identities

if 0  % noise-free one clip test
  Y = spikemod(wf, pe, Nt);
  figure; plot(Y', '-','color', .5*[1 1 1], 'linewidth',2);
  noi = setup_noisemodel(d,Nt,0.1);
  maxNs = 5;
  Nc = 1;    % one clip
  Tc = Nt;
  [p Ns Jbest info R] = multifitgreedy(wf,Y,Tc,noi,maxNs);
  f = Y-R; % get model outp
  hold on; plot(f', '-'); drawnow;
  fprintf('no noise: Jbest=%.3f. Params:\n',Jbest); pe, p
  fprintf('success = %d\n',spikesetmatch(pe,p))
  info.Jhist
end

if 1   % do multiple same-size clips w/ noise, to assess performance
  samep = 0;   % 1 to repeat above true params, 0 to use new rand ones
  maxNs = 6; etafit = 18; etagen = 25; % make a bit larger. from EJ data
  Nc = 1e4;  % # clips (was 1e5)
  Y = nan(M,Nc*Nt);
  noi = setup_noisemodel(d,Nt,etagen,0.0002); % for generation
  if samep
    Y1 = spikemod(wf, pe, Nt); Nsetot = numel(pe.t)*Nc;
    for c=1:Nc, Y(:,Nt*(c-1)+(1:Nt)) = Y1 + noisesample(noi); end % make Y
  else
    I = 5;   % how many spikes per clip, fixed
    pe = []; Nse=I*ones(1,Nc); %Nse = randi(I,1,Nc);% how many spikes (up to I)
    Nsetot = sum(Nse);
    for c=1:Nc
      pe(c).l = randi(K,1,Nse(c)); pe(c).t = 5+20*rand(1,Nse(c)); pe(c).a = ones(1,Nse(c)); % labels, times, ampls
      Y1 = spikemod(wf, pe(c), Nt);
      Y(:,Nt*(c-1)+(1:Nt)) = Y1 + noisesample(noi); % make Y
    end    
  end
  noi = setup_noisemodel(d,Nt,etafit);  % fitting (can't use eta=0!)
  Tc = Nt*ones(1,Nc);   % all same length
  fprintf('fitting... '); tic
  omfg=[]; %omfg.locflag = 0;   % fitting opts (for short clips 0 is best)
  [p Ns Jbest info R] = multifitgreedy(wf,Y,Tc,noi,maxNs,omfg);
  fprintf('done in %.3f s: %.3g clips/s (%.3g spikes/s)\n',toc,Nc/toc,sum(Ns)/toc);
  suc = nan(1,Nc);  % # successful spikes in each clip
  o.terr = 2;  % matching time error allowed
  for c=1:Nc
    if samep, pec = pe; else, pec = pe(c); end  % exact params
    [~,suc(c)] = spikesetmatch(pec,p(c),o);       % count # spikes matched
  end
  fprintf('fraction of true spikes matched = %.4g\n',sum(suc)/Nsetot)
  fprintf('reduction from mean J0=%.3f to mean Jbest=%.3f\n',mean(info.Jhist(1,:)),mean(Jbest))
  fprintf('histogram of Ns...  ')
  histc(Ns,0:maxNs)
  plot_spike_shapes(wf.W,'W',sc);
  n=50; plot_spike_shapes(reshape(Y(:,1:n*Nt),[M Nt n]),'Y',sc);
  plot_spike_shapes(reshape(R(:,1:n*Nt),[M Nt n]),'R',sc);
  figure; plot(1:n, [info.Jhist(1,1:n); Jbest(1:n)], '+');
  title('J0 and Jb vs c');
  %plot_spike_shapes(Y,'Y',sc); plot_spike_shapes(R,'R',sc) % all
  if 0, c = 35; pe(c),p(c) % check out one clip
    Yc = Y(:,Nt*(c-1)+(1:Nt)); plot_spike_shapes(Yc,'Y_c',sc);
    Rc = R(:,Nt*(c-1)+(1:Nt)); plot_spike_shapes(Rc,'R_c',sc);
    o.verb=2; fitgreedyspikes(wf,Yc,noi,maxNs,o)
  end
end

if 0  % single-clip reliability exploration
  o = []; pec=[]; pec.l=[2 2 2]; pec.t = [6 13 22]; pec.a = [1 1 1] % nice fail!
  noi = setup_noisemodel(d,Nt,0.1,0.0002);
  Yc = spikemod(wf,pec,Nt); Yc = Yc + noisesample(noi);
  plot_spike_shapes(Yc,'Y_c',1.2);
  pc = fitgreedyspikes(wf,Yc,noi,maxNs,o)  % can have up to 5 spikes, many l=7
  plot_spike_shapes(spikemod(wf,pc,Nt),'F_c',1.2,pc.t) % show model
  pt=[]; pt.l=2; pt.t=6; pt.a=1; % but if feed it correct 1st spike, success...
  pc = fitgreedyspikes(wf,Yc-spikemod(wf,pt,Nt),noi,maxNs,o)
end

% todo: debug why not counting as successes. output p from spikefit.


% laptop beth timings: (Nt=30)
% For large sets (Nc 1e4 or more):
%   Ns=3: 700 clips/sec one core i7; 3000 clips/sec on 4 cores i7 (8 threads)
%   Ns=1: 6400 clips/sec on 4 cores i7.
%   Ns=3, Ns=1e5: 2800 clips/s (4 core i7, 8 threads)
% Once R requested for Nc 1e3 or less, only
%   Ns=3: 1700 clips/s; Ns=1; 3500 clips/s, ie half-speed!
% Mystery how get 800% CPU (hyperthreading) - since caching keeps in each thread
% ??
%  samep=0: 2300 clips/s i7

% desktop scda003 timings:
% Ns=3:
% samep=0: 

