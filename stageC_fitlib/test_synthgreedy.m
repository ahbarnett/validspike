% Test multi-spike clip-based greedy fitting accuracy on synth clips.
% Barnett 2/28/15, based on test_multifitgreedy.m. 6/10/15 use default waveforms
%
% todo: update to use synth_Poissonspiketrain in each clip; std accuracy measure
clear;
wf = loaddefaultwaveforms; [M,T,K] = size(wf.W); d=wf.d;   % setup wf, # channels
sc = 1.5*max(abs(wf.W(:)));  % for plotting only
Nt = 40;   % Tc for fitting
tpad = 4;   % edge padding for alignment times

maxNs = 6; etafit = 18; etagen = 25; % make a bit larger. from EJ data
Is = 1:5;  % # spikes per clip to test
Nc = 1e3;  % # clips (was 1e5 for prelimreport)
msI = Is; fpI = Is;
Y = nan(M,Nc*Nt);
noi = setup_noisemodel(d,Nt,etagen,0.0002); % for generation
for i=1:numel(Is), I = Is(i); % ----------- loop over how many spikes
  pe = []; Nse=I*ones(1,Nc);
  Nsetot = sum(Nse); fprintf('I=%d spikes per clip: generating...\n',I);
  for c=1:Nc
    pe(c).l = randi(K,1,Nse(c)); pe(c).t = tpad+(Nt-2*tpad)*rand(1,Nse(c)); pe(c).a = ones(1,Nse(c)); % labels, times, ampls
    Y1 = spikemod(wf, pe(c), Nt);
    Y(:,Nt*(c-1)+(1:Nt)) = Y1 + noisesample(noi); % make Y
  end
  noi = setup_noisemodel(d,Nt,etafit);  % fitting (can't use eta=0!)
  Tc = Nt*ones(1,Nc);   % all same length
  fprintf('fitting... '); tic
  [p Ns Jbest info R] = multifitgreedy(wf,Y,Tc,noi,maxNs);
  fprintf('done in %.3f s: %.3g clips/s (%.3g spikes/s)\n',toc,Nc/toc,sum(Ns)/toc);
  missed = nan(1,Nc); falsepos = missed;
  o.terr = 4;  % matching time error allowed (0.2 ms)
  for c=1:Nc
    [~,~,ii] = spikesetmatch(pe(c),p(c),o);  % count # spikes matched
    missed(c) = numel(ii.pjmiss); falsepos(c) = numel(ii.qjmiss);
  end
  msI(i) = sum(missed)/Nsetot; fpI(i) = sum(falsepos)/Nsetot;
  fprintf('frac true spikes matched = %.4g\n',1-msI(i))
  fprintf('frac false positive spikes = %.4g\n',fpI(i));
  fprintf('reduction from mean J0=%.3f to mean Jbest=%.3f\n',mean(info.Jhist(1,:)),mean(Jbest))
  %fprintf('histogram of Ns...  '), histc(Ns,0:maxNs)
end

'summary of missed and false positive rates vs I:'
format short g
[Is; msI; fpI]
