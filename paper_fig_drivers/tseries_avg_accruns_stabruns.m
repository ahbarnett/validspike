% averaging over synth-data runs for accuracy, and other runs for stability.
% Barnett 2/3/16

clear; clear all classes
load data_valid/wf_for_synth_matching_fig7b.mat    % load "wftrue" source wf
wftrue = wf; clear wf; K = size(wftrue.W,3);
%figure(1); plot_spike_shapes(wftrue.W,'truth W',[],[],0); drawnow

nsa = 1; %100;   % # indep synth data expts for acc
nss = 10;   % # indep synth data expts for stab (10 for add_r10, 50 for rev)

so = [];  % Synth opts -----
N = 2.4e6; %1e6;      % length in samples - should be 2.4e6 to match EJ 120s?
rates = wftrue.freqs;  % use the same rates as in data...
tpad = 10;           % end padding in samples
so.truePois = 1;
amplsig = 0.2;      % ampl variation they wanted: relative ampl std dev
so.ampl = amplsig;
noi.eta = 20;      % noise size, taken from fitting EJ (meth j)
noi.Nt = N; noi.M = size(wftrue.W,1);   % basic params for noise model

% Sorter opts ---------
% waveform extraction (clustering) opts... (happen to be same as used above)
co = []; co.cmethod='k++'; co.K = K; co.Kkeep = co.K; co.thresh=100; co.verb=1;
% fitting opts...
so = []; so.verb = 1; so.skip = 5; so.nlps = 10;
S = @(Y) spikesort_timeseries(Y,wftrue.d.samplefreq,co,[],so); % interface

% Accuracy opts --------
ao = []; ao.Nt = 60; ao.max_matching_offset = 10;

% ACCURACY RUNS =======================================
for i=1:nsa, fprintf('\nSYNTH REALIZATION FOR ACC #%d:\n\n',i)
  [Y ptrue] = synth_Poissonspiketrain(wftrue,N,rates,noi,tpad,[],so);
  fprintf('synth done, %d spikes\n',numel(ptrue.t));
  [T L] = S(Y);
  [perm{i} P{i} acc{i}] = times_labels_accuracy(ptrue.t,ptrue.l,T,L,ao);
end
% mean...
accm = 0*acc{1}.p; for i=1:nsa, accm = accm+acc{i}.p; end, accm = accm/nsa;
v = vertcat(acc{:}); accsam = vertcat(v.p);   % make nsa * K acc samples array
sso.blobcolor=[.9 .4 0]; show_stabilities(accm,accsam,sso); title(sprintf('ampl=%g: accs a_k',amplsig)); drawnow
%clear Y; save accsam_nsa30_N2400000_eta20_ampl0.mat


% Stability metric options ---------- num_runs = 10 (20 default for add)
o = []; o.Nt = 60; o.max_matching_offset=10;
o.meth = 'add'; o.ratescale = 0.25; o.num_runs = 10; o.verb = 5; % 5 lots info
%o.meth = 'rev'; o.verb = 0;

% STABILITY RUNS =======================================
for i=1:nss, fprintf('\nSYNTH REALIZATION FOR STAB #%d:\n\n',i)
  [Y ptrue] = synth_Poissonspiketrain(wftrue,N,rates,noi,tpad,[],so);
  fprintf('synth done, %d spikes\n',numel(ptrue.t));
  [fahat{i},fasam{i},info{i}] = eval_stability_tseriesbased(S, Y, o);
  [sperm{i} sP{i}] = times_labels_accuracy(ptrue.t,ptrue.l,info{i}.T,info{i}.L,ao);  % get the perm from 1st stab run to truth order (sP is truth-by-1ststabrun)
  [~,iperm] = sort(sperm{i});  % get ready to perm labels to truth order...
  Wp{i} = info{i}.wf.W(:,:,iperm); % 1st-run permed W
  fahat{i} = fahat{i}(iperm); fasam{i} = fasam{i}(:,iperm); % stab -> truth order
end  % recall info{i}.L not permed yet
% means...
fahatm = 0*fahat{1}; for i=1:nss, fahatm = fahatm+fahat{i}; end
fahatm = fahatm/nss;
show_stabilities(fahatm,vertcat(fahat{:})); title(sprintf('ampl=%g: stabs f_k',amplsig));

%clear Y; save data_valid/accstabsam_add_r10_nss10_N2400000_eta20_ampl0.mat
%save data_valid/accstabsam_add_wfap_N2400000_eta20_ampl02.mat

figure; plot(accm,fahatm,'.','markersize',20);
for k=1:K, text(accm(k),fahatm(k),sprintf('%d',k),'fontsize',20); end
xlabel('mean acc_k'); ylabel('mean f_k hat'); axis([0 1 min(fahatm) 1]);
title(sprintf('ampl=%g',amplsig));
