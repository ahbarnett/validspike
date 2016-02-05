% averaging over synth-data runs of comparing accuracy to stability.
% Barnett 2/3/16

clear; close all
load data_valid/wf_for_synth_matching_fig7b.mat    % load "wftrue" source wf
wftrue = wf; clear wf; K = size(wftrue.W,3);
figure(1); plot_spike_shapes(wftrue.W,'truth W',[],[],0); drawnow

ns = 10;   % # indep synth data expts

so = [];  % Synth opts -----
N = 2.4e6; %1e6;      % length in samples - should be 2.4e6 to match EJ 120s?
rates = wftrue.freqs;  % use the same rates as in data...
tpad = 10;           % end padding in samples
so.truePois = 1;
so.ampl = 0;     % ampl variation they wanted: relative ampl std dev
noi.eta = 20;      % noise size, taken from fitting EJ (meth j)
noi.Nt = N; noi.M = size(wftrue.W,1);   % basic params for noise model

% Sorter opts ---------
% waveform extraction (clustering) opts... (happen to be same as used above)
co = []; co.cmethod='k++'; co.K = K; co.Kkeep = co.K; co.thresh=100; co.verb=1;
% fitting opts...
so = []; so.verb = 1; so.skip = 5; so.nlps = 10;
S = @(Y) spikesort_timeseries(Y,wftrue.d.samplefreq,co,[],so); % interface

% Stability metric options ---------- (verb=3 allows paper outputs)
o = []; o.Nt = 60; o.max_matching_offset=10;
o.meth = 'add'; o.ratescale = 0.25; o.num_runs = 20; o.verb = 3;
%o.meth = 'rev'; o.verb = 1;

% Accuracy opts --------
ao = []; ao.Nt = 60; ao.max_matching_offset = 10;

for i=1:ns, fprintf('\nSRUN #%d:\n\n',i) % =========== main loop over synth runs
  [Y ptrue] = synth_Poissonspiketrain(wftrue,N,rates,noi,tpad,[],so);
  fprintf('synth done, %d spikes\n',numel(ptrue.t));
  [fahat{i},fasam{i},infoa{i}] = eval_stability_tseriesbased(S, Y, o);
  [perm{i} P{i} acc{i} times{i}] = times_labels_accuracy(ptrue.t,ptrue.l,infoa{i}.T,infoa{i}.L,ao);  % accuracy labeled in truth order
  [~,iperm] = sort(perm{i}); Wp = infoa{i}.wf.W(:,:,iperm);
  figure(i+1); set(gcf,'position',[1000,250*(i-1),1000,250]);
  subplot(1,4,1); plot_spike_shapes(Wp,sprintf('srun %d: 1st stab run W permed to truth',i),[],[],0); drawnow
  if strcmp(o.meth,'rev'), Ts=infoa{i}.Tnr; Ls=infoa{i}.Lnr;  % get pert spikes
    Ls = perm{i}(infoa{i}.permL2(Ls));  % best perm of 1st run, permed to truth
    wfsp = pull_waveforms_from_tseries(infoa{i}.Ynr,Ts,Ls,o.Nt,o);
    %else, Ts=infoa{i}.Ta; Ls=infoa{i}.La; end
    subplot(1,4,2); plot_spike_shapes(wfsp.W,sprintf('srun %d: 2nd stab run W permed to truth',i),[],[],0); drawnow
  end
  fahat{i} = fahat{i}(iperm); fasam{i} = fasam{i}(:,iperm); % stab -> truth order
  sso=[]; sso.ylab='f_k'; sso.fig=0;        % start stability vs accuracy fig...
  subplot(1,4,3); show_stabilities(fahat{i},fasam{i},sso); hold on;
  sso.blobcolor=[.9 .4 0]; sso.fig = 0; % overlay accuracies in orange
  show_stabilities(acc{i}.p,[],sso);
  sso=[]; sso.fig = 0; show_stabilities(fahat{i},[],sso);  % re-overlay black blobs
  title(sprintf('srun %d: orange = acc, black = stability', i)); drawnow
  subplot(1,4,4); plot(acc{i}.p,fahat{i},'.','markersize',20);
  for k=1:K, text(acc{i}.p(k),fahat{i}(k),sprintf('%d',k),'fontsize',20); end
  xlabel('acc'); ylabel('f_k hat'); axis([0 1 min(fahat{i}) 1]); drawnow
end                                     % ============
% compute means over sruns
fahatm = 0*fahat{1}; accm = 0*acc{1}.p;
for i=1:ns, fahatm = fahatm+fahat{i}; accm = accm+acc{i}.p; end
fahatm = fahatm/ns; accm = accm/ns;
% ... try medians too
figure(ns+2); plot(accm,fahatm,'.','markersize',20);
for k=1:K, text(accm(k),fahatm(k),sprintf('%d',k),'fontsize',20); end
xlabel('mean acc_k'); ylabel('mean f_k hat'); axis([0 1 min(fahatm) 1]);
%title(sprintf('%d sruns, %d reps, amplsig=%.g\n',ns,o.num_runs,so.ampl));
hold on; plot([0 1],[0 1],'r-'); drawnow
%hgsave mean_acc_vs_stab_N2400000_ns2_add_r3_ampl0.fig
%clear Y; save mean_acc_vs_stab_ns10_addr20_ampl0.mat

figure(ns+3); cs = rand(K,3); % plot all runs
for k=1:K, for i=1:ns, text(acc{i}.p(k),fahat{i}(k),sprintf('%d',k),'fontsize',15,'color', cs(k,:)); end, end
xlabel('acc_k'); ylabel('f_k hat');
