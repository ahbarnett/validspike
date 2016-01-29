% Starting the stability metric vs accuracy comparison for JNM resubmission.
% Barnett 11/20/15
% For Jeremy to edit..!

clear
%load data_valid/synth_accuracy_gndtruth_45swap_vs_paper.mat

if 0  % if true, extract wf from a fresh fit (see synth_tseries_as_paperfig7b)
  
% 1. Extract waveforms by doing spikesorting (20 sec laptop) ----------------
clear; d = loaddata('e');  % EJ 2005 elec359 timeseries from data_external/
d.A = freqfilter(d.A,d.samplefreq,300,[]);
% waveform extraction opts... (these are as in fig_tseries_meths.m)
co.cmethod='k++'; co.K = 10; co.Kkeep = 10; co.thresh=100; co.verb = 2;
% fitting opts...
so.skip = 5; so.nlps = 10;
[~,~,p,wftrue] = spikesort_timeseries(d.A,d.samplefreq,co,[],so);
wftrue.freqs = histc(p.l,1:co.Kkeep)/d.T;   % mean firing rates (Hz) of each label
noi = empiricalnoise(d);            % get noi.eta = noise std deviation.
% the point was to get a wf (wavefunc object) and noi (noise model); clear rest
clear p co so

else
  load data_valid/wf_for_synth_matching_fig7b.mat
  wftrue = wf; clear wf
  noi.eta = 18.3;      % matches fit
end

% 2. Create fresh synthetic data (post-frequency-filtered) -----------------
N = 1e6;             % length in samples - should be 2.4e6 to match EJ 120s?
rates = wftrue.freqs;  % use the same rates as in data...
%rates(9) = 0;   % ...but kill noise cluster? didn't rearrange wf yet
tpad = 10;           % end padding in samples
so = [];  % synth opts...
so.truePois = 1;
so.ampl = 0.2;     % *** firing variation they wanted: relative ampl std dev
% but we should vary more
noi.Nt = N; noi.M = size(wftrue.W,1);   % basic params for noise model
[Y ptrue] = synth_Poissonspiketrain(wftrue,N,rates,noi,tpad,[],so);
fprintf('synth done, %d spikes\n',numel(ptrue.t));
% *** make a code edited from synth_Poissonspiketrain if you want to change the ptrue.t timing synthesis, eg to make correlated firings, etc!
% Or, fill your own ptrue.t .l .a and just call
%    Y = spikemod(wf, ptrue, N);
%    Y = Y + noisesample(noi);
% to run fwd model and add iid noise
% *** adding correlated noise would be better, by doing instead
%    Y = Y + freqfilter(noisesample(noi),wftrue.d.samplefreq,300,[]));
% We can discuss this.

% NB only wftrue.d used from wftrue:
%save data_valid/synth_accuracy_gndtruth.mat Y ptrue so noi wftrue

plot_spike_shapes(wftrue.W,'true W');
%spikespy({Y,ptrue.t,ptrue.l,'synth'});   % show what we made
%show_crosscorr(ptrue.l,ptrue.t);  % no refractory holes, of course, since
                                  % purely Poisson spike trains for now


% Setup inline sorter function S for steps 3 and 4 below:
% waveform extraction (clustering) opts... (happen to be same as used above)
co = []; co.cmethod='k++'; co.K = 10; co.Kkeep = co.K; co.thresh=100; co.verb=1;
% fitting opts...
so = []; so.verb = 1; so.skip = 5; so.nlps = 10;
S = @(Y) spikesort_timeseries(Y,wftrue.d.samplefreq,co,[],so); % interface


% 3. Test sorting accuracy --------------------------------------------
[t,l,~,wf] = S(Y);   % sort (few sec)
o.max_matching_offset=10;
[permL2 P acc times]=times_labels_accuracy(ptrue.t,ptrue.l,t,l,o); % a few sec
% acc.p are the accuracy f_k values (not acc.f, careful!)
figure; imagesc(P); colormap(goodbw); axis equal tight; title('accuracy Q'); colorbar; xlabel('sorted k (permed)'); ylabel('truth k');
[~,invp] = sort(permL2); Wp = wf.W(:,:,invp);     % apply best perm to W
plot_spike_shapes(Wp,'sorted W (permed)');

spikespy({Y,ptrue.t,ptrue.l,'Y+truth'}, {Y,t,permL2(l),'Y+sorting(permed)'});

stop


% 4. Stability metric for sorting ----------------------------------
% general stability metric options (verb=3 allows paper outputs)
o = []; o.Nt = 60; o.max_matching_offset=10; o.verb = 3;
% stability opts specific to addition metric...
%o.meth = 'add'; o.ratescale = 0.25; o.num_runs = 20; % 20 runs takes 3-4 mins
o.meth = 'rev';
[fahat,fasam,infoa] = eval_stability_tseriesbased(S, Y, o);

%save data_valid/acc_vs_stab_ampl.4_45swap.mat

% 5. just an attempt to plot for now: we may not even need a fig -----------
sso=[]; sso.ylab='f_k';               % start stability vs accuracy fig...
show_stabilities(fahat,fasam,sso);
sso.blobcolor=[.9 .4 0]; sso.fig = gcf; % overlay accuracies in orange
show_stabilities(acc.p,[],sso);
sso=[]; sso.fig = gcf; show_stabilities(fahat,[],sso);  % re-overlay black blobs
title('orange = acc, black = stability')
set(gcf,'paperposition',[0 0 5 5]);
% you can see good correlation, but accuracy is 2-3x closer to 1 than stability.

%print -depsc2 acc_vs_stab_ampl.4_Poissonfiring.eps

% Leslie discussed stating corr coeff of the two...
figure; plot(acc.p, fahat, '+'); xlabel('acc'); ylabel('f_k hat');
title('scatter for ampl std dev = 0.4, K=10');
%print -depsc2 acc_vs_stab_ampl.4_Poissonfiring_scatter.eps

% ie acc.p and fahat
