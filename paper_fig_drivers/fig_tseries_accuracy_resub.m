% Starting the stability metric vs accuracy comparison for JNM resubmission.
% Barnett 11/20/15
% For Jeremy to edit..!

% 1. Extract waveforms by doing spikesorting (20 sec laptop) ----------------
clear; d = loaddata('e');  % EJ 2005 elec359 timeseries from data_external/
d.A = freqfilter(d.A,d.samplefreq,300,[]);
% waveform extraction opts... (these are as in fig_tseries_meths.m)
co.cmethod='k++'; co.K = 10; co.Kkeep = 10; co.thresh=100; co.verb = 2;
% fitting opts...
so.skip = 5; so.nlps = 10;
[~,~,p,wf] = spikesort_timeseries(d.A,d.samplefreq,co,[],so);
wf.freqs = histc(p.l,1:co.Kkeep)/d.T;   % mean firing rates (Hz) of each label
noi = empiricalnoise(d);            % get noi.eta = noise std deviation.
% the point was to get a wf (wavefunc object) and noi (noise model); clear rest
clear p co so


% 2. Create fresh synthetic data (post-frequency-filtered) -----------------
N = 1e6;             % length in samples
rates = wf.freqs; % use the same rates as in data...
%rates(9) = 0;    % .. but kill noise cluster
tpad = 10;           % end padding in samples
so = [];  % synth opts...
so.truePois = 1;
so.ampl = 0.2;     % *** firing variation they wanted: relative ampl std dev
noi.Nt = N; noi.M = size(wf.W,1);   % basic params for noise model
[Y ptrue] = synth_Poissonspiketrain(wf,N,rates,noi,tpad,[],so);
fprintf('synth done, %d spikes\n',numel(ptrue.t));
% *** make a code edited from synth_Poissonspiketrain if you want to change the ptrue.t timing synthesis, eg to make correlated firings, etc!
% Or, fill your own ptrue.t .l .a and just call
%    Y = spikemod(wf, ptrue, N);
%    Y = Y + noisesample(noi);
% to run fwd model and add iid noise
% *** adding correlated noise would be better, by doing instead
%    Y = Y + freqfilter(noisesample(noi),wf.d.samplefreq,300,[]));
% We can discuss this.

spikespy({Y,ptrue.t,ptrue.l,'synth'});   % show what we made
show_crosscorr(ptrue.l,ptrue.t);  % no refractory holes, of course, since
                                  % purely Poisson spike trains for now


% Setup inline sorter function S for steps 3 and 4 below:
% waveform extraction (clustering) opts... (happen to be same as used above)
co = []; co.cmethod='k++'; co.K = 10; co.Kkeep = co.K; co.thresh=100; co.verb=1;
% fitting opts...
so = []; so.verb = 1; so.skip = 5; so.nlps = 10;
S = @(Y) spikesort_timeseries(Y,wf.d.samplefreq,co,[],so); % interface: t,l out


% 3. Test sorting accuracy --------------------------------------------
[t l] = S(Y);   % sort (few sec)
o.max_matching_offset=10;
[permL2 P acc times]=times_labels_accuracy(ptrue.t,ptrue.l,t,l,o); % a few sec
% acc.p are the accuracy f_k values


% 4. Stability metric for sorting ----------------------------------
% general stability metric options (verb=3 allows paper outputs)
o = []; o.Nt = 60; o.max_matching_offset=10; o.verb = 3;
% stability opts specific to addition metric...
o.meth = 'add'; o.ratescale = 0.25; o.num_runs = 20;
[fahat,fasam,infoa] = eval_stability_tseriesbased(S, Y, o);


% 5. just an attempt to plot for now: we may not even need a fig -----------
sso=[]; sso.ylab='f_k';               % start stability vs accuracy fig...
show_stabilities(fahat,fasam,sso);
sso.blobcolor=[.9 .4 0]; sso.fig = gcf; % overlay accuracies in orange
show_stabilities(acc.p,[],sso);
sso=[]; sso.fig = gcf; show_stabilities(fahat,[],sso);  % re-overlay black blobs
title('orange = acc, black = stability')
set(gcf,'paperposition',[0 0 5 5]);
% you can see good correlation, but accuracy is 2-3x closer to 1 than stability.

% Leslie discussed stating corr coeff of the two...
