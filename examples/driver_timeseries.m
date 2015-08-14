% example script for spike-sorting a time-series
% Barnett 6/10/15

clear
d = loaddemodata;                  % load some timeseries data into d.A array

d.A = freqfilter(d.A,d.samplefreq,300,[]);  % filter and noise-unmix channels
d = channelprewhiten(d,[]);

[t l p wf R] = spikesort_timeseries(d.A,d.samplefreq);  % do the alg

plot_spike_shapes(wf.W,'waveforms found'); drawnow;
if exist('spikespy'), spikespy({d.A,t,l,d.name},{R,'residual'}); end

% assess accuracy, since we have luxury of ground-truth...
g = load('EC_default_synth_groundtruth.mat');
times_labels_accuracy(g.p.t,g.p.l,t,l);  % good since synthetic data
