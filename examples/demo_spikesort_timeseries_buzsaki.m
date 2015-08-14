% demo script for spike-sorting a time-series from Buzsaki data
% Barnett 6/11/15

clear
%d = loaddata('B'); fprintf('loaded\n')   % long-ish (8 min), 10 channel, filt
%d = loaddata('b'); fprintf('loaded\n')   % 50 sec, 10 channel
d = loaddata('bb'); fprintf('loaded\n')   % 2.5 min, 10 channel
%d = loaddata('b1'); fprintf('loaded\n')   % 2.5 min, 10 channel, gp1
d.A = freqfilter(d.A,d.samplefreq,300,[]);  % filter

d = channelprewhiten(d,[]);  % unmix

% sort in one go from just the timeseries...
co.thresh = 110; co.verb = 2; % classifying opts
so.skip=6; so.gamma=12;   % some sorting opts
[t l p wf R] = spikesort_timeseries(d.A,d.samplefreq,co,[],so);

plot_spike_shapes(wf.W,'waveforms found'); drawnow;
if exist('spikespy'), spikespy({d.A,t,l,d.name},{R,'residual'}); end
o.dtau = 20; o.taumax = 500; show_crosscorr(l,t,[],o);

% todo: show saving waveforms, reloading and sorting from them...
