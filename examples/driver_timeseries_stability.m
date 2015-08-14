% example script for measuring stability of spike-sorting on a time-series.
% Based on driver_timeseries.m
% Barnett 7/14/15

% *** todo

clear
d = loaddemodata;                  % load some timeseries data into d.A array
d.A = freqfilter(d.A,d.samplefreq,300,[]);  % filter and noise-unmix channels
d = channelprewhiten(d,[]);

% need to pick std interface to SS
%***
ss = @(Y,fs) spikesort_timeseries(Y,fs);

%[t l p wf R] = spikesort_timeseries(d.A,d.samplefreq);  % do the alg

% ***  this is all inline funcs... better to redo clips stab in same way??
[fbar,fsam] = eval_stability_tseriesbased(@(Y) ss(Y,d.samplefreq,copts,smeth,copts), X, o, oo);
% todo

