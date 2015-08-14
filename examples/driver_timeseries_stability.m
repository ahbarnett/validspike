% example script for measuring stability of spike-sorting on a time-series.
% Based on driver_timeseries.m
% Barnett 8/14/15

clear
d = loaddemodata;                  % load some timeseries data into d.A array

%d.A = freqfilter(d.A,d.samplefreq,300,[]);  % filter and noise-unmix channels
%d = channelprewhiten(d,[]);

% clustering opts and sorting opts...
co.cmethod='k++'; co.K = 7; co.thresh=120; co.verb=1;
so.verb = 1; so.skip = 5; so.nlps = 10;
S = @(Y) spikesort_timeseries(Y,d.samplefreq,co,[],so);  % sorter interface

vo.meth = 'add'; vo.ratescale = 0.25; vo.num_runs = 5;   % validation opts
vo.verb = 2;
[fbar,fsam,info] = eval_stability_tseriesbased(S, d.A, vo);  % do it (<1 min)

% with current demo data, this shows exceptionally high stabilities
