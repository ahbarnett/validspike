% figures for time-series noise-reversal validation, EJ 2005 elec359 dataset
% Barnett 7/17/15

clear; d = loaddata('e');  % gets from data_external/
d.A = freqfilter(d.A,d.samplefreq,300,[]);

% waveform extraction (clustering) opts:
co = []; co.cmethod='k++'; co.K = 10; co.Kkeep = 10; co.thresh=90; co.verb = 2;
% timeseries fitting opts:
so = []; so.verb = 1; so.skip = 5; so.nlps = 10;
S = @(Y) spikesort_timeseries(Y,d.samplefreq,co,[],so);  % interface: tj, kj out

% stability options (verb=3 for paper)
o.meth = 'rev'; o.Nt = 60; o.max_matching_offset=10; o.verb = 3;
[fhat,fsam,info] = eval_stability_tseriesbased(S, d.A, o); info

show_crosscorr(info.L,info.T);  % physiological (refractory hole) check, run 1
stop
sso=[]; sso.ylab='f_k'; sso.blobcolor=[.9 .4 0]; show_stabilities(fhat,fsam,sso);
hgsave data_valid/tseries_rev_K8keep8.fig
%d.A = []; save data_valid/tseries_rev_K9keep9.mat

% figs: first grab W by hand...
title ''
set(gcf,'paperposition',[0 0 4.5 4.5]);
print -depsc2 ~/spikesorting/validpaper/tseries_K9_W.eps

show_crosscorr(info.L,info.T);  % physiological (refractory hole) check, run 1
set(gcf,'paperposition',[0 0 5 5]);                                   
print -depsc2 ../spikesorting/validpaper/tseries_rev_K9keep9_xcorr.eps


o.pops = info.pops; show_stabilities(fhat,fsam,o);
set(gcf,'paperposition',[0 0 5 5]);
print -depsc2 ../spikesorting/validpaper/tseries_rev_K9keep9_fk.eps
