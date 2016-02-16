% figures for time-series validation, EJ 2005 elec359 dataset
% Barnett 7/22/15

clear; d = loaddata('e');  % gets from data_external/
d.A = freqfilter(d.A,d.samplefreq,300,[]);

% Setup inline sorter function S:
% waveform extraction (clustering) opts...
co = []; co.cmethod='k++'; co.K = 10; co.Kkeep = co.K; co.thresh=100; co.verb=1;
so = []; so.verb = 1; so.skip = 5; so.nlps = 10; % fitting opts
S = @(Y) spikesort_timeseries(Y,d.samplefreq,co,[],so);  % interface: tj, kj out

% general stability metric options (verb=3 allows paper outputs)
o.Nt = 60; o.max_matching_offset=10; o.verb = 3;

[o.T o.L p wf] = S(d.A);  % cheat by using a single unpert run, grab its outputs

nk = histc(o.L,1:co.K); K = co.K;
sc = 350; plot_spike_shapes(wf.W,[],sc);
for k=1:K, text((k-1)*(size(wf.W,2)+4),-100,sprintf('%d',nk(k))); end
text(K*(size(wf.W,2)+4),-100,'$$n_k$$','interpreter','latex')
set(gcf,'paperposition',[0 0 5.5 4.5]);
print -depsc2 ~/spikesorting/validpaper/tseries_W.eps

show_crosscorr(o.L,o.T);  % physiological (refractory hole) check on unpert run
set(gcf,'paperposition',[0 0 5 5]);                                   
print -depsc2 ../spikesorting/validpaper/tseries_xcorr.eps

o.meth = 'rev';                                   % ------- noise rev
[fhat,fsam,info] = eval_stability_tseriesbased(S, d.A, o);

[Tnr Lnr wfnr] = S(info.Ynr);  % redo the rev SS run (assume same as stab...)
wfnr = pull_waveforms_from_tseries(info.Ynr,info.Tnr,info.Lnr,o.Nt,o); % rev W

[~,invp] = sort(info.permL2); Wnr = wfnr.W(:,:,invp); % reorder Wnr to match
nknr = histc(info.permL2(info.Lnr),1:K);  % reordered pops
plot_spike_shapes(Wnr,[],sc);  % permuted ones
%plot_spike_shapes(wfnr.W,[],sc); nknr = histc(info.Lnr,1:K);
for k=1:K, text((k-1)*(size(wfnr.W,2)+4),-100,sprintf('%d',nknr(k))); end
text(K*(size(wfnr.W,2)+4),-100,'$$n_k$$','interpreter','latex')
set(gcf,'paperposition',[0 0 5.5 4.5]);
print -depsc2 ~/spikesorting/validpaper/tseries_Wnr.eps

o.meth = 'add'; o.ratescale = 0.25; o.num_runs = 20;   % ----- addition
[fahat,fasam,infoa] = eval_stability_tseriesbased(S, d.A, o);

sso=[]; sso.ylab='f_k';                              % start stability fig...
show_stabilities(fahat,fasam,sso);
sso.blobcolor=[.9 .4 0]; sso.fig = gcf; % overlay orange rev results
show_stabilities(fhat,fsam,sso);
sso=[]; sso.fig = gcf; show_stabilities(fahat,[],sso);  % re-overlay black blobs
set(gcf,'paperposition',[0 0 5 5]);
print -depsc2 ../spikesorting/validpaper/tseries_fk.eps

%save data_valid/tseries_thresh120K9keep9.mat
%save data_valid/tseries_thresh100K10keep10.mat

if 0  % regen the f_k plot...
  load data_valid/tseries_thresh100K10keep10.mat
  sso=[]; sso.ylab='f_k';                              % start stability fig...
  show_stabilities(fahat,fasam,sso);
  sso.blobcolor=[.9 .4 0]; sso.fig = gcf; % overlay orange rev results
  show_stabilities(fhat,fsam,sso);
  sso=[]; sso.fig = gcf; show_stabilities(fahat,[],sso);  % re-overlay black blobs
  text(.7,.24,'noise-reversal','color',[.9 .4 0]);
  text(.7,.14,'spike addition','color',[0 0 0]);
  set(gcf,'paperposition',[0 0 5 5]);
  print -depsc2 ../spikesorting/validpaper/tseries_fk.eps
end

if 0 % regen the t-series plot...
  load data_valid/tseries_thresh100K10keep10.mat
  spikespy({d.A,o.T,o.L,'Y, unperturbed run'},{info.Ynr,info.Tnr,info.permL2(info.Lnr),'Ynr, noise-reversed run'});
  % get scale correct:
  min(d.A(2,1:40*20))  % -355 sets vert scale in plot
end

if 0 % regen the xcorr plot...
  load data_valid/tseries_thresh100K10keep10.mat
  show_crosscorr(o.L,o.T);  % physiological (refractory hole) check on unpert run
  set(gcf,'paperposition',[0 0 5 5]);                                   
  print -depsc2 ../spikesorting/validpaper/tseries_xcorr.eps
end

if 0 % report the number of "overlapping" (nearby) spikes...
  t = sort(o.T);
  fprintf('tot spikes/sec = %.3g\n',numel(t)/(d.dt*max(t)))
  fprintf('frac spikes within 0.5 ms of prev = %.3g\n',numel(find(diff(t)<=10))/numel(t))
  fprintf('frac spikes within 0.1 ms of prev = %.3g\n',numel(find(diff(t)<=2))/numel(t))
  fprintf('frac spikes within 0.5 ms of another = %.3g\n',numel(find(diff(t(1:end-1))<=10 | diff(t(2:end))<=10))/numel(t))
  fprintf('frac spikes within 0.1 ms of prev = %.3g\n',numel(find(diff(t(1:end-1))<=2 | diff(t(2:end))<=2))/numel(t))
end
