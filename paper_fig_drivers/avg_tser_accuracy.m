% explore average accuracy for synth data w/ gnd truth.
% See fig_tseries_accuracy_resub.m  for making the gndtruth data.
% Barnett 1/21/16

clear
load data_valid/synth_accuracy_gndtruth_45swap_vs_paper.mat
wf.d.samplefreq = 2e4;

% waveform extraction (clustering) opts... (happen to be same as used above)
co = []; co.cmethod='k++'; co.K = 10; co.Kkeep = co.K; co.thresh=100; co.verb=1;
% fitting opts...
so = []; so.verb = 1; so.skip = 5; so.nlps = 10;
S = @(Y) spikesort_timeseries(Y,wf.d.samplefreq,co,[],so); % interface: t,l out

nr = 10;
% 3. Test sorting accuracy --------------------------------------------
for r=1:nr
  fprintf('run %d:\n',r)
  [t{r} l{r}] = S(Y);   % sort (few sec)
  o.max_matching_offset=10;
  [permL2 P{r} acc{r} times{r}]=times_labels_accuracy(ptrue.t,ptrue.l,t{r},l{r},o); % a few sec
  % acc.p are the accuracy f_k values (not acc.f, careful!)
  figure; imagesc(P{r}); colormap(goodbw); title(sprintf('run %d',r)); colorbar
end
% compute avg acc:
accbar = 0*acc{1}.p; for r=1:nr, accbar = accbar+acc{r}.p; end
accbar = accbar/nr;


figure
for r=1:nr
  subplot(3,4,r); imagesc(P{r}); colormap(goodbw);title(sprintf('run %d',r)); colorbar; drawnow
end
set(gcf,'name','for each run, truth is down and sorting across');
