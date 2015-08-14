% figure for validation by simple rerunning: r-convergence
% Barnett 6/26/15
%
% todo: add r=... labels and (a) etc to figs

clear
load data_valid/clips_bb_short_th120_3ms_fac3.mat
fprintf('loaded\n')

if 0 % warmup
  o.cmethod='k++'; o.K = 8; o.num_trials = 10;  % ss alg opts (K,r)
  oo.meth='rerun'; oo.num_runs = 20; oo.verb = 2; oo.ylab='f_k^{rerun}'; % stability metric opts
  [fhat,fsam] = eval_stability_clipbased(@spikesort_clips, X, o, oo);
end

% fixed K, loop over various r
rs = [1 5 10 20 100];  % num_trials of k-means++.  Note 100 takes 8 sec per run
for i=1:numel(rs)
  o.cmethod='k++'; o.K = 8; o.num_trials = rs(i);  % ss alg opts (K,r)
  oo.meth='rerun'; oo.num_runs = 20; oo.verb = 2; oo.ylab='f_k^{rerun}'; % stability metric opts
  [fhat,fsam] = eval_stability_clipbased(@spikesort_clips, X, o, oo);
  title(sprintf('r=%d',rs(i)))
  set(gcf,'paperposition',[0 0 4 4]);
  print('-depsc2',sprintf('~/spikesorting/validpaper/clips_rerun_bb_K%d_r%d.eps',o.K,rs(i)))
end
