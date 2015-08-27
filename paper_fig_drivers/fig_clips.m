% validation paper f:clips
% Barnett 6/24/15. Legend blobs bigger 8/26/15

clear
% ------------ buzsaki:
load data_valid/clips_bb_short_th120_3ms_fac3.mat

o = []; o.cmethod = 'k++'; o.K = 8; % choose clustering method (this one needs K)
o.verb = 1; [l W z cinfo] = spikesort_clips(X,o);   % 8 sec
plot_labeled_pts(z,l);
set(gcf,'paperposition',[0 0 6 4]);
%h=get(gca,'children'); set(h,'markersize',20)
h=findobj(get(findobj(gcf,'tag','legend'),'children'),'type','line');
set(h,'markersize',20);   % make legend blobs only bigger
% now move the legend by hand before outputting...
print -depsc2 ~/spikesorting/validpaper/clipsz.eps

sc = 400;
overview_sorted_clips(X,l,'',[],struct('equal',1,'vs',sc),[],struct('nums',0)); % show sim # clips of each k, amplified a bit
set(gcf,'paperposition',[0 0 16 4]);
print -depsc2 ~/spikesorting/validpaper/clips.eps

plot_spike_shapes(W,[],sc); K = o.K;
nk = histc(l,1:K);
for k=1:K, text((k-1)*(size(Wk,2)+4),-100,sprintf('%d',nk(k))); end
text(K*(size(Wk,2)+4),-100,'$$n_k$$','interpreter','latex')
set(gcf,'paperposition',[0 0 4.5 4.5]);
print -depsc2 ~/spikesorting/validpaper/clipsW.eps

% following from fig_rerun_rconv w/ r = 10 only:
rs = [10]; % num_trials of k-means++
for i=1:numel(rs)
  o.cmethod='k++'; o.K = 8; o.num_trials = rs(i);  % ss alg opts (K,r)
  oo.meth='rerun'; oo.num_runs = 40; oo.verb = 2; oo.ylab='f_k^{rerun}'; % stability metric opts
  [fhat,fsam] = eval_stability_clipbased(@spikesort_clips, X, o, oo);
  set(gcf,'paperposition',[0 0 4 4]);
  print('-depsc2',sprintf('~/spikesorting/validpaper/clips_rerun_bb_K%d_r%d.eps',o.K,rs(i)))
end
