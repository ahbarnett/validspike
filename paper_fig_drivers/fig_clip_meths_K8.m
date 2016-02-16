% figures for clips-based validation on standard Buzsaki clips at K=8.
% Used in f:cv, f:blur
% Barnett 7/1/15

clear; load data_valid/clips_bb_short_th120_3ms_fac3.mat; fprintf('loaded\n')

% SS alg to validate...
o=[]; o.verb=1; o.cmethod='k++'; o.K = 8; o.num_trials = 100; % Nfea is default
[l W] = spikesort_clips(X,o); pops = histc(l,1:o.K);  % run it to get pops

if 0  %  ----------- separate figures
meths={'cv3','blur','rev'};
ylabs = {'CV','blur','rev'};    % corresp y-labels (only used for Nfac=1)
vo=[]; vo.num_runs = 20; vo.verb = 1;  % validation global opts

for m=1:3          % ...... methods
  vo.meth = meths{m};
  [fhat,fsam,info] = eval_stability_clipbased(@spikesort_clips, X, o, vo);
  so = []; if strcmp(vo.meth,'cv3'), so.pops = pops; end   % show pops at bottom
  so.ylab=sprintf('f_k^{%s}',ylabs{m}); show_stabilities(fhat,fsam,so);
  title(sprintf('r=%d',o.num_trials)); drawnow
  set(gcf,'paperposition',[0 0 4 4]);
  print('-depsc2',sprintf('~/spikesorting/validpaper/clips_%s_K%d_r%d.eps',vo.meth,o.K,o.num_trials))
end                % ......
end

if 1 % ----------or, make joint figure for blur & rev (7/8/15):
vo.meth = 'rev';
[fhatr,fsamr,infor] = eval_stability_clipbased(@spikesort_clips, X, o, vo);
vo.meth = 'blur';
[fhat,fsam,info] = eval_stability_clipbased(@spikesort_clips, X, o, vo);
so=[]; so.ylab='f_k'; so.blobcolor = [.9 .4 0]; show_stabilities(fhatr,fsamr,so);
so = []; so.fig = gcf; show_stabilities(fhat,fsam,so);
title(sprintf('r=%d',o.num_trials));
text(.7,.24,'noise-reversal','color',[.9 .4 0]);
text(.7,.14,'self-blurring','color',[0 0 0]);
drawnow
set(gcf,'paperposition',[0 0 4 4]);
print('-depsc2',sprintf('~/spikesorting/validpaper/clips_blur+rev_K%d_r%d.eps',o.K,o.num_trials))
end

if 0 % ---------- K=10 waveforms and confusion (7/8/15):
o=[]; o.verb=1; o.cmethod='k++'; o.K = 10; o.num_trials = 100; % Nfea is default
[l W] = spikesort_clips(X,o); pops = histc(l,1:o.K);  % run it to get pops
plot_spike_shapes(W,[],400); K = o.K;
nk = histc(l,1:K);
for k=1:K, text((k-1)*(size(W,2)+4),-100,sprintf('%d',nk(k))); end
text(K*(size(W,2)+4),-100,'$$n_k$$','interpreter','latex')
set(gcf,'paperposition',[0 0 5 2.5]);
print -depsc2 ~/spikesorting/validpaper/clips_K10_W.eps

vo=[]; vo.num_runs = 5;   % only 5 runs, but 20 is the same
vo.verb = 1; vo.meth = 'blur';  vo.gamma = 0.5;
[fhat,fsam,info] = eval_stability_clipbased(@spikesort_clips, X, o, vo);
Q = info.Qs{1}; nr = numel(info.Qs);
for i=2:nr, Q = Q+ info.Qs{i}; end, Q=Q/nr; % mean Q
figure; %imagesc(Q(1:K,1:K));
imagesc(diag(1./nk)*Q(1:K,1:K)); % relative (scaled) confusion Q_{kl}/n_k matrix
colormap(goodbw); colorbar
text(4,1,'self-blurring (\gamma=0.5)')
text(4,1.8,'scaled best confusion matrix')
text(4,2.5,'$$\hat{Q}_{k\tilde{k}}/n_k$$','interpreter','latex');
ylabel('$$k$$','interpreter','latex'); xlabel('$$\tilde{k}$$','interpreter','latex');
set(gca,'xtick',1:K,'ytick',1:K);
set(gcf,'paperposition',[0 0 5 3]);
print -depsc2 ~/spikesorting/validpaper/clips_K10_Q.eps
end
