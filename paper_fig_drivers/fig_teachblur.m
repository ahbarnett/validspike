% teaching (toy and 1d data) for self-blur and noise-rev.
% Barnett 7/2/15

clear; load data_valid/clips_bb_short_th120_3ms_fac3.mat; disp('clips loaded')

% make toy subset from types 1,2,3 only:  (same as fig_toy_meths)
o=[]; o.verb=1; o.cmethod='k++'; o.K = 8; o.num_trials = 100; % ss alg opts (K,r)
[l W] = spikesort_clips(X,o);
i = find(l==1 | l==2 | l==3); Xtoy = X(:,:,i);
fprintf('N=%d\n',numel(i))

% switch to toy sorting params...
o.num_fea = 2;   % makes more stable, viewable: makes Voronoi make sense.
o.K = 4; o.verb = 1; [L Wk z] = spikesort_clips(Xtoy,o); % check SS works

% ---- blur/noise teaching figures (needs Xtoy, L, Wk) ------------------------
K = o.K;
v = [-3000 300 -2200 1000];  % a consistent view on feature space
figure;
%subplot(2,3,1);   % was annoying to get spacing correct
axes('position',[0.01 0 .25 1]);
plot_labeled_pts(z(1:2,:),L,struct('fig',1,'nodotted',1)); % is 2d now anyway
view(2); axis(v); legend off
zcen = nan(size(z,1),K); for k=1:K, zcen(:,k) = mean(z(:,L==k),2); end
hold on; voronoi(zcen(1,:),zcen(2,:)); % nice!
h=scatter(zcen(1,:),zcen(2,:), 100, [1 0 0 1]'*[1 1 1], '+');
set(h,'linewidth',2);
box on; set(gca,'xtick',[],'ytick',[]); xlabel ''; ylabel ''
title(sprintf('(a) original toy data, K=%d',K))
Xp = Xtoy; o.gamma = 1.0;  % do a self-blurring (from eval_stability_clipbased)
for k=1:K
  j = find(L==k);
  d = bsxfun(@minus, Xtoy(:,:,j), Wk(:,:,k));     % differences from centroids
  Xp(:,:,j) = Xp(:,:,j) + o.gamma*d(:,:,randperm(numel(j)));    % do blur
end
[Lp Wp zp] = spikesort_clips(Xp,o);
%subplot(2,3,2);
axes('position',[.35 0 .25 1]);
zcenp = nan(size(zp,1),K); for k=1:K, zcenp(:,k) = mean(zp(:,Lp==k),2); end
voronoi(zcenp(1,:),zcenp(2,:)); hold on;  % new centroids
plot_labeled_pts(zp(1:2,:),L,struct('fig',1,'nodotted',1)); view(2); axis(v);
legend off
box on; set(gca,'xtick',[],'ytick',[]); xlabel ''; ylabel ''
title(sprintf('(b) self-blurred, \\gamma=%g',o.gamma))
Xp = Xtoy; % do a noise-reversal (from eval_stability_clipbased)
for k=1:K
  j = find(L==k);
  Xp(:,:,j) = bsxfun(@minus, 2*Wk(:,:,k), Xtoy(:,:,j));
end
[Lp Wp zp] = spikesort_clips(Xp,o);
%subplot(2,3,3);
axes('position',[.69 0 .25 1]);
zcenp = nan(size(zp,1),K); for k=1:K, zcenp(:,k) = mean(zp(:,Lp==k),2); end
voronoi(zcenp(1,:),zcenp(2,:)); hold on;  % new centroids
plot_labeled_pts(zp(1:2,:),L,struct('fig',1,'nodotted',1)); view(2); axis(v);
legend off
%h = findobj(gcf,'Type','axes','Tag','legend'); set(h,'location','northeast');
box on; set(gca,'xtick',[],'ytick',[]); xlabel ''; ylabel ''
title('(c) noise reversed')

if 0, set(gcf,'paperposition',[0 0 10 5]);  % sensitive
  print -depsc2 -painters ~/spikesorting/validpaper/teachblur.eps
end
