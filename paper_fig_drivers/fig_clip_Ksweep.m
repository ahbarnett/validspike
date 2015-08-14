% figures for clips-based validation on standard Buzsaki clips, varying K
% Barnett 7/4/15
clear; load data_valid/clips_bb_short_th120_3ms_fac3.mat; fprintf('loaded\n')

% SS alg to validate...
o=[]; o.verb=1; o.cmethod='k++'; o.K = 8; o.num_trials = 100; % Nfea is default

vo=[]; vo.num_runs = 20; vo.verb = 1;
%vo.meth = 'blur'; vo.gamma = 0.5;
vo.meth = 'rev';
Ks = 2:15;
n = numel(Ks);
for i=1:n, o.K = Ks(i);    % sweep K
  [fhat{i},fsam{i},info{i}] = eval_stability_clipbased(@spikesort_clips, X, o, vo);
end

clear X; %save data_valid/clip_rev_Ksweep.mat
%load data_valid/clip_rev_Ksweep.mat    % .. to reoutput figs
 
Ks = 2:12; n = numel(Ks);  % if want to truncate plot
F = nan(n,max(Ks));
for i=1:n, for j=1:Ks(i), F(i,j) = fhat{i}(j); end, end % upper-tri array
figure; imagesc(1:max(Ks),Ks,F); colormap(goodbw); colorbar
axis xy; xlabel('k'); ylabel('K');
%set(gcf,'paperposition',[0 0 4 5]); % vertical format
set(gcf,'paperposition',[0 0 6 2.5]);
if strcmp(vo.meth,'rev')
  caxis([.75 1]);
  text(6,3,'f_k^{rev}, various K');
  print('-depsc2',sprintf('~/spikesorting/validpaper/clips_rev_Ksweep_r%d.eps',o.num_trials))
else                 % blur
  caxis([0 1]);
  text(6,3,'mean f_k^{blur}, \gamma=0.5, various K');
  print('-depsc2',sprintf('~/spikesorting/validpaper/clips_blur_gam05_Ksweep_r%d.eps',o.num_trials))
end
