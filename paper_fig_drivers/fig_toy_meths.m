% all figures for clips-based validation on "toy" clips dataset, K=4, Nfea=3
% Barnett 7/1/15

clear; load data_valid/clips_bb_short_th120_3ms_fac3.mat; disp('clips loaded')

% make toy subset from types 1,2,3 only:
o=[]; o.verb=1; o.cmethod='k++'; o.K = 8; o.num_trials = 100; % ss alg opts (K,r)
[l W] = spikesort_clips(X,o);
i = find(l==1 | l==2 | l==3); Xtoy = X(:,:,i);
fprintf('N=%d\n',numel(i))

% switch to toy sorting params...
o.num_fea = 2;   % makes more stable
o.K = 4; o.verb = 1; [L Wk z] = spikesort_clips(Xtoy,o); % check SS works

if 0 % make figs which don't depend on validation method...
  plot_labeled_pts(z,L);
  view(2); v = axis; v(2) = 100; axis(v);
  title(sprintf('toy dataset N_{fea}=2, K=%d',o.K))
  h = findobj(gcf,'Type','axes','Tag','legend'); set(h,'location','northeast');
  set(gcf,'paperposition',[0 0 3 3]);  
  print -depsc2 ~/spikesorting/validpaper/toy_K4_z1z2.eps
  
  plot_spike_shapes(Wk);         % should show a split waveform
  K = o.K; nk = histc(L,1:K);    % pops
  for k=1:K, text((k-1)*(size(Wk,2)+4),-100,sprintf('%d',nk(k))); end
  text(K*(size(Wk,2)+4),-100,'$$n_k$$','interpreter','latex')
  title(sprintf('toy dataset, waveforms, K=%d',o.K))
  set(gcf,'paperposition',[0 0 2 4]);
  print -depsc2 ~/spikesorting/validpaper/toy_K4_W.eps
end
  
if 0 % Validation figs:
  vo.num_runs = 20; vo.verb = 1;  % validation global opts
  meths = {'cv3','blur','rev'};   % methods to try
  ylabs = {'CV','blur','rev'};    % corresp y-labels (only used for Nfac=1)
  noi = setup_noisemodel(d,size(X,2),30,0.0005); % fit from make_clips_file_bb
  
  o.verb = 1; % don't make z fig every SS alg run

  for Nfac = [1 100] % ---------- sweep size: normal & huge data set size
                     % note: Nfac = 100 is slow, of order 10 mins (for Nfea=3)
    Y = repmat(Xtoy,[1 1 Nfac]); N = size(Y,3);  % duplicate clips (2 GB RAM)
    if Nfac>1  % add noise
      for j=1:N, Y(:,:,j) = Y(:,:,j) + noisesample(noi); end % is slow (1 min)
    end
    for m=[3 1 2]              % ..... loop methods
      vo.meth = meths{m};
      [fhat,fsam,info] = eval_stability_clipbased(@spikesort_clips, Y, o, vo);
      so.ylab=sprintf('f_k^{%s}',ylabs{m}); show_stabilities(fhat,fsam,so);
      set(gcf,'paperposition',[0 0 2.2 4]);    
      if Nfac>1, ylabel ''; set(gcf,'paperposition',[0 0 1.9 4]);
      end    % kill the ylabel for more horiz space
      title(sprintf('toy: r=%d N=%d',o.num_trials,N)); drawnow
      print('-depsc2',sprintf('~/spikesorting/validpaper/toy_%s_K%d_r%d_N%d.eps',vo.meth,o.K,o.num_trials,N))
    end
  end
end

