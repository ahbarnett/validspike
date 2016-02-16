% plot the data made by tseries_avg_accruns_stabruns.m for resubmission
% Barnett 2/5/16. Big rects removed, arrows in matlab, 2/10/16. Run from ../
% Uses: errorbar_h.m

clear; figdir='../spikesorting/validpaper/';
qs = [.25 .75];    % quantiles defining boxes

for meth={'rev','add'}
  m=meth{1};  % un-celled string
  figure
  if strcmp(m,'rev'), r=load('data_valid/accstabsam_rev_nss50_N2400000_eta20_ampl.2.mat');
    h=patch([1 1 .98 .7],[.2 .95 .95 .2], [1 .92 .92]);  % danger region
    set(h,'edgecolor','none')
  else, r=load('data_valid/accstabsam_add_r10_nss10_N2400000_eta20_ampl.2.mat'); 
    % load 2nd stab expt too:
    r2 = load('data_valid/accstabsam_add_wfap_N2400000_eta20_ampl02.mat');
    r.fahat = r2.fahat; % use 2nd expt
    h=patch([1 1 .98 .7],[-2.4 .8 .8 -2.4], [1 .92 .92]);
    set(h,'edgecolor','none');
  end
  K = r.K;
  fhats = vertcat(r.fahat{:});  % ns-by-K stability values
  fhatm = mean(fhats,1); accm = mean(r.accsam,1);   % means
  fhatv = var(fhats,1); accv = var(r.accsam,1);   % estimated variances
  fhate = sqrt(fhatv/size(fhats,1)); acce = sqrt(accv/size(r.accsam,1)); % std dev in est of means
  hold on
  for k=1:K
    x = accm(k)+[-1 1]*acce(k); y = fhatm(k)+[-1 1]*fhate(k);
    h = patch(x([1 2 2 1]), y([1 1 2 2]), [.8 .8 1]); % std err of mu estimators
    set(h,'edgecolor','none')
  end
  qx = quantile(r.accsam,qs); qy = quantile(fhats,qs);  % quantiles down cols
  h=errorbar_x(accm,fhatm,accm-qx(1,:),qx(2,:)-accm,'k.');
  set(h(1),'color',[.6 .6 .6]);
  h=errorbar(accm,fhatm,fhatm-qy(1,:),qy(2,:)-fhatm,'k.');
  set(h(1),'color',[.6 .6 .6]);
  plot(accm,fhatm,'k.','markersize',15);
  tfm = fhatm; tam = accm;               % text locations, now hack them...
  if strcmp(m,'rev'), tam(1:5) = 1.02; tfm(1:5) = 0.8-.06*(0:4);
    %drawArrow([1.03 .85],[1 0.97]);
  else, tam([1 2 3 5]) = 0.5+0.05*(0:3); tfm([1 2 3 5]) = 1.0;
    %drawArrow([.85 1],[.97 1]);
  end          % end hack
  for k=1:K, text(tam(k)+.01,tfm(k),sprintf('%d',k),'fontsize',14); end
  xlabel('accuracy $\overline{a}_k$ (mean over realizations)','interpreter','latex'); ylabel('$\overline{f}_k$ (mean over realizations)','interpreter','latex');
  if strcmp(m,'rev'), title('synthetic time-series: noise-reversal')
  else, title('synthetic time-series: spike addition'); end
  box off       % who set this on?
  drawnow
  set(gcf,'paperposition',[0 0 4 3]);
  if strcmp(m,'rev')
    axis([0 1 .2 1]);
    print('-depsc2',[figdir 'accstab_rev_ampl02.eps'])
  else
    axis([0 1 -2.4 1]);
    print('-depsc2',[figdir 'accstab_add_ampl02.eps'])
  end
end
k