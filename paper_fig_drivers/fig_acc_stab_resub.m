% plot the data made by tseries_avg_accruns_stabruns.m for resubmission
% Barnett 2/5/16. Run from ../

clear; figdir='../spikesorting/validpaper/';
qs = [.25 .75];    % quantiles defining boxes
lw = 5;  % quantile line width
for meth={'rev','add'}
  m=meth{1};  % un-celled string
  if strcmp(m,'rev'), r=load('data_valid/accstabsam_rev_nss50_N2400000_eta20_ampl.2.mat');
  else, r=load('data_valid/accstabsam_add_r10_nss10_N2400000_eta20_ampl.2.mat'); 
    % load 2nd stab expt too:
    r2 = load('data_valid/accstabsam_add_wfap_N2400000_eta20_ampl02.mat');
    r.fahat = r2.fahat; % use 2nd expt
  end
  K = r.K;
  figure;
  fhats = vertcat(r.fahat{:});  % ns-by-K stability values
  fhatm = mean(fhats,1); accm = mean(r.accsam,1);   % means
  fhatv = var(fhats,1); accv = var(r.accsam,1);   % estimated variances
  fhate = sqrt(fhatv/size(fhats,1)); acce = sqrt(accv/size(r.accsam,1)); % std dev in est of means
  for k=1:K
    x = quantile(r.accsam(:,k),qs);
    y = quantile(fhats(:,k),qs);
    % transparency ruins the plotting order of blobs over lines: !
    %h=patch(x([1 2 2 1]), y([1 1 2 2]), [.8 .8 1],'facealpha',.3); % quantiles
    h=patch(x([1 2 2 1]), y([1 1 2 2]), [.8 .8 1]); % quantiles
    hold on; plot(accm(k)+[-1 1]*acce(k),fhatm(k)*[1 1],'-','color',[.5 .5 1],'linewidth',lw); % std errors in means
    plot(accm(k)*[1 1],fhatm(k)+[-1 1]*fhate(k),'-','color',[.5 .5 1],'linewidth',lw);
  end
  hold on; plot(accm,fhatm,'k.','markersize',15);  % show all means
  tfm = fhatm; tam = accm;               % text locations, now hack them...
  if strcmp(m,'rev'), tam(1:5) = 1.02; tfm(1:5) = 0.8-.06*(0:4);
    %arrow([1.03 1],[.85 0.97]);
  else, tam([1 2 3 5]) = 0.5+0.05*(0:3); tfm([1 2 3 5]) = 1.0;
    %tam(1:5) = 0.7+0.04*(0:4); tfm(1:5) = 0.55;  % for 2nd stab expt
    %arrow([.85 1],[.85 0.97]);
  end          % end hack
  for k=1:K, text(tam(k)+.01,tfm(k),sprintf('%d',k),'fontsize',14); end
  axis([0 1 0 1]); xlabel('accuracy $\overline{a}_k$ (mean over realizations)','interpreter','latex'); ylabel('$\overline{f}_k$ (mean over realizations)','interpreter','latex');
  if strcmp(m,'rev'), title('synthetic time-series: noise-reversal')
  else, title('synthetic time-series: spike addition'); end
  drawnow
  set(gcf,'paperposition',[0 0 4 3]);
  if strcmp(m,'rev'), axis([0 1 .2 1]);
  print('-depsc2',[figdir 'accstab_rev_ampl02.eps'])
  else, axis([0 1 -2.4 1]);
  print('-depsc2',[figdir 'accstab_add_ampl02.eps'])
  end
end

% now add an arrow to each fig in xfig, output as above _lab.eps
