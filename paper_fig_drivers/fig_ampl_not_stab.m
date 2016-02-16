% show waveform amplitude not simply stability.
% Barnett 2/5/16. run from ../

clear; figdir='../spikesorting/validpaper/';
qs = [.25 .75];    % quantiles defining boxes
lw = 5;  % quantile line width
figure;     % start a joint plot
for meth={'rev','add'}, m=meth{1};  % un-celled string
  if strcmp(m,'rev'), load data_valid/accstabsam_rev_nss50_N2400000_eta20_ampl.2.mat
    col = [.9 .4 0];  % orange
  else, load data_valid/accstabsam_add_r10_nss10_N2400000_eta20_ampl.2.mat
    r2 = load('data_valid/accstabsam_add_wfap_N2400000_eta20_ampl02.mat');
    fahat = r2.fahat; % use 2nd expt
    col = [0 0 0];   % black
  end
  fhats = vertcat(fahat{:});  % ns-by-K stability values
  fhatm = mean(fhats,1); accm = mean(accsam,1);   % means
  fhatv = var(fhats,1); accv = var(accsam,1);   % estimated variances
  fhate = sqrt(fhatv/size(fhats,1)); acce = sqrt(accv/size(accsam,1)); % std dev in est of means
  ampl = squeeze(max(max(abs(wftrue.W),[],1),[],2));    % peak abs ampl
  %ampl = sqrt(squeeze(sum(sum(wftrue.W.^2,1),2)));  % l2 norm
  plot(ampl, fhatm, 'k.','markersize',20,'color',col); hold on;
  tfm = fhatm; ta = ampl;               % text locations, if need to hack
  for k=1:K
    text(ta(k)+7,tfm(k),ta(k),sprintf('%d',k),'fontsize',14,'color',col);
  end
  ylabel('$\overline{f}_k$ (mean over realizations)','interpreter','latex'); xlabel('peak amplitude $\| W^{(k)} \|_\infty$ ($\mu$V)','interpreter','latex');
  %if strcmp(m,'rev'), title('rev'); %axis([0 200 0.5 1]);
  %else, title('add'); %axis([0 200 -.6 1]); end
  legend({'noise-reversal','spike addition'},'location','east')
  drawnow
end
title('stability metrics vs waveform peak amplitude')
%noise-reversal (orange) and addition (black)
axis([0 350 -2.4 1])
set(gcf,'paperposition',[0 0 4 3])
print('-depsc2',[figdir 'amplstab_ampl02.eps']);
