% understanding accuracy and stability results
% Barnett 2/8/16
clear
%load data_valid/accstabsam_rev_nss50_N2400000_eta20_ampl.2.mat
%load data_valid/accstabsam_add_r10_nss10_N2400000_eta20_ampl.2.mat
%%%%%%%% try again w/. more info (verb=5 for 'add' metric):   2/9/16
load data_valid/accstabsam_add_wfap_N2400000_eta20_ampl0.mat

plot_spike_shapes(wftrue.W,'true W')
o  % check stability metric

i=1;   % which acc realization to look at
figure; imagesc(P{i}); colorbar; axis equal tight; title('accuracy')
disp('accuracy:'), P{i}

fasam{:}

i=9;   % which stab realization to look at (indep of above)
[~,iperm] = sort(sperm{i});
iperm = [iperm K+1];          % don't permute the no-spike index
if iscell(info{i}.Qs)           % 'add' metric
  r = 1;                       % which add run to look at - search for f_4 v neg!
  Qp = info{i}.Qs{r}(iperm,iperm);   % permute stab conf mat
  plot_spike_shapes(info{i}.wfap{r}.W(:,:,iperm(1:K)),'permed W, rth run')
else                            % 'rev' metric
  Qp = info{i}.Qs(iperm,iperm);   % permute stab conf mat
end
figure; imagesc(Qp); colorbar; axis equal tight; title('perm conf mat')
disp('permuted Q conf mat:'), Qp
plot_spike_shapes(Wp{i},'permed W from 1st stab run');




