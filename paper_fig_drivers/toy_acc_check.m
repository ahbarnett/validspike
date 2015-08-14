% explore accuracy for known toy dataset. Barnett 7/3/15

clear; load data_valid/clips_bb_short_th120_3ms_fac3.mat; disp('clips loaded')

% make toy subset from types 1,2,3 only:
o=[]; o.verb=1; o.cmethod='k++'; o.K = 8; o.num_trials = 100; % ss alg opts (K,r)
[l W] = spikesort_clips(X,o);
i = find(l==1 | l==2 | l==3); Xtoy = X(:,:,i); Ltoy = l(i);
fprintf('N=%d\n',numel(i))

% switch to toy sorting params...
o.num_fea = 2;   % makes more stable
%o.K = 3; o.verb = 1; [L3 Wk3 z3] = spikesort_clips(Xtoy,o); % same as Ltoy
o.K = 4; o.verb = 1; [L Wk z] = spikesort_clips(Xtoy,o);
labels_accuracy(Ltoy,L);