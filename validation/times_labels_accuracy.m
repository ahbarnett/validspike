function [permL2 P acc times]=times_labels_accuracy(T1,L1,T2,L2,opts)
% TIMES_LABELS_ACCURACY - output comparison stats between two {t,l} event lists
%
% [permL2 C acc times] = times_labels_accuracy(T1,L1,T2,L2,opts)
%  gives a report card on accuracy of times T2 and labels L2 against the "truth"
%  dataset of times T1 and labels L1.
%  Gives a warning if number of problem times doesn't match number of off-diag
%  entries in confusion matrix.
%
% Inputs:
%  T1,L1 - list of times (real) and labels (in 1..K1) for firing event list 1
%  T2,L2 - list of times (real) and labels (in 1..K2) for firing event list 2
%  opts  - (optional) controls various parameters such as:
%          o.verb - verbosity 0,1,...
%          Other opts are passed to times_labels_confusion_matrix
% Outputs:
%  permL2 - best permutation of labels in L2 to match T1 and L1.
%  C  - extended (L2-permuted) confusion matrix, size (K1+1) by (K2+1)
%  acc - struct with per-label accuracy statistics:
%          acc.tI - type-I error rates, ie missed
%          acc.m - type-I error rates only counting unclassified as missed
%          acc.tII - type-II error rates, ie false positive
%          acc.f - type-II error rates only counting false positive with respect
%                  to being unclassified
%          acc.p - stability-type indicator
%  times - struct with times of problems:
%          times.tfals, times.tmiss, times.twrng - lists of times of
%          false-positive, missed, and wrong labelled events respectively
%          (useful when plotted; see SORT_TIMESERIES self-test)
%
% Notes: there is no provision for unclassified events (entries of L1 or L2 that
%  are not in 1...K1 or 1....K2 respectively).
%
% See also: TIMES_LABELS_CONFUSION_MATRIX

% Barnett 6/3/15, 6/12/15. New output interface 6/16/15

if nargin==0, test_times_labels_accuracy; return; end
if nargin<5, opts = []; end
if ~isfield(opts,'verb'), opts.verb=(nargout==0); end

fprintf('doing times_labels_confusion_matrix...'), t1=tic;
[P times.tfals times.tmiss times.twrng permL2] = times_labels_confusion_matrix(T1,L1,T2,L2,opts); % the meat
fprintf(' done in %.3g s\n',toc(t1));
Nf=numel(times.tfals); Nm=numel(times.tmiss); Nw=numel(times.twrng);
Ntot = Nf+Nm+Nw;
fprintf('t_l_compare: %d mistakes (%d false pos, %d missed, %d wrong label)\n',Ntot,Nf,Nm,Nw)

[P acc] = labels_similarity(P);

K = min(size(P))-1;  % smaller of K1 and K2
Pod = P; for i=1:K, Pod(i,i) = Pod(i,i) - P(i,i); end  % off-diag part of P
if Ntot~=sum(Pod(:))
  warning(sprintf('sum of off-diag P entries (%d) differs from total # mistakes from time lists (%d)!', sum(Pod(:)), Ntot));
end

%%%%%%%%%%%%%%%
function test_times_labels_accuracy
o = [];
T=[10 20 30]; L=[1 2 3]; S=[10 20]; M=[1 2];
'each should give one mistake, of different type:'
times_labels_accuracy(T,L,S,M,o);
times_labels_accuracy(S,M,T,L,o);
times_labels_accuracy(T,[2 2 3],T,L,o);
