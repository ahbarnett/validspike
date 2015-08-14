function [permL2 Qe acc] = labels_accuracy(L1,L2,o)
% LABELS_ACCURACY  summarize best-permuted accuracy info between two label lists
%
% labels_accuracy(L1,L2) finds best matching between two same-length lists
%  of labels L1, L2, with entries in 1..K1 and 1..K2 respectively, then reports
%  error metrics. L1 is assumed to the be ground-truth, and L2 the set whose
%  accuracy is being assessed. Zeros are treated as unlabeled/unclassified.
%
% [permL2 Qe acc] = labels_accuracy(L1,L2,opts) returns the permutation of L2
%  that was best, the (L2-permuted) extended confusion matrix size (K1+1)*(K2+1)
%  where the last row and column count unclassified events in L1 and L2
%  respectively, and finally an accuracy struct containing fields:
%          acc.tI - type-I error rates, ie missed
%          acc.m - type-I error rates only counting unclassified as missed
%          acc.tII - type-II error rates, ie false positive
%          acc.f - type-II error rates only counting false positive with respect
%                  to being unclassified
%          acc.p - JFM accuracy indicators (row vector), is f_k in the paper.
%
% opts controls options such as:
%    opts.verb : verbosity 0,1,... (default 1 unless output args)
%
% See also: LABELS_SIMILARITY which defines the metrics in acc thoroughly.

% Barnett 6/11/15. Include wrong 6/12/15. Output interface 6/16/15
% replaced getbestshuffling with bestcolpermconfmat 8/14/15

if nargin==0, test_labels_accuracy; return; end
if nargin<3, o=[]; end
if ~isfield(o,'verb'), o.verb = (nargout==0); end

Q = confusion_matrix(L1,L2);       % not extended
[permL2] = bestcolpermconfmat(Q);
K2 = max(L2);
ii = (L2>=1 & L2<=K2);                % classified indices in L2
L2(ii) = permL2(L2(ii));              % reshuffle L2's classified labels
if o.verb
  disp('Labels_accuracy. best permutation of labels L2: ')
  fprintf('%5d ',permL2), fprintf('\n')
end
[Qe acc] = labels_similarity(L1,L2,o);  % now show extended (a little wasteful)
%%%%%%%%%

function test_labels_accuracy
L1 = [1 1 2 2 3 3];
% missing, false pos:
L2 = [1 2 2 3 3 3]; [pL2,Qe,acc]=labels_accuracy(L1,L2);
if sum(pL2-[1 2 3])>0, fprintf('bad\n'); else fprintf('good\n'); end
% permuation
L2 = [2 3 3 1 1 1]; [pL2,Qe,acc]=labels_accuracy(L1,L2);
if sum(pL2-[3 1 2])>0, fprintf('bad\n'); else fprintf('good\n'); end

% differing K's
L1 = [1 2 3]; L2 = [1 2 2]; labels_accuracy(L1,L2); labels_accuracy(L2,L1);
