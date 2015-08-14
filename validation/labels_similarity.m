function [Q acc] = labels_similarity(a,b,o)
% LABELS_SIMILARITY  similarity metrics between two lists of labels
%
% Q = labels_similarity(a,b) with a,b vectors of integer labels, with a's entries
%  in [1,..,Ka] (or unclassified, either <1 or NaN), and b's entries in
%  [1,..,Kb] (or unclassified), returns a (Ka+1)-by-(Kb+1) extended confusion
%  matrix Q, where the last row and column correspond to unclassified in a and b
%  respectively. Entry Q(i,j) counts the number of entries labeled i in a and
%  j in b, for i=1..Ka, and j=1..Kb.  i=Ka+1 and j=Kb+1 are unclassified.
%  a is designed to be "truth" while b the "test" labels.
%
%  Accuracy metrics per label type are reported.
%
%  No label permutation (searching for best permuation) is done.
%
% [Q acc] = labels_similarity(a,b) also returns metric vectors in the struct acc,
%   (here K = min(Ka,Kb)), namely fields:
%    p = 2Q_kk / sum_j(Q_kj+Q_jk), for k=1..Ka, Jeremy's accuracy indicator
%        (is zero for k=Kb+1..Ka, if entries exist). This is f_k in the paper.
%   tI = type-I error rates (ie missed, counting wrong and unclassified together)
%        length Ka. Entries k=Kb+1,..,Ka, if they exist, are 1.
%    m = type-I error rates counting only unclassified in the test set b (ie
%        fraction completely missed as labeled types), length Ka
%  tII = type-II error rates (ie false positives, counting wrong and unclassified
%        together), length Kb. Entries k=Ka+1,..,Kb, if they exist, are 1.
%    f = type-II error rates counting only unclassified in the truth set a,
%        ie events that were not classified in the truth set. Length Kb.
%    goodfrac = overall fraction of events that are correctly matched
%        (ratio of diagonal sum of Q to total # events)
%
% Thus m counts a subset of tI; f counts a subset of tII.
%
% [...] = labels_similarity(Q) treats input as an extended confusion matrix,
%   and proceeds as above.
%
% [...] = labels_similarity(a,b,o)
% [...] = labels_similarity(Q,[],o)
%  controls various options such as:
%      o.verb - verbosity (0=silent, 1=text, default)
%
% See also: LABELS_ACCURACY which serves as a test for this function

% Barnett 12/11/14. rectangular case with unclass also handled 1/9/15
% 6/12/15 extended for confusion matrix. 6/16/15 redo definitions
% 8/14/15 used confusion_matrix, renamed A to Q

if nargin<3, o=[]; end
if ~isfield(o,'verb'), o.verb = 1; end
if nargin==1 || isempty(b)       % conf mat supplied
  Q = a; [Ka,Kb]=size(Q); Ka=Ka-1; Kb=Kb-1; Ns = sum(Q(:)); % Ns = # events
else
  Ns = numel(a);
  if numel(b)~=Ns, error('lists a and b must have same lengths'); end
  Ka = max(a); Kb = max(b);
  a(find(a<1 | isnan(a))) = Ka+1;  % standardize unclassified entries of a
  b(find(b<1 | isnan(b))) = Kb+1;  % ditto b
  Q = confusion_matrix(a,b,[Ka+1,Kb+1]);   % compute extended confusion matrix
end
clear a b                          % everything dep on conf mat Q alone now
K = min(Ka,Kb);                    % largest poss number of matches
na = sum(Q,2); nb = sum(Q,1);      % populations
acc.tI = ones(1,Ka); acc.m = acc.tI;           % missing
for k=1:K, acc.tI(k) = 1 - Q(k,k)/na(k); acc.m(k) = Q(k,Kb+1)/na(k); end
acc.tII = ones(1,Kb); acc.f = acc.tII;         % false positives
for k=1:K, acc.tII(k) = 1 - Q(k,k)/nb(k); acc.f(k) = Q(Ka+1,k)/nb(k); end
acc.p = zeros(1,Ka);                   % Jeremy's
for k=1:K, acc.p(k) = 2*Q(k,k) / (na(k)+nb(k)); end
acc.goodfrac = sum(diag(Q(1:K,1:K))) / Ns;

if o.verb % verbose
  fprintf('extended confusion matrix:\n')
  for i=1:size(Q,1), fprintf('%d\t',Q(i,:)), fprintf('\n'); end
  disp('type-I err rates (missing, wrong+unclass):'), fprintf('%8.3g ', acc.tI),fprintf('\n')
  if max(acc.m)>0
    disp('  m err rates (missing, unclass only):'), fprintf('%8.3g ', acc.m),fprintf('\n')
  end
  disp('type-II err rates (false pos, wrong+unclass):'),fprintf('%8.3g ', acc.tII),fprintf('\n')
  if max(acc.m)>0
    disp('  f err rates (false pos, unclass only):'), fprintf('%8.3g ', acc.f),fprintf('\n')
  end
  disp('accuracies f_k for true labels:'), fprintf('%8.3g ', acc.p),fprintf('\n')
  fprintf('fraction correctly matching of total events: %.3g\n',acc.goodfrac)
end
