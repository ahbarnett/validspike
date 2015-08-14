function [colperm Qperm] = bestcolpermconfmat(Q)
% BESTCOLPERMCONFMAT  column permutation making confusion matrix most diagonal.
%
% [colperm Qperm] = bestcolpermconfmat(Q) given a K1-by-K2 matrix Q of
%   non-negative entries returns colperm, the best permutation of the column
%   labels of Q that makes the matrix have the largest sum of diagonal entries,
%   and the resulting permuted matrix Qperm = Q(:,invcolperm) where invcolperm
%   is the inverse of the permutation. Note that colperm is the permutation
%   to apply to the 2nd list of labels that generated the confusion matrix.
%
% Notes: 1) The Hungarian algorithm is used, which is polynomial in K1 and K2
% 2) There can be zero rows or columns, ie, unused label classes.
% 3) A true permutation of labels 1...K2 is always returned.
%
% If you wish to find a row permutation, send in Q' and do Qperm=Qperm' after.
%
% With no arguments, a self-test is done.

% Barnett 8/14/15, supercedes getbestshuffling by Magland/Barnett
if nargin==0, test_bestcolpermconfmat; return; end
[K1 K2] = size(Q);
M = Hungarian(-Q);      % matching matrix of zeros and ones
colperm = zeros(1,K2);  
next_unused_label = K1+1;                      % trick from Jeremy
for j=1:K2
  row = find(M(:,j)==1);
  if isempty(row)
    colperm(j) = next_unused_label;
    next_unused_label = next_unused_label + 1;
  else
    colperm(j) = row;
  end
end
[~,invcolperm] = sort(colperm);
Qperm = Q(:,invcolperm);

%%%%%%
function test_bestcolpermconfmat

disp('warm-up tests...')
[p Qp] = bestcolpermconfmat([0 1 0; 0 0 1; 1 0 0])  % K1=K2
if sum(p-[3 1 2])==0 && sum(diag(Qp))==3, fprintf('good\n');
else fprintf('bad\n'); end
[p Qp] = bestcolpermconfmat([0 1 0 0 ;0 0 1 0])  % K2>K1
[p Qp] = bestcolpermconfmat([0 1 0 0;0 0 1 0]')  % K1>K2

% variants of tests inherited from getbestshuffling
disp('testing K1=K2...')
N = 1e3;  % number of cases
K = 5;    % how many label types
a = randi(K,1,N);
je = randperm(K);   % the true shuffling
b = je(a);
j = bestcolpermconfmat(confusion_matrix(b,a));
fprintf('noise-free case: errors in the shuffling = %d\n',numel(find(j~=je)))

i = rand(1,N)<0.05; % some indices to mess up in b
b(i) = randi(K,1,numel(find(i))); % mess them up
j = bestcolpermconfmat(confusion_matrix(b,a));
fprintf('noisy case: errors in the shuffling = %d\n',numel(find(j~=je)))

disp('testing K1>K2...')
[p Qp]=bestcolpermconfmat(confusion_matrix([1,1,1,2,2,2,3,3,3,4,4,4,5,5,5],[1,1,1,2,2,3,3,3,3,4,4,4,2,2,2]))
'note 2 is mapped to 5 (K2+1)'
