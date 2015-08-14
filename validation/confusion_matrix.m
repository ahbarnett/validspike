function Q = confusion_matrix(a,b,siz)
% CONFUSION_MATRIX  fill K1-by-K2 (not-extended) confusion matrix given labels.
%
% Q = confusion_matrix(a,b) where a and b are row or column vectors of equal
%  length N containing natural number entries, returns a K1-by-K2 matrix,
%  where K1=max(a) and K2=max(b), whose ij entry gives the number of values k
%  in 1...N such that a(k)=i and b(k)=j.
%
% Q = confusion_matrix(a,b,[K1 K2]) overrides K1 and K2, fixing the size of Q.
%
% With no arguments, a self-test is done.

% Barnett pulled out 8/14/15
if nargin==0, test_confusion_matrix; return; end
if nargin<3
  K1 = max(a); K2 = max(b);
else
  K1 = siz(1); K2 = siz(2);
end
Q = zeros(K1,K2);
for i=1:K1
  for j=1:K2
    Q(i,j) = sum((a(:)==i).*(b(:)==j));
  end
end

%%%%%%%
function test_confusion_matrix
Q = confusion_matrix([1],[1 2])
Q = confusion_matrix([1 2],[1])
Q = confusion_matrix([2 1 1],[1 2 1])
Q = confusion_matrix([2 1 1],[1 2 1],[3 4])
