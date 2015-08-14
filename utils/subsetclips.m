function n = subsetclips(m,i)
% SUBSETCLIPS - make variable-length clips object from subset of another
%
% n = subsetclips(m,i) makes n a variable-length clips object from the indices
%  i in variable-length clip object m.
%
% See also (for format of m and n): mergeclips, detectevents

% Barnett 2/19/15

n = m;   % copies over unknown fields
n.Ts = m.Ts(i);
n.Ns = numel(n.Ts);
n.Ttot = sum(n.Ts);
n.tptr = [1 1+cumsum(n.Ts(1:end-1))];
j = nan(1,n.Ns);    % will index which cols to get from m.X
if islogical(i), i = find(i); end   % insure a list of indices
p = 0;     % pointer
for c=1:numel(i)       % make index list
  k=1:m.Ts(i(c)); j(p+k) = m.tptr(i(c))+k; p = p + m.Ts(i(c));
end
n.X = m.X(:,j);   % move the data via index list
