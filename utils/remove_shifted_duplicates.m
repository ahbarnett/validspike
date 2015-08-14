function i = remove_shifted_duplicates(W,o)
% REMOVE_SHIFTED_DUPLICATES - identify and remove which k in W are duplicates
%
% i = remove_shifted_duplicates(W,o)
%  returns index list to keep. Prefers small indices (since we assume the W_k
%  were in decreasing population order)
% Inputs:
%  W - M*T*K waveform array
%  o - (optional) controls following:
%      o.eps (default 0.3) max allowed error relative to norm(W_k)
%              to count as duplicate
% Output:
%  i - index list (or logical over 1...K)
%      such that W(:,:,i) is then a cleaned-up waveform array
%
% Barnett 5/19/15

if nargin<2, o = []; end
if ~isfield(o,'eps'), o.eps = 0.3; end   % adjust this better!

[M T K] = size(W);
trange = -round(T/2):round(T/2);  % range of timeshifts to use
i=1;
for j=2:K
  wj = W(:,:,j); errmax = o.eps * norm(wj(:));
  killit = 0;
  for k=i         % check against all previous unique k and timeshifts
    for tsh=trange
      ii = max(0,-tsh)+1:min(0,-tsh)+T;
      dw = W(:,ii,k) - W(:,ii+tsh,j);
      err = norm(dw(:));
      if err<errmax, killit=1; end
    end
  end
  if ~killit, i = [i j]; end % append to accepted index list
end
