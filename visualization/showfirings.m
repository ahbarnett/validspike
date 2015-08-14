function showfirings(t,l)
% SHOWFIRINGS - plot spikes as vertical lines to show timing
%
% showfirings(t,l)
% Inputs:
%  t - vector of firing times
%  l - vector of firing labels in 1..K
%
% Barnett 2/23/15

t = t(:)'; N=numel(t); % for plotting
K = max(l);
co = get(gca,'colororder'); Nco = size(co,1); % standard color ordering
for k=1:K
  c = co(mod(k-1,Nco)+1,:); % 1x3 color vector
  tk = t(find(l==k));
  plot([tk;tk], repmat(k+[-.5;.5],[1 numel(tk)]), '-','color',c);
  hold on
end
set(gca,'ytick',[1:K]);
set(gca, 'Ticklength', [0 0])
axis tight
