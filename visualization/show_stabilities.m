function show_stabilities(fhat,fsam,o)
% SHOW_STABILITIES - plot stability metric distributions per spike type
%
% function show_stabilities(fhat,fsam) plots fhat for each label k=1...K as
%  large blobs, and the full set of samples and quantiles in fsam, if present.
%
% fhat (1-by-K)
% fsam should have samples down each column, with K columns (if empty, no
%   samples or quantiles 
% n is populations
%
% function show_stabilities(fhat,fsam,o)
%  o.ylab
%  o.pops - if given (length-K vector), move the f_k bar labels up and add
%    populations
%  o.fig - if present, add to this figure handle
%  o.blobcolor - (1x3) if present, use as color for blobs plotting

% Barnett 6/26/15
if nargin<2, fsam=[]; end
if nargin<3, o=[]; end
if ~isfield(o,'ylab'), o.ylab = 'f_k'; end
if ~isfield(o,'blobcolor'), o.blobcolor = [0 0 0]; end

K = numel(fhat);
if ~isempty(fsam) & size(fsam,2)~=K
  if size(fsam,1)==K, fsam = fsam'; else error('fsam wrong size'); end
end
if ~isempty(fsam) & size(fsam,1)>1
  q = quantile(fsam,[.25 .75]);         % outputs 2 by K
else, q = []; end

if isfield(o,'fig'), figure(o.fig); hold on; else figure; end
wid = .25; col = [.7 1 1];     % half-width and color of quantile bars
if ~isempty(q)
  h=patch(repmat(wid*[-1 1 1 -1]',[1 K])+kron(ones(4,1),1:K), [q(1,:);q(1,:);q(2,:);q(2,:)], col);
  hold on
end
if ~isempty(fsam), plot(1:K, fsam, 'b.'); end
plot(1:K, fhat, '.','markersize',30,'color',o.blobcolor);
xlabel('k'); ylabel(o.ylab); axis([.5 K+.5 0 1]);
if o.blobcolor==[0 0 0], ysh = 0; else, ysh=0.1; end         % y text offset

if isfield(o,'pops')  % show both f_k bar and pops...
  for k=1:K
    text(k-0.3,0.2,sprintf('%d%%',round(100*fhat(k))));
    text(k-0.3,0.1,sprintf('%d',round(o.pops(k))));
  end
  text(K+1-0.3,0.2,'$$\bar{f}_k$$','interpreter','latex')
  text(K+1-0.3,0.1,'$$n_k$$','interpreter','latex')
else                  % just f_k bar...
  for k=1:K
    text(k-0.3,0.1+ysh,sprintf('%d%%',round(100*fhat(k))),'color',o.blobcolor);
  end
end
set(gca,'xtick',1:K); box off
