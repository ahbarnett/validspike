function overview_sorted_clips(X,L,title,n,o,varargin)
% OVERVIEW_SORTED_CLIPS - show sensible samples of clips of each label
%
% overview_sorted_clips(X,L)
% overview_sorted_clips(X,L,title)
% overview_sorted_clips(X,L,title,n,fig,opts)
%   controls opts such as opts.equal = 1 (equalize populations), 0 (sample
%                         according to pops)
%                         opts.vs = vertical scale (default [])
% overview_sorted_clips(X,L,title,n,fig,opts,pssopts)
%   controls n the number of clips shown, the figure handle to plot in,
%   and other options passed to plot_spike_shapes.

% todo: doc
% todo: overlay multiple spikes of same type, as most other ppl do?
%
% Barnett 6/10/15. equalized pops 6/14/15

if nargin<3, title=''; end
if nargin<4 || isempty(n), n = 100; end  % default total # clips to show
if nargin<5, o = []; end
if ~isfield(o,'equal'), o.equal = 1; end  % equalize pops
if ~isfield(o,'vs'), o.vs = []; end

[M T N] = size(X);
if o.equal      % equalize pops
  j = []; K = max(L); nk = diff(round((0:K)/K*n)); % approx equal pops
  for k=1:K
    i = find(L==k); l = randperm(min(numel(i),nk(k)),nk(k));
    j = [j i(l)];
  end
  L = L(j);
else            % sample according to pops
  j = randperm(N,n);
  [L,i] = sort(L(j)); j = j(i);            % sort labels in order
end
p = []; for i=1:n, p(i).l = L(i); end    % build param struct array
plot_spike_shapes(X(:,:,j), title, o.vs, p, varargin{:});
