function [t l Y] = gluttonstuffme(wf,Y,noi,nlps,o)
% GLUTTONSTUFFME  single-core version of multiglutton (MEX, devel)
%
% function [t l Y] = gluttonstuffme(wf,Y,noi,nlps,o)
%
% todo: doc. add S output
%
% See for usage: sort_filtereddata.m
%
% Note: a single-core development MEX interface

if wf.fac~=floor(wf.fac), error('wf.fac must be positive integer!'); end
W = wf.W; fac = wf.fac;
[M T K] = size(W);
if mod(T,2)==0, warning('waveform T must be odd'); end
[M Nt] = size(Y);
eta = noi.eta;
if nargin<4 || isempty(nlps), nlps = zeros(1,K); end
if nargin<5, o = []; end
if isfield(o,'skip'), skip = o.skip; else skip = fac; end  % conservative
if isfield(o,'tpad'), tpad = o.tpad; else tpad = 2.0; end
if isfield(o,'gamma'), gamma = o.gamma; else gamma = 10.0; end
if isfield(o,'verb'), verb = o.verb; else verb = 0; end

maxNs = Nt;  % alloc outputs to max size

mex_id_ = 'gluttonstuffme(i double[], i int, i int, i int, i int, io double[], i int, i double, i double, i int, i double, i double[], o double[x], o int[x], o double[x], o int[x], i int)';
[Y, t, l, a, Ns] = gf(mex_id_, W, M, T, K, fac, Y, Nt, tpad, eta, skip, gamma, nlps, verb, maxNs, maxNs, maxNs, 1);

t = t(1:Ns); l = l(1:Ns); a = a(1:Ns);
t = t(:)'; l = l(:)'; a = a(:)';  % make row vecs

% ------------------------------------------------------------------------------
