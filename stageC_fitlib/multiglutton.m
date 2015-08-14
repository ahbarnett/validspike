function [t l R] = multiglutton(wf,Y,noi,nlps,o)
% MULTIGLUTTON  thresholdless fit spike times, labels in timeseries (MEX/OpenMP)
%
% function [t l] = multiglutton(wf,Y,noi,nlps,o) returns times t and labels l of
%  all spikes fitted in multi-channel time-series Y, given waveform object wf
%  (containing wf.W waveforms, etc), noise parameters in noi, prior firing
%  probabilies in nlps, and possible options in o.
%
% function [t l R] = multiglutton(wf,Y,noi,nlps,o) also returns residual R;
%  however this residual has errors within T of the chunk breakpoints. Therefore,
%  not recommended!
%
% See for usage: sort_filtereddata

% todo: better doc
% Barnett May 2015

% same code as gluttonstuffme...
if wf.fac~=floor(wf.fac), error('wf.fac must be positive integer!'); end
W = wf.W; fac = wf.fac;
[M T K] = size(W);
if mod(T,2)==0, warning('waveform T must be odd'); end
[M Nt] = size(Y);
Y = double(Y);        % essential in case comes in as float or int
eta = noi.eta;
if nargin<4 || isempty(nlps), nlps = zeros(1,K); end
if nargin<5, o = []; end
if isfield(o,'skip'), skip = o.skip; else skip = fac; end  % conservative
if isfield(o,'tpad'), tpad = o.tpad; else tpad = 2.0; end
if isfield(o,'gamma'), gamma = o.gamma; else gamma = 10.0; end
if isfield(o,'verb'), verb = o.verb; else verb = 0; end

maxNs = Nt; % alloc outputs to max size
wantR = double(nargout>2);
if wantR, warning('residual R from multiglutton has errors at chunk ends!'); end

mex_id_ = 'multiglutton(i double[], i int, i int, i int, i int, io double[], i int, i double, i double, i int, i double, i double[], o double[x], o int[x], o double[x], o int[x], o double[xx], i int, i int)';
[Y, t, l, a, Ns, R] = gf(mex_id_, W, M, T, K, fac, Y, Nt, tpad, eta, skip, gamma, nlps, wantR, verb, maxNs, maxNs, maxNs, 1, M, Nt);

t = t(1:Ns); l = l(1:Ns); a = a(1:Ns);
t = t(:)'; l = l(:)'; a = a(:)';  % make row vecs
