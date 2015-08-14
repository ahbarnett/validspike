function S = fillscore(wf,Y,tsh,noi,o)
% FILLSCORE  fill S score matrix given timeseries and waveforms (MEX, devel)
%
% S = fillscore(wf,Y,tshnoi,opts)
% Inputs:
%  wf  - waveform object containing W (M*T*K) and fac (beta upsampling param)
%  Y   - data matrix (M*Nt)
%  tsh - list of timeshifts in signal units, relative to first time point in Y
%  noi - noise model with property noi.eta setting iid Gaussian std deviation.
%  opts - (optional) controls options such as
%        opts.skip - 1,2,.. calls smartscore if >1 (overrides shflags)
%        opts.shflags - boolean giving which cols of S to calc, or set to Nans
% Outputs:
%  S   - K*numel(tsh) score matrix
%
% Note: this is single-core, and a development MEX interface

if wf.fac~=floor(wf.fac), error('wf.fac must be positive integer!'); end
W = wf.W; fac = wf.fac;
[M T K] = size(W);
if mod(T,2)==0, warning('waveform T must be odd'); end
[M Nt] = size(Y);
eta = noi.eta;
Nsh = numel(tsh); numelS = K*Nsh;  % size of S
if nargin<5, o = []; end
if isfield(o,'skip'), skip = o.skip; else skip = fac; end     % default skip

if skip>1
  mex_id_ = 'smartscore(i double[], i int, i int, i int, i int, i double[], i int, i double[], i int, o double[xx], i double, i int)';
[S] = gf(mex_id_, W, M, T, K, fac, Y, Nt, tsh, Nsh, eta, skip, K, Nsh);
else    % obsolete interface, for checking fillscore...
  if isfield(o,'shflags'), shflags = double(o.shflags); else shflags=ones(1,Nsh); end % default is calc for all shift indices
  mex_id_ = 'fillscore(i double[], i int, i int, i int, i int, i double[], i int, i double[], i int, o double[xx], i double, i int[])';
[S] = gf(mex_id_, W, M, T, K, fac, Y, Nt, tsh, Nsh, eta, shflags, K, Nsh);
end		     


% ------------------------------------------------------------------------------
