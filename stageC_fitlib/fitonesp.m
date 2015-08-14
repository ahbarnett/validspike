function [p Jbest Fb] = fitonesp(wf,Y,noi,opts)
% FITONESP  maximum likelihood fit for single spike in single clip (MEX, devel)
%
% [p Jbest Fb] = fitonesp(wf,Y,noi,opts)
%
% Inputs:
%  wf - waveforms struct containing: W = M*T*K waveforms
%                                    d = raw EC data struct (d.A not needed)
%                                    fac = upsampling ratio
%                                    freqs = 1*K firing freq estimates
%  Y - M*Nt data array to fit (M channels, Nt time points)
%  noi - noise model struct, minimally noi.eta = sqrt(variance)
%  opts :  (optional) struct controlling parameters:
%          opts.nlps (default 20 for all types) neg log prior for each spike
%                         type 1..K
%          opts.locflag (default 0): 0 - use std NLL calc, 1 - local update.
%
% Outputs:
%  p - best-fit param struct, eg p.t time, p.a ampl, p.m = identity
%  Jbest - best-fit objective func = neg log lik
%  Fb - model output at best-fit params p
%
% See also: fitonespike.m, for which this is a drop-in replacement

if nargin<4, opts=[]; end
if wf.fac~=floor(wf.fac), error('wf.fac must be positive integer!'); end
W = wf.W; fac = wf.fac; % enforcing int32 for ints in C breaks MEX interface!
[M T K] = size(W);
if mod(T,2)==0, warning('waveform T must be odd'); end
[MY Nt] = size(Y);
if MY~=M, error('Y must have same number of channels M as W'); end
eta = noi.eta;
if isfield(noi,'covar'), warning('covar noise model not avail in MEX!'); end
tpad = 1.0;         % default edge padding
if isfield(opts,'nlps'), nlps = opts.nlps, else nlps = zeros(1,K); end
if numel(nlps)==1, nlps = nlps*ones(1,K); end   % duplicate across spike types
if isfield(opts,'locflag'), locflag = opts.locflag; else locflag = 0; end
srt = sum(Y.^2,1);        % initialize squared resid per time-pt

TK = T*K; W = reshape(W,[M TK]);  % mwrap can't do 3d arrays

% gave up using C function output; reverted to void type & output via pointer.
% Seems like can't have math expressions in sizes
mex_id_ = 'fitonesp(i double[], i int, i int, i int, i int, o int[x], o double[x], o double[x], i int, i double[], i double, i double, o double[xx], o double[x], i double[], i int, io double[])';
[lb, tb, ab, Fb, Jbest, srt] = sf(mex_id_, W, M, T, K, fac, Nt, Y, eta, tpad, nlps, locflag, srt, 1, 1, 1, M, Nt, 1);

p.l = lb; p.t = tb; p.a = ab;   % make output param struct


%=============================================================================
