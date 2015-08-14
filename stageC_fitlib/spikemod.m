function F = spikemod(wf, p, Nt)
% SPIKEMOD  forward model for single clip given spike labels, times, ampls (MEX)
%
% F = spikemod(wf, p, Nt)
%
% Inputs:
%  wf - waveform object:
%    W     = M*NT*K upsampled centered waveforms (NT assumed odd)
%    fac   = factor to downsample (1, or factor which upsampling was done)
%  p - parameters object:
%    l     = 1*Ns labels in 1,..,K (Ns is the number of spike events)
%    t     = 1*Ns time alignments peak, measured in output samps rel to 1st samp
%    a     = 1*Ns amplitudes (typically close to 1; if empty, all set to 1)
%  Nt    = # output time samples
%
% Outputs:
%  F     = M*Nt signal window
%
% This is same interface and function as synthesis/spikemodel.m mostly for
% testing of MEX.

if wf.fac~=floor(wf.fac), error('wf.fac must be positive integer!'); end
W = wf.W; fac = wf.fac;
[M T K] = size(W);
if mod(T,2)==0, warning('waveform T must be odd'); end
if Nt<1, error('Nt must be positive!'); end
if isempty(p) || numel(p.l)==0, Ns=0;        % how many spikes in this event?
else Ns = numel(p.l);
end
if ~isfield(p,'a') || isempty(p.a), p.a = 1+0*p.t; end       % default ampl
if numel(p.t)~=Ns || numel(p.a)~=Ns, error('t and a must be same size as l');end
l = p.l; t = p.t; a = p.a;  % since mwrap can't handle struct elements
subF = 0;   % start with F zero (subF=1 would subtract from whatever Y was)
iran = [0 0];   % for support of F index range output

mex_id_ = 'spikemod(i double[], i int, i int, i int, i int, i int, i int[], i double[], i double[], i int, o double[xx], i int, io int[x])';
[F, iran] = sf(mex_id_, W, M, T, K, fac, Ns, l, t, a, Nt, subF, iran, M, Nt, 2);



%=============================================================================
