function [jt l s] = locvalidmins(S,nlps)
% LOCVALIDMINS  find local valid minima in the S score matrix (MEX, devel)
%
% [jt l s] = locvalidmins(S,nlps)
%
% outputs:
%  jt - time indices (not times) of local minima found
%  l  - labels of same
%  s  - S values of same (without nlps contrib)
%
% Note: a single-core development MEX interface

[K Nsh] = size(S);
if nargin<2, nlps = zeros(1,K); end

jt = zeros(1,Nsh); l=jt; s=jt;  % alloc: assume <=1 minimum per time shift pt

mex_id_ = 'locvalidmins(i double[], i int, i int, i double[], io int[], io int[], io double[], o int[x])';
[jt, l, s, Nmin] = gf(mex_id_, S, K, Nsh, nlps, jt, l, s, 1);
% int[1] to get 1 integer out

jt = jt(1:Nmin)+1;  % 1-indexed for matlab
l = l(1:Nmin); s = s(1:Nmin); % truncate to correct size


% ------------------------------------------------------------------------------
