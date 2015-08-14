function keep = minisgammaloc(jt,l,val,tsh,j,gamma)
% MINISGAMMALOC  test if local minimum is global over an interval (MEX, devel)
%
% keep = minisgammaloc(jt,l,val,tsh,j,gamma)
% Inputs:
%  jt - list of 1-indexed indices to tsh timeshifts
%  l - list of types (each element is in 1..K)
%  val - list of NLP-corrected S values
%  tsh - list of timeshifts
%  j - index in 1...numel(jt) to check (1-indexed)
%  gamma - width in time-samples of region to test global minimum over
%
% Note: a single-core development MEX interface

j = j-1; jt=jt-1; % convert to 0-indexing
N = numel(jt);

mex_id_ = 'o int = minisgammaloc(i int[], i int[], i double[], i double[], i int, i int, i double)';
[keep] = gf(mex_id_, jt, l, val, tsh, N, j, gamma);


% ------------------------------------------------------------------------------
