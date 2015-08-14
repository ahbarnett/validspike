function m = mergeclips(A,B)
% MERGECLIPS - merge two sets of variable- or fixed-length clip arrays X
%
% m = mergeclips(A,B) produces a variable-length clip struct m from A & B,
%  in that order.
%  A, B may each be either 3d (fixed-length) clips or variable-length clip
%  structs.
%
% The fixed-length format is
%    X (double M*Nt*Ns) - signal array. M channels, Nt time points
%        per clip, and Ns clips.
%
% The variable-length clip struct m (or inputs A and/or B) has fields:
%    X (double M*Ttot) - concatenated signal array. M channels
%    Ttot (int) - total time points = sum(Ts)
%    Ts (int 1*Ns) - clip lengths in time points
%    tptr (int 1*Ns) - indices (hence 1-indexed) of starts of each clip
%    Ns (int) - number of clips
%
% m = mergeclips(X) produces variable-length clip struct from 3D array X.

% todo: think about if A, B huge in RAM.
%
% Barnett 2/19/15

if nargin<2, B.X = []; B.Ts = []; end  % dummy B so A converted to m

if isstruct(A) && isstruct(B)
  m.X = [A.X B.X]; m.Ts = [A.Ts B.Ts];
elseif isstruct(A) && ~isstruct(B)
  [M Nt Ns] = size(B);
  m.X = [A.X reshape(B,[M Nt*Ns])];
  m.Ts = [A.Ts Nt*ones(1,Ns)];
elseif ~isstruct(A) && isstruct(B)
  [M Nt Ns] = size(A);
  m.X = [reshape(A,[M Nt*Ns]), B.X];
  m.Ts = [Nt*ones(1,Ns) B.Ts];
else
  [M NtA NsA] = size(A); [M NtB NsB] = size(B);
  m.X = [A B]; m.Ts = [NtA*ones(1,NsA) NtB*ones(1,NsB)];
end
m.Ttot = sum(m.Ts);                 % less to go wrong than if case-by-case
m.Ns = numel(m.Ts);
m.tptr = [1 1+cumsum(m.Ts(1:m.Ns-1))];
