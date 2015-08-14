function [m1 m2] = signal_moments(X)
% SIGNAL_MOMENTS - extract 1st, 2nd moments of each of multichannel events
%
% function [m1 m2] = signal_moments(X) computes centroid along 2nd index,
%  weighting 1st index (channels) equally, for each 2D array indexed by 3rd
%  index.
%
% Inputs:
%  X - M*T*N signal array (M channels, T time points, N events)
% Outputs:
%  m1 - (1*N) mean time of each event, in time index units (1 to T)
%  m2 - (1*N) time variance of each event about m1, in time index units (1 to T)

% Barnett 2/10/15

if nargin==0, test_signal_moments; return; end
[M T N] = size(X);
m1 = nan(1,N); m2 = m1;
for n=1:N
  m0 = sum(sum(X(:,:,n)));
  if m0~=0
    m1(n) = sum(sum(bsxfun(@times, X(:,:,n), 1:T))) / m0;
    m2(n) = sum(sum(bsxfun(@times, X(:,:,n), ((1:T)-m1(n)).^2))) / m0;
  end
end

function test_signal_moments
N=10;
M = 3;
T = 30;
X = nan(M,T,N);
% make successively wider square pulses. note n=1 gives zero sig...
for n=1:N, X(:,:,n) = kron(abs((1:T)-T/2)<n-1,0.7*ones(M,1)); end
[m1 m2] = signal_moments(X);
fprintf('should both be zero:\n')
norm(m1(2:end)-T/2)   % compare true means
norm(m2(2:end)-(0:N-2).*(1:N-1)/3)  % compare true variances n(n+1)/3
