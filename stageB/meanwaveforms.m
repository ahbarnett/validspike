function W = meanwaveforms(X,L)
% MEANWAVEFORMS - compute mean waveforms W from clips X and labels L
%
% W = meanwaveforms(X,L) computes mean waveforms W (M-by-T-by-K) given X clips
%   (M-by-T-by-N) and L labels (1-by-N). K is determined by the max label.

% Barnett 7/1/15

[M T N] = size(X);
K = max(L);
W = nan(M,T,K);
for j=1:K, W(:,:,j) = mean(X(:,:,L==j),3); end
