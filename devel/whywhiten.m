% testing prewhitening ideas on synthetic noise data.
% Barnett + Magland 1/27/16

clear all; close all

N = 1e6;  % samples
M = 10;   % channels

Y = randn(M,N);

%Y = Y + repmat(randn(1,N),[M 1]);  % cm noise
Y = Y + repmat(cos(pi*(1:N)),[M 1]);  % common-mode noise, known osc func.

Y(1,100:150) = Y(1,100:150) + 5;   % add a 5-sigma "spike" to see if preserved

spikespy({Y(:,1:300),'Y'});
Y = Y - repmat(mean(Y,2),[1 N]);

Sig = Y*Y'/N;   % est covar
figure; imagesc(Sig); colorbar; caxis([0 max(Sig(:))]); title('Sig')

[U D] = eig(Sig); D = sqrt(D);  % make D singular values

Qj = U*inv(D)*U';    % jeremy, symm, preserves channel identities
%Qj = inv(D)*U';    % tweaked jeremy that rotates into pca bases (mixes chans)
Qa = inv(chol(Sig))';   % note the transpose is crucial here!
figure; imagesc(Qj); colorbar; title('Qj prewhitening matrix');
figure; imagesc(Qa); colorbar; title('Qa prewhitening matrix');

Ytj = Qj*Y;
Yta = Qa*Y;

% should all be Id:
Sigtj = Ytj*Ytj'/N
figure; imagesc(Sigtj); colorbar; caxis([0 max(Sigtj(:))]); title('Sigtj')
Sigta = Yta*Yta'/N
figure; imagesc(Sigta); colorbar; caxis([0 max(Sigta(:))]); title('Sigta')

%U*D^2*U'
spikespy({Yta(:,1:300),'Yta'});
spikespy({Ytj(:,1:300),'Ytj'});

% Conclusions: Qa method isn't symmetric across channels, still mixes most of
% them.
% Qj = U*inv(D)*U' seems best, shrinks common-mode noise by 1/M, preserves
% chan idents.
% Qj = inv(D)*U' moves common-mode noise to a single channel, but mixes.

% check x-t corr matrices after doing this, in spike-free clips & w/ whole Y.

