function powerspec(d)
% POWERSPEC - show power spectrum of raw EC dataset across all channels
%
% powerspec(d)

% Barnett 2/26/15

if nargin<1, test_powerspec; return; end

[M N] = size(d.A);
T = N*d.dt;  % total time
df = 1/T;    % freq grid...
if mod(N,2)==0, f=df*[0:N/2 -N/2+1:-1];else, f=df*[0:(N-1)/2 -(N-1)/2:-1]; end
i = f>=0;  % indices to keep

F = fft(d.A,[],2)/N;
F = sum(abs(F).^2,1);
figure; plot(f(i),F(i),'-'); xlabel('f (Hz)'); ylabel('spectral density')
title(d.name)
%%%

function test_powerspec
d = loaddemodata;
powerspec(d);
