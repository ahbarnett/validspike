function X = freqfilter(X,fs,flo,fhi)
% FREQFILTER - filter rows of X using smooth roll-offs in Fourier space
%
% X = freqfilter(X,fs,flo,fhi) where X is M*N matrix returns a matrix of same
%  size where each row has been filtered in a bandpass filter with soft rolloffs
%  from flo to fhi (in Hz). fs is the sampling freq in Hz. If flo is empty,
%  a lo-pass is done; if fhi is empty, a hi-pass.
%
% Note: MATLAB's fft is single-core even though claims multicore, for fft(X,[],2)
%  only. fft(X) is multicore.
%
% Hidden parameters: fwid - width of roll-off (Hz). todo: make options.
%
% todo: This could act on EC data object instead?
%
% Barnett 11/14/14.
% 3/11/15: transpose to fft cols (is multicore), blocking (was slower!)

if nargin<1, test_freqfilter; return; end
if nargin<3, flo = []; end
if nargin<4, fhi = []; end
if ~isempty(fhi) & ~isempty(flo) & fhi<=flo, warning('fhi<=flo: are you sure you want this??'); end

[M N] = size(X);

Nbig = inf; %2^22;  % do it in blocks (power of 2 = efficient for FFT)
if N>2*Nbig
  pad = 1e3;   % only good if filters localized in time (smooth in k space)
  B = Nbig-2*pad;  % block size
  Xpre = zeros(M,pad);
  ptr = 0;
  while ptr+B<N  % so last block can be at most 2B wide
    if ptr+2*B<N, j = ptr+(1:B+pad);
    else, j = ptr+1:N; end  % final block is to end of data
    Y = [Xpre X(:,j)];  % for all but last time Y has Nbig cols, fast
    Xpre = X(:,ptr+B+(-pad+1:0));  % before gets overwritten get next left-pad
    size(Y)
    Yf = freqfilter(Y,fs,flo,fhi);
    X(:,j) = Yf(:,pad+1:end); % right-pad region will get overwritten, fine
    ptr = ptr + B;
  end
  return              % !
end                   % todo: figure out why this was slower than unblocked
                      % even for Nbig = 2^22

T = N/fs;  % total time
df = 1/T; % freq grid...
if mod(N,2)==0, f=df*[0:N/2 -N/2+1:-1];else, f=df*[0:(N-1)/2 -(N-1)/2:-1]; end

a = ones(size(f));
fwidlo = 100;       % roll-off width (Hz). Sets ringing timescale << 10 ms
if ~isempty(flo), a = a .* (1+tanh((abs(f)-flo)/fwidlo))/2; end
fwidhi = 1000;       % roll-off width (Hz). Sets ringing timescale << 1 ms
if ~isempty(fhi), a = a .* (1-tanh((abs(f)-fhi)/fwidhi))/2; end
X = ifft(bsxfun(@times, fft(X'), a'))'; % filter: FFT fast along fast
% storage direction, transposing at input & output

%%%%%
function test_freqfilter
M=1; N=1e6; X = zeros(M,N); X(:,N/2) = 1; %X = randn(M,N); % dummy data spike
fs = 2e4;
T = N/fs; df = 1/T; f=df*[0:N/2 -N/2+1:-1]; j = 1:N/2+1; % which freqs to plot
Y = freqfilter(X,fs,300,5000);
Xh = fft(X); Yh = fft(Y);
figure; subplot(2,1,1); t = (1:N)/fs; plot(t,X,'+-',t,Y,'o-');
axis([T/2-.002 T/2+.002 -.5 1]); title('impulse response');
subplot(2,1,2); plot(f(j),abs(Xh(:,j)).^2,'-');
hold on; plot(f(j),abs(Yh(:,j)).^2,'r-');
title('unfiltered and filtered |Fhat|');

% test applying to actual data:
%d = loaddata; d.A = freqfilter(d.A,d.samplefreq,300,10000); viewraw(d);
