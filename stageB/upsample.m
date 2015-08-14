function [Y t] = upsample(X,fac,kerpars)
% UPSAMPLE - interpolate rows of matrix onto possibly finer time grid
%
% Y = upsample(X,fac,kerpars) interpolates rows of matrix X into rows of Y,
%  using a time grid a factor fac finer. The output time grid is contained
%  within the input grid, allowing a buffer region of kerpars.Tf-1 at either
%  end. kerpars controls the kernel (filter) parameters.
%
% [Y t] = upsample(X,fac,kerpars) also gives the output time grid t.
%
% Inputs:
%   X = signal row vector, matrix, or 3d array.
%   fac = factor to upsample by (1,2,...)
%   kerpars = kernel parameters struct with fields:
%             Tf = half-width of filter in units of sample timesteps
% Outputs:
%   Y = matrix with same 1st (and 3rd) dimensions as X, but size(Y,2)=numel(t)
%
% Given assumed input time grid of 1,..,T, the output time grid t is
%  Tf, Tf+1/fac, Tf+2/fac, ..., T-Tf+1-1/fac, T-Tf+1.

% Barnett 11/26/14. default kerpars 6/10/15. faster 7/20/15

if nargin<1, test_upsample; return; end
if nargin<2 || isempty(fac), fac = 5; end   % defaults
if nargin<3, kerpars.Tf = 5; end

T = size(X,2);                              % # input time samples
t = (fac*kerpars.Tf:fac*(T-kerpars.Tf+1))/fac;  % output time grid
dt = bsxfun(@minus, t, (1:T)');             % matrix of time differences
PT = kernel(dt, kerpars);                   % transpose interp matrix
if ndims(X)==2                              % row vec or regular matrix
  Y = X*PT;                                 % right-mult to interp rows of X
else                                        % 3D case; loop it seems fast...
  Ns = size(X,3);
  Y = nan(size(X,1),numel(t),Ns);
  if 0  % loop over 3rd dim (slow if Ns huge)
    for s=1:Ns
      Y(:,:,s) = X(:,:,s)*PT;                 % right-mult to interp rows of X
    end
  else  % faster: avoid any loop
    M = size(X,1); X = reshape(permute(X,[1 3 2]), [M*Ns T]);
    Y = X*PT;                         % single matvec to interp all rows of X
    Y = ipermute(reshape(Y, [M Ns numel(t)]),[1 3 2]);
  end
end
%%%%%%%%%%%%%%

function f = kernel(t,pars)
% KERNEL - evaluate interpolation kernel at set of ordinates t-t'
%
% f = kernel(t,pars) where t is an array of time displacements and pars
% a parameter struct with fields:
%   pars.Tf  - half-width of support of kernel in original (input) samples
%
% Hann-windowed sinc for now
f = sin(pi*t)./(pi*t) .* cos((pi/2/pars.Tf)*t).^2; % compute even outside win
f(abs(t)>=pars.Tf) = 0;  % enforce compact support
f(t==0) = 1.0;           % fix nans in sinc at origin
%%%%%%%%%%%%%%

function test_kernel
pars.Tf = 7;
t = -10:0.1:10;
figure; plot(t,kernel(t,pars),'-'); hold on; t = -10:10; plot(t,t==0,'.');
title('interpolation kernel for Tf=7');


function test_upsample % three types of test
% i)
test_kernel

if 1 % ii) accuracy for sine wave at close to Nyquist
T = 40;   % # original samples
fac = 5;  % upsampling factor
Tf = 6; kerpars.Tf = Tf;   % kernel support half-width
t = 1:T;
samplerate = 1e4;
f = 4e3; g = @(t) sin(2*pi*f*t/samplerate); % test sine wave (freq < Nyquist)
u = g(t);  % row vector
[us ts] = upsample(u,fac,kerpars);
figure; plot(ts,us,'g.-',t,u,'k*'); xlabel('t'); title('test upsample');
Ns = numel(ts);
fprintf('rms error on output pts = %.3g\n',norm(us-g(ts))/sqrt(Ns))
% should be merely algebraic decay eg 1/Tf^4 or something
end

% iii) speed for 3d array case
clear; d = loaddemodata;
d.A = freqfilter(d.A,d.samplefreq,300,[]); [M N] = size(d.A);
[s t m info] = detectevents(d); Nt = size(s,2); Ns = size(s,3); % run it
fac = 5; kerpars.Tf = 5;
tic; [y tu] = upsample(s,fac,kerpars); toc
figure; subplot(1,2,1); imagesc(reshape(permute(s,[2 1 3]),[M*Nt Ns]));
colorbar; c=caxis; title('raw single events'); subplot(1,2,2);
imagesc(reshape(permute(y,[2 1 3]),[M*numel(tu) Ns])); caxis(c); colorbar;
title('upsampled single events');

