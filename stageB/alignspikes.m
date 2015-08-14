function [Y t] = alignspikes(X,f,fit,ta)
% ALIGNSPIKES  Time-align upsampled clips (eg by negative peak).
%
% [Y t] = alignspikes(X,f,fit,t)
%  produces a 3D array Y of shifted signals, aligned to their nearest subsample
%  (sampled a factor f more frequently than original sample rate), and the
%  time shifts t that were needed. The output time grid (size No) is shrunk in
%  number relative to the input, to pad a little and make it odd-numbered.
%  The peak is (unambiguously) aligned to the center element of the output grid.
%  
% Inputs:
%  X = M (# channels) by NT (# time subsamples) by Ns (# spikes) signal array
%  f = upsampling factor of the signal X from the original sample rate
%  fit = struct for L2-fitting alignment to fit.W with identities fit.L
%     (optional; obsolete). Leave as [] to skip such L2 fitting.
%  t = specified alignment times, prevents fitting of t (optional)
%
% Outputs:
%  Y = M (# channels) by No (# output time subsamples) by Ns (# spikes) signal
%  t = 1 by Ns list of real time shifts, in units of original samples
%
% Note that No is chosen by alignspikes

% todo:  change detectevents so it uses same criterion for a peak as here?
%    (this would guarantee shifts to be within one original sample)
% todo: fix L2 fit realignment so it's in data space not upsampled space
%         (speeds up, low priority).
% todo: maxfitslide an opt?
%
% Barnett 11/26/14. L2-fit 1/22/15. input t override 1/30/15

if nargin<1, test_alignspikes; return; end
fitalign = nargin>2 && ~isempty(fit);
givent = nargin>3;

[M,Nt,Ns] = size(X);
No = Nt - ceil(2.0 * f);  % # subsamples to output, allowing some sliding. opt?
if mod(No,2)==0, No=No-1; end   % make odd so has a center entry
maxfitslide = 10.0; maxsh = ceil(f*maxfitslide);  % max dt fit in orig samples
indsh = [-maxsh:maxsh];         % set of index shifts to explore in fit
Y = zeros(M,No,Ns); t = nan(1,Ns);  % allocate (zero since not all copied over)
for s=1:Ns          % loop over single-spike events
  A = X(:,:,s);     % current spike data (M * Nt)
  if givent         % use t to set index jo, don't fit
    jo = round(f*ta(s) + 1);
  elseif ~fitalign      % use plain across-channel min voltage
    j = find(A==min(A(:))); [~,jo] = ind2sub(size(A),j); % t index in A of min
  else              % L2 fit to existing waveforms with labels
    lj = fit.L(s); if lj==0, lj=1; end % hack so unlabeled pts aligned to sth!
    w = fit.W(:,:,lj); % appropriate waveform for this event, M-by-Nt
    errs = 0*indsh;  % will store mean squared L2 errs
    for a=1:numel(indsh), d = indsh(a);  % loop over index shifts
      i = 1+max(0,-d):No+min(0,-d);       % indices in w array to get err at
      errs(a) = sum(sum((w(:,i) - A(:,i+d)).^2)) / numel(i);
    end
    [~,abest] = min(errs);
    jo = indsh(abest) + (No+1)/2;    % central index
    %fprintf('%d\n',indsh(abest))
    % todo: ampl fit now with this best index fixed
  end
  t(s) = (jo-1)/f;   % peak time relative to first upsampled input
  j = max(1,jo-(No-1)/2):min(Nt,jo+(No-1)/2); % source index list (max len No)
  i = j - (jo-(No-1)/2-1);                    % targ index list, 2nd dim of Y
  if givent, ii = find(i>0 & i<=No); i=i(ii); j=j(ii); end  % safety clip
  Y(:,i,s) = X(:,j,s);                        % copy in
end

%%%%%%%%%%%%%

function test_alignspikes      % demo alignment on real data; visual test only
fitalign = 0; % choose if test L2 fitting or just plain peak
d = loaddemodata;
d.A = freqfilter(d.A,d.samplefreq,300,10000); [M N] = size(d.A); % Stage A
[s t m info] = detectevents(d); Nt = size(s,2); Ns = size(s,3);
fac = 5; kerpars.Tf = 5;  % upsampling
if ~fitalign
  fprintf('upsampled peak real data alignment test...\n'), fit = [];
else
  fprintf('upsampled L2 fit real data alignment test...\n')
  o.K = 4; % zero-th pass Stage B just to get fit waveforms & labels...
  o.upsampfac = fac; [fit.L fit.W] = ssalg_feature_cluster(s,o);
end
[y tu] = upsample(s,fac,kerpars);
tic; [Y ta] = alignspikes(y,fac,fit); toc; No = size(Y,2);
t = t - 1 + kerpars.Tf + ta; % firing times (t was coarse window start indices)

figure; plot(1:N,d.A','-'); axis([2000 4500 -250 100]); % first 6 spikes in 'h'
hold on; plot(t,0*t,'k*');  % check alignment of supposed minima
for i=1:6 % overplot some upsampled signals to check matches data...
  plot(t(i)+(-(No-1)/2:(No-1)/2)/fac, Y(:,:,i)', '.'); % assumes No odd
end
ta(1:6) % w/ buszaki data spike 5 is v dependent on fitalign: 2 nearby spikes

figure; subplot(1,2,1); imagesc(reshape(permute(y,[2 1 3]),[M*numel(tu) Ns]));
colorbar; c=caxis; title('upsampled single spikes'); subplot(1,2,2);
imagesc(reshape(permute(Y,[2 1 3]),[M*size(Y,2) Ns])); caxis(c); colorbar;
title('aligned single spikes');
