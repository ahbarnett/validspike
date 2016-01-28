function d = channelprewhiten(d,thresh,o)
% channelprewhiten - spatially pre-whiten a timeseries using C from noise clips
%
% function d = channelprewhiten(d,thresh,o)
%  d is a EC dataset object containing d.A timeseries data, d.samplefreq
%  On output d.A is changed. thresh sets the threshold below which a clip is
%  considered noise; if thresh=[] then a noise level is chosen automatically.
%  o - controls options such as:
%       o.verb - 0,1... verbosity
%       o.meth - 'c' Q=inv(chol(S))' Alex original (S = estimated covar mat)
%                'u' Q = U*inv(D)*U' where S = U*D^2*U' is eigendecomp of S,
%                    (Jeremy, 1/27/16)
%       o.rownorm = 1 row-normalizes Q (keeps overall variance same), 0 doesn't

% todo: self-test
%
% Barnett 6/4/15.
% Jeremy's Q and my variant (see: devel/compare_prewhite.m) 1/27/16

if nargin<3, o=[]; end
if ~isfield(o,'verb'), o.verb = 0; end
if ~isfield(o,'meth'), o.meth = 'c'; end
if ~isfield(o,'rownorm'), o.rownorm = 1; end

if isinf(thresh), disp('corr est via full sig')
  [M N] = size(d.A);
  S = (d.A*d.A') / N;                     % full signal covar est (faster)
else, disp('corr est via noise clips only')
  C = empiricalspacetimecorr(d,thresh,o);      % estimate noise auto-covar
  if o.verb, set(gcf,'name','before Q'); end
  S = squeeze(C(:,:,(size(C,3)+1)/2));      % zero-time-shift covar estimate
end

[U D] = eig(S); D = sqrt(D); % get left sing vectors, sing vals, of data mat
if strcmp(o.meth,'c')   % unmixing matrix choice... 
  Q = inv(chol(S)');
elseif strcmp(o.meth,'u')
  Q = U*inv(D)*U';
elseif strcmp(o.meth,'s')    % alex's variant of jeremy
  Q = inv(D)*U';
end
if o.rownorm
  for i=1:size(Q,1), Q(i,:) = Q(i,:)/norm(Q(i,:)); end
end      % weirdly, for meth='s', this produces very differing channel variances
if o.verb, figure; imagesc(Q); axis equal; title('Q channel prewhitening mat');
  caxis([-1 1]*max(abs(caxis))); colorbar; drawnow; end
d.A=Q*d.A;     % do the linear data transform
if o.verb>1   % check outgoing x-t corr
  C = empiricalspacetimecorr(d,thresh,o); set(gcf,'name','noise x-t corr, after Q');
end
